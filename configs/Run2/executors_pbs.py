"""
pocket-coffea executor for PBS Pro: one batch job per chunk of files, each running the futures
executor on all cores of its node inside apptainer, plus a merge job.

    pocket-coffea run --cfg inference_conf_signal.py -o outputs/nominal/signal \
        --executor pbs --executor-custom-setup executors_pbs.py -ro run_options_pbs.yaml

Splitting, the per-job configurator pickles, jobs_config.yaml and the job.sh wrapper are
pocket-coffea's manual-jobs machinery (same as the HTCondor executor); only the submission is
PBS-specific (pbs_tools.py). Inside the apptainer image qsub is not available: the jobs are then
only prepared, and submitted from the host with `python pbs_tools.py submit <outputdir>`
(run_pipeline.py --backend pbs does both steps).
"""
import json
import os
import shutil
import sys

import cloudpickle
import yaml

from pocket_coffea.executors.executors_base import IterativeExecutorFactory, FuturesExecutorFactory
from pocket_coffea.executors.executors_lxplus import ExecutorFactoryCondorCERN, build_job_script, get_worker_env
from pocket_coffea.executors.executors_manual_jobs import write_inner_run_options, INNER_RUN_OPTIONS_WHITELIST

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pbs_tools  # noqa: E402


def restrict_to_split(config, split):
    """Drop the `bysample` entries (columns, weights, variations) of samples that are not in this
    job. The per-job configurator only knows the samples of its own split, and
    Configurator.load() rejects entries for other samples ("Wrong columns configuration").
    Works on copies: the cloned configurators share these dicts with the full configuration."""
    samples = {ds["metadata"]["sample"] for ds in split.values()}

    def keep(key):
        return key in samples or key.split("__")[0] in samples

    def pruned(cfg):
        if not cfg or "bysample" not in cfg:
            return cfg
        cfg = dict(cfg)
        cfg["bysample"] = {k: v for k, v in cfg["bysample"].items() if keep(k)}
        return cfg

    config.columns_cfg = pruned(config.columns_cfg)
    config.weights_cfg = pruned(config.weights_cfg)
    if config.variations_cfg:
        config.variations_cfg = {kind: pruned(cfg) for kind, cfg in config.variations_cfg.items()}


class ExecutorFactoryPBS(ExecutorFactoryCondorCERN):
    """Reuses the HTCondor executor's splitting (prepare_splitting); replaces the submission with
    a PBS Pro job array."""

    def prepare_jobs(self, splits):
        """Same as the HTCondor executor (one pickled configurator per job + jobs_config.yaml),
        with each job's configurator restricted to the samples of its split."""
        jobs_config = {
            "job_name": self.job_name,
            "job_dir": os.path.abspath(self.jobs_dir),
            "output_dir": os.path.abspath(self.outputdir),
            "split_by_category": self.run_options["split-by-category"],
            "config_pkl_total": f"{os.path.abspath(self.outputdir)}/configurator.pkl",
            "jobs_list": {},
        }
        # Postprocessing happens once, in the merge job
        self.config.do_postprocessing = False
        config_files = []
        for i, split in enumerate(splits):
            partial_config = self.config.clone()
            restrict_to_split(partial_config, split)
            partial_config.set_filesets_manually(split)
            path = f"{self.jobs_dir}/config_job_{i}.pkl"
            with open(path, "wb") as f:
                cloudpickle.dump(partial_config, f)
            config_files.append(path)
            jobs_config["jobs_list"][f"job_{i}"] = {
                "filesets": split,
                "config_file": path,
                "output_file": f"{os.path.abspath(self.outputdir)}/output_job_{i}.coffea",
            }
        with open(f"{self.jobs_dir}/jobs_config.yaml", "w") as f:
            yaml.dump(jobs_config, f)
        return config_files

    def option(self, key):
        value = self.run_options.get(key, None)
        return pbs_tools.DEFAULTS[key] if value is None and key in pbs_tools.DEFAULTS else value

    def submit_jobs(self, jobs_config):
        jobs_dir = os.path.abspath(self.jobs_dir)
        output_dir = os.path.abspath(self.outputdir)
        analysis_dir = os.path.abspath(self.run_options.get("analysis-dir") or os.getcwd())
        os.makedirs(os.path.join(jobs_dir, "logs"), exist_ok=True)

        # Worker environment: the analysis modules (cand_helper, genmatcher, ...) must be importable
        env = get_worker_env(self.run_options, getattr(self, "x509_path", None), "condor")
        env.append(f"export PYTHONPATH={analysis_dir}:$PYTHONPATH")

        # The per-job `pocket-coffea run` must also know whether a grid proxy is required
        inner_yaml = write_inner_run_options(self.jobs_dir, self.run_options,
                                             whitelist=INNER_RUN_OPTIONS_WHITELIST + ("ignore-grid-certificate",))
        script = build_job_script(
            env_extras="\n".join(env),
            abs_jobdir_path=jobs_dir,
            abs_output_path=output_dir,
            copy_command="cp",
            runnercmd="pocket-coffea run",
            inner_yaml_basename=os.path.basename(inner_yaml),
            split_by_category=self.run_options["split-by-category"],
            cores_per_worker=self.option("cores-per-worker"),
        )
        with open(os.path.join(jobs_dir, "job.sh"), "w") as f:
            f.write(script)
        os.chmod(os.path.join(jobs_dir, "job.sh"), 0o755)

        chunksize_cfg = self.run_options["chunksize"]
        self._validate_chunksize_keys(chunksize_cfg, self.filesets)
        chunksizes = [int(self._resolve_chunksize_for_job(chunksize_cfg, split)) for split in self._splits]

        name = self.run_options.get("pbs-name") or "_".join(output_dir.rstrip("/").split("/")[-2:])
        settings = {key: self.option(key) for key in pbs_tools.DEFAULTS}
        settings["apptainer-binds"] = list(settings["apptainer-binds"])
        settings.update({
            "name": "".join(c if c.isalnum() or c == "_" else "_" for c in name)[:15],
            "analysis-dir": analysis_dir,
            "output-dir": output_dir,
            "inner-yaml": os.path.basename(inner_yaml),
            # proxy used by the jobs (refreshed by pbs_tools at every submission); None = no proxy
            "x509-path": None if self.run_options["ignore-grid-certificate"] else self.x509_path,
            "chunksizes": chunksizes,
        })
        with open(os.path.join(jobs_dir, pbs_tools.SETTINGS_FILE), "w") as f:
            json.dump(settings, f, indent=1)
        print(f"Prepared {len(chunksizes)} PBS job(s) in {jobs_dir}")

        if self.run_options.get("dry-run", False) or os.environ.get(pbs_tools.NO_SUBMIT_ENV):
            print(f"Dry run: submit with  python {pbs_tools.__file__} submit {output_dir}")
        elif shutil.which("qsub") is None:
            print("qsub is not available here (inside the container?). Submit from the host with\n"
                  f"    python {pbs_tools.__file__} submit {output_dir}")
        else:
            pbs_tools.submit(output_dir)

    def recreate_jobs(self, jobs_to_recreate):
        sys.exit(f"For PBS use:  python {pbs_tools.__file__} resubmit {os.path.abspath(self.outputdir)}")


def get_executor_factory(executor_name, **kwargs):
    if executor_name == "iterative":
        return IterativeExecutorFactory(**kwargs)
    elif executor_name == "futures":
        return FuturesExecutorFactory(**kwargs)
    elif executor_name == "pbs":
        return ExecutorFactoryPBS(**kwargs)
    raise ValueError(f"Executor '{executor_name}' not available in executors_pbs.py (iterative, futures, pbs)")
