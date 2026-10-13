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

from pocket_coffea.executors.executors_base import IterativeExecutorFactory, FuturesExecutorFactory
from pocket_coffea.executors.executors_lxplus import ExecutorFactoryCondorCERN, build_job_script, get_worker_env
from pocket_coffea.executors.executors_manual_jobs import write_inner_run_options

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pbs_tools  # noqa: E402


class ExecutorFactoryPBS(ExecutorFactoryCondorCERN):
    """Reuses the HTCondor executor's splitting (prepare_splitting) and per-job configurator
    pickles (prepare_jobs); replaces the submission with a PBS Pro job array."""

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

        inner_yaml = write_inner_run_options(self.jobs_dir, self.run_options)
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
