#!/usr/bin/env python3
"""
Merge the outputs of one group of pocket-coffea jobs (one dataset, see pbs_tools.merge_groups) into
<outputdir>/merged/<group>.coffea. Run by the per-dataset merge array of pbs_tools.py inside the
apptainer image (needs pocket-coffea's merge-outputs and PyYAML).

    python merge_group.py <jobs_dir> <groups file> <line index> [-m MAX_MEM_GB]

The groups file has one line per group: '<name> <job id>,<job id>,...'. A jobs_config restricted to
the group's jobs is written to <jobs_dir>/merge/jobs_config_<name>.yaml and passed to
`merge-outputs -jc`, so the usual checks apply to the group: every job's output must exist (jobs
with a missing output are marked .failed), then the outputs are merged and postprocessed. Each
group has its own merge cache (<outputdir>/merge_cache/<name>), so the groups can run in parallel.
"""
import argparse
import os
import subprocess
import sys

import yaml

MERGED_DIR = "merged"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("jobs_dir")
    ap.add_argument("groups_file")
    ap.add_argument("index", type=int, help="0-based line of the groups file (PBS_ARRAY_INDEX)")
    ap.add_argument("-m", "--max-mem-gb", type=float, default=None)
    args = ap.parse_args()

    with open(args.groups_file) as f:
        lines = [l.split() for l in f if l.strip()]
    if not 0 <= args.index < len(lines):
        sys.exit(f"index {args.index} outside the {len(lines)} groups of {args.groups_file}")
    name, ids = lines[args.index][0], [int(i) for i in lines[args.index][1].split(",")]

    with open(os.path.join(args.jobs_dir, "jobs_config.yaml")) as f:
        cfg = yaml.safe_load(f)
    jobs = {f"job_{i}": cfg["jobs_list"][f"job_{i}"] for i in ids}
    cfg["jobs_list"] = jobs
    os.makedirs(os.path.join(args.jobs_dir, "merge"), exist_ok=True)
    sub_cfg = os.path.join(args.jobs_dir, "merge", f"jobs_config_{name}.yaml")
    with open(sub_cfg, "w") as f:
        yaml.dump(cfg, f)

    out_dir = os.path.join(cfg["output_dir"], MERGED_DIR)
    os.makedirs(out_dir, exist_ok=True)
    cmd = ["merge-outputs", "-jc", sub_cfg, "-o", os.path.join(out_dir, f"{name}.coffea"), "--mark-failed", "-f",
           "-cd", os.path.join(cfg["output_dir"], "merge_cache", name)]
    if args.max_mem_gb:
        cmd += ["-m", f"{args.max_mem_gb:.0f}"]
    datasets = sorted({d for j in jobs.values() for d in (j.get("filesets") or {})})
    print(f"Group {name}: jobs {ids}, datasets {datasets}")
    print(" ".join(cmd), flush=True)
    sys.exit(subprocess.run(cmd).returncode)


if __name__ == "__main__":
    main()
