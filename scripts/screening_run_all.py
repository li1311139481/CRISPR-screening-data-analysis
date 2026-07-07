#!/usr/bin/env python3
"""Generate and optionally submit pairwise MAGeCK comparison commands."""

import argparse
import os
import shlex
import subprocess as sb
import sys
from itertools import permutations

import pandas as pd


def q(value: str) -> str:
    return shlex.quote(str(value))


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate MAGeCK pairwise comparison commands")
    parser.add_argument("-d", "--DIR", required=True, help="Project working directory")
    parser.add_argument("-l", "--library", required=True, help="sgRNA library annotation file")
    parser.add_argument("-c", "--controlid", required=True, help="Control sgRNA ID file")
    parser.add_argument("-n", "--prefix", required=True, help="Output prefix, e.g. TF")
    parser.add_argument("-p", "--parallel", type=int, required=True, help="Number of parallel MAGeCK jobs")
    parser.add_argument("-m", "--metadata", default=None, help="Path to meta.data.csv (default: DIR/meta.data.csv)")
    parser.add_argument("--submit", action="store_true", help="Submit run.slurm automatically")
    parser.add_argument("--nodelist", default=None, help="Slurm nodelist (e.g. node2), omit if not needed")
    parser.add_argument("--python", default=sys.executable, help="Python executable")
    parser.add_argument("--mageck-wrapper", default=None, help="Path to mageck.py wrapper")
    args = parser.parse_args()

    DIR = os.path.abspath(args.DIR)
    library = os.path.abspath(args.library)
    controlid = os.path.abspath(args.controlid)
    prefix = args.prefix
    parallel = args.parallel

    script_dir = os.path.dirname(os.path.abspath(__file__))
    mageck_wrapper = args.mageck_wrapper or os.path.join(script_dir, "mageck.py")

    print(f"library is: {library}")
    print(f"control_id is: {controlid}")

    metadata_path = args.metadata or os.path.join(DIR, "meta.data.csv")
    print(f"metadata is: {metadata_path}")
    meta_data = pd.read_csv(metadata_path)

    required_cols = {"sample_name", "sample_number", "group_name"}
    missing = required_cols - set(meta_data.columns)
    if missing:
        raise ValueError(f"{metadata_path} is missing required columns: {sorted(missing)}")

    group_list = meta_data["group_name"].dropna().unique()
    print("group list is:")
    print(list(group_list))

    group_permutations = list(permutations(group_list, 2))
    cutadapt_path = os.path.join(DIR, "cutadapt")
    command_file = os.path.join(DIR, "commands.sh")

    if os.path.exists(command_file):
        os.remove(command_file)

    generated = 0
    for group_control, group_treatment in group_permutations:
        control_samples = meta_data.loc[meta_data["group_name"] == group_control, "sample_name"].tolist()
        treatment_samples = meta_data.loc[meta_data["group_name"] == group_treatment, "sample_name"].tolist()

        name_control = ",".join(control_samples)
        name_treatment = ",".join(treatment_samples)
        comparison_name_for_labels = f"{name_control}--{name_treatment}"

        control_fastqs = meta_data.loc[meta_data["group_name"] == group_control, "sample_number"].apply(
            lambda x: os.path.join(cutadapt_path, f"{x}_1.trm.fq.gz")
        ).tolist()
        treatment_fastqs = meta_data.loc[meta_data["group_name"] == group_treatment, "sample_number"].apply(
            lambda x: os.path.join(cutadapt_path, f"{x}_1.trm.fq.gz")
        ).tolist()
        all_fastqs = ",".join(control_fastqs + treatment_fastqs)

        output_dir = os.path.join(DIR, f"{prefix}_{group_treatment}.vs.{group_control}")
        volcano_file = os.path.join(output_dir, "VolcanoView.pdf")

        if not os.path.exists(volcano_file):
            cmd = " ".join([
                q(args.python),
                q(mageck_wrapper),
                "-p", q(prefix),
                "-l", q(library),
                "-c", q(controlid),
                "-d", q(DIR),
                "-g", q(f"{group_control},{group_treatment}"),
                "-n", q(comparison_name_for_labels),
                "-i", q(all_fastqs),
                ">", q(os.path.join(DIR, f"{group_treatment}.vs.{group_control}.log")), "2>&1",
            ])
            with open(command_file, "a", encoding="utf-8") as cmd_file:
                cmd_file.write(cmd + "\n")
            generated += 1
            print(f"Generated command for {group_treatment} vs {group_control}")

    if generated == 0:
        print(f"\n所有 {len(group_permutations)} 对比较的 VolcanoView.pdf 已存在，无需生成新命令。")
        print(f"如需重新运行，请删除对应的输出目录后重试。")
        return

    nodelist_line = f"#SBATCH --nodelist={args.nodelist}\n" if args.nodelist else ""
    slurm_content = f"""#!/bin/bash
#SBATCH -J {prefix}
#SBATCH -o {DIR}/%j.%x.out
#SBATCH -e {DIR}/%j.%x.err
#SBATCH --partition=cpu
#SBATCH --cpus-per-task={parallel}
{nodelist_line}
set -euo pipefail

parallel --jobs {parallel} < {command_file}
"""

    slurm_file = os.path.join(DIR, "run.slurm")
    with open(slurm_file, "w", encoding="utf-8") as f:
        f.write(slurm_content)

    os.chmod(command_file, 0o750)
    os.chmod(slurm_file, 0o750)

    print(f"\nGenerated {generated} commands")
    print(f"command script generated at:\n{command_file}\n")
    print(f"Slurm job script generated at:\n{slurm_file}\n")

    if args.submit:
        print("Submitting Slurm job")
        sb.run(["sbatch", slurm_file], check=True)
    else:
        print("Slurm script generated but not submitted. Add --submit to submit automatically.")


if __name__ == "__main__":
    main()
