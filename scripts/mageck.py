#!/usr/bin/env python3
"""Run MAGeCK count/test and optional MAGeCKFlute reporting for one comparison."""

import argparse
import os
import subprocess as sb
from typing import List


def run_command(cmd: List[str]) -> None:
    print("[INFO] Running command:")
    print(" ".join(cmd))
    sb.run(cmd, check=True)


def main() -> None:
    parser = argparse.ArgumentParser(description="Run MAGeCK analysis for one comparison")
    parser.add_argument("-g", "--group_name", required=True, help="Control and treatment group names, e.g. Input,Day7_Total")
    parser.add_argument("-n", "--sample_name", required=True, help="Control and treatment sample labels, e.g. Input_1,Input_2--Day7_Total_1,Day7_Total_2")
    parser.add_argument("-i", "--sample_file", required=True, help="Comma-separated trimmed FASTQ files")
    parser.add_argument("-d", "--DIR", default=None, help="Project working directory")
    parser.add_argument("-p", "--prefix", default="mageck", help="Output folder prefix")
    parser.add_argument("-l", "--library", required=True, help="sgRNA library annotation file")
    parser.add_argument("-c", "--controlid", required=True, help="Control sgRNA ID file")
    parser.add_argument("--mageck", default="mageck", help="MAGeCK executable")
    parser.add_argument("--rscript", default="Rscript", help="Rscript executable")
    parser.add_argument("--skip-flute", action="store_true", help="Skip MAGeCKFlute reporting")
    args = parser.parse_args()

    DIR = os.path.abspath(args.DIR or os.getcwd())
    script_dir = os.path.dirname(os.path.abspath(__file__))
    mageck_flute_script = os.path.join(script_dir, "mageck_flute.R")

    group_name = args.group_name.split(",")
    sample_name = args.sample_name.split("--")
    sample_files = [x for x in args.sample_file.split(",") if x]

    if len(group_name) != 2:
        raise ValueError("--group_name must contain exactly two comma-separated groups")
    if len(sample_name) != 2:
        raise ValueError("--sample_name must contain control and treatment labels separated by --")

    control_group, treatment_group = group_name
    control_samples, treatment_samples = sample_name
    sample_labels = f"{control_samples},{treatment_samples}"

    print("group_name sep = ,   sample_name sep = --   sample_file sep = ,")
    print(f"[INFO] group names: {group_name}")
    print(f"[INFO] sample names: {sample_name}")
    print("[INFO] sample files:")
    print("\n".join(sample_files))
    print(f"[INFO] library file: {args.library}")
    print(f"[INFO] control sgRNA file: {args.controlid}")
    print(f"[INFO] output folder root: {DIR}")

    mageck_folder_name = os.path.join(DIR, f"{args.prefix}_{treatment_group}.vs.{control_group}")
    comparison_name = f"{treatment_group}.vs.{control_group}"
    os.makedirs(mageck_folder_name, exist_ok=True)

    count_prefix = os.path.join(mageck_folder_name, f"{comparison_name}_count")
    test_prefix = os.path.join(mageck_folder_name, f"{comparison_name}_test")
    count_table = f"{count_prefix}.count.txt"

    cmd_count = [
        args.mageck, "count",
        "-l", args.library,
        "-n", count_prefix,
        "--sample-label", sample_labels,
        "--fastq", *sample_files,
        "--control-sgrna", args.controlid,
    ]
    run_command(cmd_count)

    cmd_test = [
        args.mageck, "test",
        "-k", count_table,
        "-t", treatment_samples,
        "-c", control_samples,
        "-n", test_prefix,
        "--gene-lfc-method", "mean",
        "--control-sgrna", args.controlid,
    ]
    run_command(cmd_test)

    if not args.skip_flute:
        cmd_flute = [
            args.rscript,
            mageck_flute_script,
            mageck_folder_name,
            sample_labels,
            comparison_name,
        ]
        run_command(cmd_flute)


if __name__ == "__main__":
    main()
