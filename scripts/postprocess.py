#!/usr/bin/env python3
"""Post-processing after MAGeCK count/test: merge count tables, run per-sample
library QC, and run MAGeCKFlute report for every comparison.

Designed to be appended to the end of run.slurm (after parallel MAGeCK jobs finish),
so all count/test outputs are guaranteed to exist.
"""

import argparse
import glob
import os
import subprocess as sb

import pandas as pd


def q(value: str) -> str:
    return "'" + str(value).replace("'", "'\\''") + "'"


def find_count_tables(DIR: str, prefix: str):
    """Return all {prefix}_*.vs.*/*_count.count.txt under DIR."""
    pattern = os.path.join(DIR, f"{prefix}_*.vs.*", "*_count.count.txt")
    return sorted(glob.glob(pattern))


def merge_counts(count_tables: list, out_path: str) -> None:
    """Merge multiple count tables (sgRNA, Gene, sample columns) by sgRNA into one."""
    # Each MAGeCK count table contains the full library (same sgRNA set),
    # but only the samples used in that comparison. Merge by outer join on sgRNA.
    merged = None
    for table in count_tables:
        df = pd.read_csv(table, sep="\t")
        if "Gene" not in df.columns:
            df["Gene"] = ""
        if merged is None:
            # First table: keep sgRNA + Gene, then its sample columns.
            sample_cols = [c for c in df.columns if c not in ("sgRNA", "Gene")]
            merged = df.copy()
            all_sample_cols = list(sample_cols)
        else:
            # Gene annotation: back-fill for any sgRNA not yet present.
            gene_map = df.set_index("sgRNA")["Gene"].to_dict()
            existing = set(merged["sgRNA"])
            new_gene_rows = [{"sgRNA": g, "Gene": a}
                             for g, a in gene_map.items() if g not in existing and str(a) != "nan"]
            if new_gene_rows:
                merged = pd.concat([merged, pd.DataFrame(new_gene_rows)], ignore_index=True)
            # Add new sample columns via outer merge on sgRNA.
            new_cols = [c for c in df.columns if c not in ("sgRNA", "Gene") and c not in merged.columns]
            all_sample_cols.extend(new_cols)
            if new_cols:
                merged = merged.merge(df[["sgRNA"] + new_cols], on="sgRNA", how="outer")
    # Clean up.
    merged = merged.drop_duplicates(subset="sgRNA", keep="first")
    merged = merged.fillna({**{c: 0 for c in all_sample_cols}})
    merged = merged[["sgRNA", "Gene"] + all_sample_cols].sort_values("sgRNA")
    merged.to_csv(out_path, sep="\t", index=False)
    print(f"[postprocess] Merged {len(count_tables)} count tables -> {out_path}")


def run_qc(DIR: str, count_table: str, sample_cols: list, py: str, qc_script: str) -> None:
    """Run plot_library_qc.py --all on the merged table (regenerates/overwrites PDFs)."""
    cmd = [py, qc_script, "-i", count_table, "--all", "-o", DIR]
    print("[postprocess] Running library QC:")
    print("  " + " ".join(cmd))
    sb.run(cmd, check=True)


def run_flute_all(DIR: str, prefix: str) -> None:
    """bash every comparison dir's run_flute.sh (skip if missing)."""
    comp_dirs = sorted(glob.glob(os.path.join(DIR, f"{prefix}_*.vs.*")))
    for comp_dir in comp_dirs:
        flute_script = os.path.join(comp_dir, "run_flute.sh")
        if not os.path.exists(flute_script):
            continue
        print(f"[postprocess] Running MAGeCKFlute report in {comp_dir}")
        sb.run(["bash", flute_script], check=True)


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Post-process MAGeCK outputs: merge count tables by default; "
                    "optionally also run library QC (--qc) and/or MAGeCKFlute reports (--flute).")
    ap.add_argument("-d", "--DIR", required=True, help="Project working directory")
    ap.add_argument("-n", "--prefix", default="TF", help="Output folder prefix, e.g. TF")
    ap.add_argument("--python", default="python",
                    help="Python executable to run plot_library_qc.py (default 'python', resolved from PATH)")
    ap.add_argument("-s", "--scripts", default=None, help="Path to scripts dir (default: alongside this file)")
    ap.add_argument("--qc", action="store_true",
                    help="Also run per-sample library QC (plot_library_qc.py --all) after merging")
    ap.add_argument("--flute", action="store_true",
                    help="Also run MAGeCKFlute report for every comparison after merging")
    ap.add_argument("--flute-only", action="store_true",
                    help="Only run MAGeCKFlute reports for every comparison (skip merging/QC)")
    args = ap.parse_args()

    DIR = os.path.abspath(args.DIR)
    script_dir = os.path.abspath(args.scripts or os.path.dirname(os.path.abspath(__file__)))
    qc_script = os.path.join(script_dir, "plot_library_qc.py")

    count_tables = find_count_tables(DIR, args.prefix)

    # --flute-only: regenerate every comparison's MAGeCKFlute report only.
    if args.flute_only:
        if not count_tables:
            print(f"[postprocess] No count tables found under {DIR}/{args.prefix}_*.vs.* — nothing to do.")
            return
        run_flute_all(DIR, args.prefix)
        return

    if not count_tables:
        print(f"[postprocess] No count tables found under {DIR}/{args.prefix}_*.vs.* — nothing to do.")
        return

    # Default behaviour: merge per-comparison count tables into all_samples.count.txt.
    all_counts = os.path.join(DIR, "all_samples.count.txt")
    merge_counts(count_tables, all_counts)

    # Optional per-sample library QC (regenerates library_qc_<sample>.pdf at project root).
    if args.qc:
        run_qc(DIR, all_counts, [], args.python, qc_script)

    # Optional MAGeCKFlute report for every comparison.
    if args.flute:
        run_flute_all(DIR, args.prefix)


if __name__ == "__main__":
    main()
