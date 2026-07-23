#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
CRISPR 文库质控三图绘制 (per-sample, 带样本名)
  Fig3  sgRNA Readcount 核密度图 (KDE)          -> 2.2.1
  Fig6  均一性斜率图 (Skew Ratio)               -> 2.3.2
  Fig7  sgRNA read counts 分布 AUC 图 (Lorenz)  -> 2.3.3

输入: MAGeCK count 格式的 readcount 表 (tab 分隔)
  列: sgRNA  Gene  <sample1>  <sample2>  ...
例:
  sgRNA    Gene    library_huang
  A1BG_1   A1BG    4321
  A1BG_2   A1BG    5689
  ...

用法:
  # 单个样本 (输出文件名与标题都带样本名)
  python plot_library_qc.py -i all.count.txt -s Day7_Total_1 -o ./qc_pdf

  # 批量: 自动遍历 readcount 表里所有样本列, 每个样本生成一个 PDF
  python plot_library_qc.py -i all.count.txt --all -o ./qc_pdf

  # 模拟数据快速预览效果
  python plot_library_qc.py --demo -o ./qc_pdf

依赖: numpy pandas scipy matplotlib
判定阈值: Skew Ratio < 10 合格 ; AUC < 0.7 合格
"""
import argparse
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from scipy.stats import gaussian_kde

# 兼容 numpy 1.x (trapz) / 2.0 (trapezoid)
_trapz = getattr(np, "trapezoid", None) or np.trapz


def load_counts(path, sample):
    df = pd.read_csv(path, sep="\t")
    cols = df.columns.tolist()
    if sample not in cols:
        raise ValueError(f"sample '{sample}' 不在文件列里: {cols}")
    return df[sample].astype(float).values


def list_samples(path):
    """返回 readcount 表中除 sgRNA/Gene 之外的所有样本列名"""
    df = pd.read_csv(path, sep="\t", nrows=1)
    skip = {"sgRNA", "Gene"}
    return [c for c in df.columns if c not in skip]


def demo_counts(n=1984, seed=42):
    """生成类质粒文库 readcount (Log2 中位 ~12, Skew ~5, AUC ~0.67)"""
    rng = np.random.default_rng(seed)
    return (rng.lognormal(mean=8.3, sigma=0.6, size=n).astype(int) + 10)


def fig3_kde(counts, sample):
    c = counts[counts > 0]
    log2c = np.log2(c)
    kde = gaussian_kde(log2c, bw_method="scott")
    xs = np.linspace(log2c.min(), log2c.max(), 500)
    ys = kde(xs)

    fig, ax = plt.subplots(figsize=(6, 4.5))
    ax.plot(xs, ys, color="steelblue", lw=1.8)
    ax.fill_between(xs, ys, alpha=0.30, color="steelblue")
    ax.set_xlabel("Log2(sgRNA Readcount)")
    ax.set_ylabel("Density")
    ax.set_title(f"sgRNA Readcount density (KDE)  |  {sample}")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    return fig


def fig6_skew(counts, sample):
    c = counts[counts > 0]
    sorted_c = np.sort(c)
    n = len(sorted_c)
    cum_frac = np.arange(1, n + 1) / n          # 经验 CDF: 累积到该 count 的 sgRNA 比例
    p10 = float(np.percentile(c, 10))
    p90 = float(np.percentile(c, 90))
    skew = p90 / p10

    fig, ax = plt.subplots(figsize=(6, 4.5))
    ax.plot(sorted_c, cum_frac, color="red", lw=1.6)
    ax.axhline(0.1, ls="--", color="grey", lw=0.8)
    ax.axhline(0.9, ls="--", color="grey", lw=0.8)
    ax.axvline(p10, ls=":", color="grey", lw=0.8)
    ax.axvline(p90, ls=":", color="grey", lw=0.8)
    ax.scatter([p10], [0.1], color="black", zorder=5, s=25)
    ax.scatter([p90], [0.9], color="black", zorder=5, s=25)
    ax.annotate(f"10%: {p10:.0f}", xy=(p10, 0.1), xytext=(p10 + max(c) * 0.02, 0.13), fontsize=9)
    ax.annotate(f"90%: {p90:.0f}", xy=(p90, 0.9), xytext=(p90 + max(c) * 0.02, 0.82), fontsize=9)
    ax.set_xlabel("Sequencing depth (reads per sgRNA)")
    ax.set_ylabel("Cumulative fraction of sgRNAs")
    status = "PASS" if skew < 10 else "FAIL"
    ax.set_title(f"Uniformity slope  |  Skew ratio = {skew:.2f}  |  {sample}")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    return fig, skew


def fig7_auc(counts, sample):
    c = counts[counts > 0]
    sorted_c = np.sort(c)[::-1]                 # 降序 -> 曲线凸向左上, AUC ∈ [0.5,1]
    cum_ab = np.concatenate([[0.0], np.cumsum(sorted_c) / sorted_c.sum()])
    x = np.concatenate([[0.0], np.arange(1, len(sorted_c) + 1) / len(sorted_c)])
    auc = float(_trapz(cum_ab, x))

    fig, ax = plt.subplots(figsize=(6, 4.5))
    ax.plot(x, cum_ab, color="red", lw=1.6, label="observed")
    ax.plot([0, 1], [0, 1], "k--", lw=1, label="ideal (uniform)")
    ax.fill_between(x, cum_ab, alpha=0.15, color="red")
    ax.set_xlabel("Cumulative fraction of sgRNAs (ranked by abundance)")
    ax.set_ylabel("Cumulative fraction of total reads")
    status = "PASS" if auc < 0.7 else "FAIL"
    ax.set_title(f"sgRNA read-count AUC  |  AUC = {auc:.4f}  |  {sample}")
    ax.legend(loc="upper left", frameon=False, fontsize=9)
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    return fig, auc


def make_pdf(counts, sample, outdir, keep_png=False):
    """生成一个样本的三图合并 PDF (文件名带样本名), 可选保留 PNG。返回 (pdf_path, skew, auc)"""
    os.makedirs(outdir, exist_ok=True)
    pdf_path = os.path.join(outdir, f"library_qc_{sample}.pdf")

    with PdfPages(pdf_path) as pdf:
        fig3 = fig3_kde(counts, sample)
        if keep_png:
            fig3.savefig(os.path.join(outdir, f"Fig3_kde_{sample}.png"), dpi=150)
        pdf.savefig(fig3, dpi=150)
        plt.close(fig3)

        fig6, skew = fig6_skew(counts, sample)
        if keep_png:
            fig6.savefig(os.path.join(outdir, f"Fig6_skew_{sample}.png"), dpi=150)
        pdf.savefig(fig6, dpi=150)
        plt.close(fig6)

        fig7, auc = fig7_auc(counts, sample)
        if keep_png:
            fig7.savefig(os.path.join(outdir, f"Fig7_auc_{sample}.png"), dpi=150)
        pdf.savefig(fig7, dpi=150)
        plt.close(fig7)

    return pdf_path, skew, auc


def main():
    ap = argparse.ArgumentParser(description="CRISPR library QC: KDE / Skew / AUC (per-sample)")
    ap.add_argument("-i", "--input", help="MAGeCK readcount 表 (tab 分隔)")
    ap.add_argument("-s", "--sample", help="样本列名 (单选); 与 --all 互斥")
    ap.add_argument("--all", action="store_true", help="遍历 readcount 表里所有样本列, 每样本一个 PDF")
    ap.add_argument("-o", "--outdir", default=".", help="输出目录")
    ap.add_argument("--png", action="store_true", help="同时保留各图 PNG")
    ap.add_argument("--demo", action="store_true", help="用模拟数据演示 (样本名=demo)")
    args = ap.parse_args()

    if args.demo:
        counts = demo_counts()
        pdf_path, skew, auc = make_pdf(counts, "demo", args.outdir, keep_png=args.png)
        print(f"demo  ->  {pdf_path}  Skew={skew:.2f}  AUC={auc:.4f}")
        return

    if not args.input:
        ap.error("需要 -i <readcount>, 或加 --demo")

    if args.all:
        samples = list_samples(args.input)
        print(f"发现 {len(samples)} 个样本, 开始批量生成 PDF -> {args.outdir}")
        for s in samples:
            counts = load_counts(args.input, s)
            pdf_path, skew, auc = make_pdf(counts, s, args.outdir, keep_png=args.png)
            flag_s = "PASS" if skew < 10 else "FAIL"
            flag_a = "PASS" if auc < 0.7 else "FAIL"
            print(f"  {s:18s} Skew={skew:7.2f} [{flag_s}]  AUC={auc:.4f} [{flag_a}]  -> {os.path.basename(pdf_path)}")
    else:
        if not args.sample:
            ap.error("需要 -s <sample>, 或加 --all 批量, 或加 --demo")
        counts = load_counts(args.input, args.sample)
        pdf_path, skew, auc = make_pdf(counts, args.sample, args.outdir, keep_png=args.png)
        flag_s = "PASS" if skew < 10 else "FAIL"
        flag_a = "PASS" if auc < 0.7 else "FAIL"
        print(f"sgRNAs (>0): {(counts > 0).sum()}")
        print(f"Skew ratio (p90/p10) = {skew:.2f}   {flag_s} (<10)")
        print(f"AUC (Lorenz)         = {auc:.4f}  {flag_a} (<0.7)")
        print(f"PDF written: {pdf_path}")


if __name__ == "__main__":
    main()
