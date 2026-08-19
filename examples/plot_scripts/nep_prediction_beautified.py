#!/usr/bin/env python3
"""HQmdstemkit beautified 2x2 NEP parity plot (energy/force/virial/stress).
Reads *_train.out and optional *_test.out in the current directory.
Output: nep_prediction_2x2.png
"""
import os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial", "DejaVu Sans", "Liberation Sans"]
plt.rcParams["axes.unicode_minus"] = False
plt.rcParams["axes.labelsize"] = 13
plt.rcParams["xtick.labelsize"] = 12
plt.rcParams["ytick.labelsize"] = 12
plt.rcParams["axes.linewidth"] = 1.5
plt.rcParams["xtick.major.width"] = 1.2
plt.rcParams["ytick.major.width"] = 1.2


def rmse(a, b):
    return float(np.sqrt(np.mean((a - b) ** 2)))


def mae(a, b):
    return float(np.mean(np.abs(a - b)))


def r2(true, pred):
    true = np.asarray(true).reshape(-1)
    pred = np.asarray(pred).reshape(-1)
    ss_res = np.sum((true - pred) ** 2)
    ss_tot = np.sum((true - np.mean(true)) ** 2)
    return float(1.0 - ss_res / ss_tot) if ss_tot != 0 else float("nan")


def limits(true, pred, pad=0.05):
    lo = min(np.min(true), np.min(pred))
    hi = max(np.max(true), np.max(pred))
    rng = (hi - lo) or 1.0
    return lo - pad * rng, hi + pad * rng


def beauty(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(axis="both", colors="black")


def marginals(ax, true, pred, color, bins=35):
    d = make_axes_locatable(ax)
    top = d.append_axes("top", size="16%", pad=0.05, sharex=ax)
    right = d.append_axes("right", size="16%", pad=0.05, sharey=ax)
    top.hist(true, bins=bins, color=color, alpha=0.55, edgecolor="gray", linewidth=0.6)
    right.hist(pred, bins=bins, orientation="horizontal", color=color, alpha=0.55,
               edgecolor="gray", linewidth=0.6)
    for a in (top, right):
        a.grid(False)
        a.spines["top"].set_visible(False)
        a.spines["right"].set_visible(False)
    top.tick_params(axis="x", labelbottom=False, bottom=False)
    top.tick_params(axis="y", left=False, labelleft=False)
    top.spines["left"].set_visible(False)
    right.tick_params(axis="y", labelleft=False, left=False)
    right.tick_params(axis="x", bottom=False, labelbottom=False)
    right.spines["bottom"].set_visible(False)


def residual_inset(ax, res, color):
    ins = ax.inset_axes([0.58, 0.18, 0.30, 0.18])
    ins.hist(res, bins=28, color=color, alpha=0.70, edgecolor="gray", linewidth=0.6)
    mean = np.mean(res)
    ins.axvline(mean, color="k", ls="--", lw=1.0)
    lo, hi = ins.get_xlim()
    ins.text(mean + 0.04 * (hi - lo), ins.get_ylim()[1] * 0.82, f"{mean:.2f}",
             color=color, fontsize=9, fontweight="bold")
    ins.set_xlabel("Residual", fontsize=9, labelpad=1)
    ins.set_yticks([])
    ins.tick_params(axis="x", labelsize=8, pad=1)
    for s in ("top", "right", "left"):
        ins.spines[s].set_visible(False)
    ins.patch.set_alpha(0.0)


def panel(ax, true, pred, xlab, ylab, color, unit, scale_meV=False):
    true = np.asarray(true).reshape(-1)
    pred = np.asarray(pred).reshape(-1)
    x0, x1 = limits(true, pred)
    ax.scatter(true, pred, s=28, c=color, alpha=0.35, edgecolors="none", rasterized=True)
    ax.plot([x0, x1], [x0, x1], color="grey", ls="--", lw=2.0)
    ax.set_xlim(x0, x1); ax.set_ylim(x0, x1)
    ax.set_xlabel(xlab); ax.set_ylabel(ylab)
    beauty(ax)
    r = r2(true, pred)
    rm = rmse(pred, true); ma = mae(pred, true)
    if scale_meV:
        rm, ma = rm * 1000.0, ma * 1000.0
        unit_mev = unit.replace("eV", "meV") if "eV" in unit else "meV"
        rm_txt, ma_txt = f"RMSE = {rm:.2f} {unit_mev}", f"MAE = {ma:.2f} {unit_mev}"
    else:
        rm_txt, ma_txt = f"RMSE = {rm:.4f} {unit}", f"MAE = {ma:.4f} {unit}"
    ax.text(0.05, 0.96, f"$R^2 = {r:.4f}$", transform=ax.transAxes, fontsize=11, va="top")
    ax.text(0.05, 0.86, ma_txt, transform=ax.transAxes, fontsize=11, va="top")
    ax.text(0.05, 0.76, rm_txt, transform=ax.transAxes, fontsize=11, va="top")
    marginals(ax, true, pred, color)
    residual_inset(ax, pred - true, color)


def load(dtype, tag):
    fn = f"{dtype}_{tag}.out"
    if not os.path.isfile(fn):
        return None
    d = np.loadtxt(fn)
    if d.ndim == 1:
        d = d.reshape(1, -1)
    half = d.shape[1] // 2
    return d[:, :half], d[:, half:]


def main():
    items = [
        ("energy", "#1f77b4", "eV/atom", True),
        ("force", "#2ca02c", r"eV/$\mathrm{\AA}$", True),
        ("virial", "#ff7f0e", "eV/atom", True),
        ("stress", "#d62728", "GPa", False),
    ]
    fig, axes = plt.subplots(2, 2, figsize=(12, 10), dpi=150)
    for ax, (dtype, color, unit, mev) in zip(axes.flat, items):
        train = load(dtype, "train")
        if train is None:
            ax.text(0.5, 0.5, f"no {dtype}_train.out", ha="center", va="center",
                    transform=ax.transAxes)
            continue
        nep, dft = train
        panel(ax, dft, nep, f"DFT {dtype} ({unit})", f"NEP {dtype} ({unit})",
              color, unit, mev)
        test = load(dtype, "test")
        if test is not None:
            ax.scatter(test[1].reshape(-1), test[0].reshape(-1), s=18,
                       facecolors="none", edgecolors="black", linewidths=0.8, alpha=0.8)
    for i, ax in enumerate(axes.flat):
        ax.text(-0.12, 1.05, f"({chr(97 + i)})", transform=ax.transAxes,
                fontsize=16, ha="left", va="bottom")
    plt.subplots_adjust(left=0.09, right=0.98, bottom=0.07, top=0.97,
                        wspace=0.28, hspace=0.30)
    out = "nep_prediction_2x2.png"
    if len(sys.argv) > 1 and sys.argv[1] == "--out":
        out = sys.argv[2]
    fig.savefig(out, dpi=300, bbox_inches="tight", facecolor="white")
    print(f"saved -> {out}")


if __name__ == "__main__":
    main()