#!/usr/bin/env python3
"""Generate status-labelled temperature and descriptor figures from v5 tables.

Use ``--cohort full16`` or ``--cohort matched3metal`` for a fixed cohort.
``--cohort available`` keeps every completed 5/100-bar pair and reports the
number of structures separately for each temperature.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parents[1]
OUT = HERE / "outputs"
TEMPERATURES = [77, 120, 160, 200, 233, 253, 273, 298]
COLORS = {"TiCo": "#C6B488", "TiMg": "#3D8064", "TiNi": "#C68E8E"}
DESCRIPTORS = [
    ("H2_POAV_cm3_g", "H$_2$ probe-occupiable pore volume (cm$^3$ g$^{-1}$)", "POAV"),
    ("H2_ASA_m2_g", "H$_2$-probe-accessible surface area (m$^2$ g$^{-1}$)", "ASA"),
    ("density_g_cm3", "Framework density (g cm$^{-3}$)", "density"),
    ("PLD_A", "Pore-limiting diameter (Å)", "PLD"),
]


def binned_medians(x: pd.Series, y: pd.Series) -> tuple[np.ndarray, np.ndarray]:
    frame = pd.DataFrame({"x": x, "y": y}).dropna().sort_values("x")
    if len(frame) < 20 or frame.x.nunique() < 4:
        return np.array([]), np.array([])
    bins = min(10, max(4, len(frame) // 25))
    frame["bin"] = pd.qcut(frame.x, q=bins, duplicates="drop")
    med = frame.groupby("bin", observed=True)[["x", "y"]].median()
    return med.x.to_numpy(), med.y.to_numpy()


def load(cohort: str) -> pd.DataFrame:
    name = {
        "available": "working_capacity_available_v5.csv",
        "full16": "working_capacity_full16_v5.csv",
        "matched3metal": "working_capacity_matched3metal_v5.csv",
    }[cohort]
    path = OUT / name
    if not path.is_file():
        raise FileNotFoundError(f"Run build_working_capacity_v5.py --cohort {cohort} first")
    return pd.read_csv(path)


def plot_temperature(frame: pd.DataFrame, cohort: str) -> None:
    target = OUT / f"figures_{cohort}"
    target.mkdir(parents=True, exist_ok=True)
    stats = []
    for t in TEMPERATURES:
        group = frame[frame.temperature_K.eq(t)]
        for quantity in ["loading_5bar_mg_g", "loading_100bar_mg_g", "working_capacity_mg_g"]:
            values = group[quantity].dropna()
            stats.append({
                "temperature_K": t, "quantity": quantity, "n": len(values),
                "q25": values.quantile(.25), "median": values.quantile(.5), "q75": values.quantile(.75),
            })
    summary = pd.DataFrame(stats)
    fig, axes = plt.subplots(1, 2, figsize=(12.2, 4.8), sharex=True)
    for quantity, color, label in [
        ("loading_100bar_mg_g", "#3D8064", "100 bar"),
        ("loading_5bar_mg_g", "#C6B488", "5 bar"),
    ]:
        q = summary[summary.quantity.eq(quantity)].set_index("temperature_K").reindex(TEMPERATURES)
        axes[0].fill_between(TEMPERATURES, q.q25, q.q75, color=color, alpha=.18)
        axes[0].plot(TEMPERATURES, q.median, marker="o", color=color, label=label)
    q = summary[summary.quantity.eq("working_capacity_mg_g")].set_index("temperature_K").reindex(TEMPERATURES)
    axes[1].fill_between(TEMPERATURES, q.q25, q.q75, color="#3D8064", alpha=.18)
    axes[1].plot(TEMPERATURES, q.median, marker="o", color="#3D8064", label="Median")
    axes[0].set_ylabel("Absolute uptake (mg H$_2$ g$^{-1}$ framework)")
    axes[1].set_ylabel("Working capacity (mg H$_2$ g$^{-1}$ framework)")
    for ax in axes:
        ax.set_xlabel("Temperature (K)")
        ax.spines[["top", "right"]].set_visible(False)
        ax.grid(axis="y", color="#E4E8E4")
        ax.legend(frameon=False)
    n_text = ", ".join(f"{t}: {int(n)}" for t, n in zip(TEMPERATURES, q.n))
    fig.suptitle(f"v5 H$_2$ adsorption; cohort={cohort}", y=.99)
    fig.text(.5, .01, f"n per temperature (working capacity): {n_text}; shading = 25th–75th percentile", ha="center", fontsize=8.5)
    fig.tight_layout(rect=(0, .06, 1, .94))
    fig.savefig(target / "temperature_response.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    summary.to_csv(target / "temperature_summary.csv", index=False)


def spearman(x: pd.Series, y: pd.Series) -> float:
    rx = x.rank(method="average").to_numpy(dtype=float)
    ry = y.rank(method="average").to_numpy(dtype=float)
    if len(rx) < 2 or np.std(rx) == 0 or np.std(ry) == 0:
        return float("nan")
    return float(np.corrcoef(rx, ry)[0, 1])


def plot_descriptors(frame: pd.DataFrame, cohort: str) -> None:
    target = OUT / f"figures_{cohort}"
    target.mkdir(parents=True, exist_ok=True)
    for column, label, stem in DESCRIPTORS:
        fig, axes = plt.subplots(2, 4, figsize=(14.2, 7.0), sharex=False, sharey=False)
        for ax, temperature in zip(axes.flat, TEMPERATURES):
            group = frame[frame.temperature_K.eq(temperature)]
            subset = group[[column, "working_capacity_mg_g", "metal_node"]].dropna()
            for metal, color in COLORS.items():
                m = subset[subset.metal_node.eq(metal)]
                ax.scatter(m[column], m.working_capacity_mg_g, s=16, alpha=.55, color=color, edgecolors="none")
            rho = spearman(subset[column], subset.working_capacity_mg_g)
            med_x, med_y = binned_medians(subset[column], subset.working_capacity_mg_g)
            if len(med_x):
                ax.plot(med_x, med_y, color="#333333", marker="o", markersize=2.8, linewidth=1.2)
            ax.text(.04, .94, f"rho = {rho:.3f}\nn = {len(subset)}", transform=ax.transAxes, va="top", fontsize=8.5,
                    bbox={"facecolor": "white", "alpha": .8, "edgecolor": "none"})
            ax.set_title(f"{temperature} K")
            ax.spines[["top", "right"]].set_visible(False)
            ax.grid(color="#E4E8E4", linewidth=.5)
        for ax in axes[-1, :]: ax.set_xlabel(label)
        for ax in axes[:, 0]: ax.set_ylabel("Working capacity (mg H$_2$ g$^{-1}$)")
        fig.suptitle(f"H$_2$ working capacity versus {label}; cohort={cohort}", y=.995)
        fig.tight_layout(rect=(0, 0, 1, .96))
        fig.savefig(target / f"{stem}_vs_working_capacity.png", dpi=300, bbox_inches="tight")
        plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cohort", choices=["available", "full16", "matched3metal"], default="full16")
    args = parser.parse_args()
    frame = load(args.cohort)
    plot_temperature(frame, args.cohort)
    plot_descriptors(frame, args.cohort)
    print(f"Wrote v5 figures for {args.cohort}: {frame.current_basename.nunique()} structures, {len(frame)} pairs")


if __name__ == "__main__":
    main()
