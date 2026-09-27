#!/usr/bin/env python3
"""Build auditable H2 working-capacity tables from the v5 dataset.

The 432-framework design library is kept separate from the completed
adsorption observations.  By default this script writes all available
5/100-bar pairs and reports the number of structures at each temperature.
It also writes two strict cohorts: structures complete at all 16 conditions,
and the subset of complete topology-linker combinations represented by all
three heterometallic nodes.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

try:  # Optional: the rank coefficient itself has a NumPy fallback.
    from scipy.stats import spearmanr as _scipy_spearmanr
except ImportError:  # pragma: no cover - exercised only in minimal environments
    _scipy_spearmanr = None

HERE = Path(__file__).resolve().parents[1]
DATA = HERE / "data"
OUT = HERE / "outputs"
SUMMARY = DATA / "KISTI_adsorption_overall_summary_v5.csv"
ZEO = DATA / "zeopp_H2_complete_584_new_names.csv"
TEMPERATURES = [77, 120, 160, 200, 233, 253, 273, 298]
DESCRIPTORS = {
    "H2_POAV_cm3_g": "H2 probe-occupiable pore volume (cm3 g-1)",
    "H2_ASA_m2_g": "H2-probe-accessible surface area (m2 g-1)",
    "density_g_cm3": "Framework density (g cm-3)",
    "PLD_A": "Pore-limiting diameter (A)",
}


def spearman_values(x: pd.Series, y: pd.Series) -> tuple[float, float]:
    """Return Spearman rho and p when SciPy is available.

    The rank coefficient has a NumPy fallback so the core tables remain
    reproducible in a minimal environment. A missing SciPy installation only
    leaves ``p_value`` as NaN; it does not change ``spearman_rho``.
    """
    if _scipy_spearmanr is not None:
        result = _scipy_spearmanr(x, y, nan_policy="omit")
        return float(result.statistic), float(result.pvalue)
    rx = pd.Series(x).rank(method="average").to_numpy(dtype=float)
    ry = pd.Series(y).rank(method="average").to_numpy(dtype=float)
    if len(rx) < 2 or np.std(rx) == 0 or np.std(ry) == 0:
        return float("nan"), float("nan")
    return float(np.corrcoef(rx, ry)[0, 1]), float("nan")


def complete_mask(frame: pd.DataFrame) -> pd.Series:
    return (
        frame["completed"].astype(str).str.lower().eq("yes")
        & frame["result_parse_status"].eq("parsed")
        & frame["absolute_loading_mg_g_framework"].notna()
    )


def build_pairs(raw: pd.DataFrame, descriptors: pd.DataFrame) -> pd.DataFrame:
    valid = raw.loc[complete_mask(raw)].copy()
    pivot = valid.pivot_table(
        index=["current_basename", "temperature_K"],
        columns="pressure_bar",
        values="absolute_loading_mg_g_framework",
        aggfunc="first",
    ).rename(columns={5: "loading_5bar_mg_g", 100: "loading_100bar_mg_g"})
    pivot = pivot.dropna(subset=["loading_5bar_mg_g", "loading_100bar_mg_g"]).reset_index()
    pivot["working_capacity_mg_g"] = pivot["loading_100bar_mg_g"] - pivot["loading_5bar_mg_g"]
    pivot["working_capacity_wt_percent"] = pivot["working_capacity_mg_g"] / 10.0

    metadata = raw.drop_duplicates("current_basename")[[
        "current_basename", "connectivity", "topology", "metal_node", "organic_linker"
    ]]
    features = descriptors.rename(columns={"basename": "current_basename"})[[
        "current_basename", "density_g_cm3", "H2_ASA_m2_g", "H2_POAV_cm3_g",
        "H2_POAV_fraction", "H2_PSD_accessible_samples", "H2_PSD_accessible_fraction",
        "H2_PSD_mode_A", "H2_PSD_mean_A", "H2_PSD_median_A", "H2_PSD_p90_A",
        "LCD_A", "PLD_A", "LFPD_A", "H2_probe_radius_A",
    ]]
    merged = pivot.merge(metadata, on="current_basename", how="left", validate="many_to_one")
    merged = merged.merge(features, on="current_basename", how="left", validate="many_to_one")
    if merged["H2_POAV_cm3_g"].isna().any():
        missing = sorted(merged.loc[merged["H2_POAV_cm3_g"].isna(), "current_basename"].unique())
        raise ValueError(f"Missing Zeo++ descriptors for {len(missing)} working-capacity structures")
    if not np.isfinite(merged["working_capacity_mg_g"]).all():
        raise ValueError("Non-finite working capacities found")
    return merged.sort_values(["temperature_K", "current_basename"]).reset_index(drop=True)


def cohort_sets(raw: pd.DataFrame) -> dict[str, set[str]]:
    complete = complete_mask(raw)
    per_structure = raw.assign(_complete=complete).groupby("current_basename")["_complete"].all()
    full16 = set(per_structure[per_structure].index)
    metadata = raw.drop_duplicates("current_basename")
    full_meta = metadata[metadata.current_basename.isin(full16)]
    combo_counts = full_meta.groupby(["topology", "organic_linker", "connectivity"])["metal_node"].nunique()
    matched_combos = set(combo_counts[combo_counts.eq(3)].index)
    matched = set(
        full_meta.loc[
            full_meta.apply(
                lambda r: (r.topology, r.organic_linker, r.connectivity) in matched_combos,
                axis=1,
            ),
            "current_basename",
        ]
    )
    return {"design_432": set(metadata.current_basename), "full16": full16, "matched3metal": matched}


def write_statistics(name: str, frame: pd.DataFrame) -> None:
    summary_rows = []
    corr_rows = []
    for temperature in TEMPERATURES:
        group = frame[frame.temperature_K.eq(temperature)]
        for quantity in ["loading_5bar_mg_g", "loading_100bar_mg_g", "working_capacity_mg_g"]:
            values = group[quantity].dropna()
            summary_rows.append({
                "temperature_K": temperature,
                "quantity": quantity,
                "n": int(values.size),
                "q25": float(values.quantile(.25)) if len(values) else np.nan,
                "median": float(values.quantile(.50)) if len(values) else np.nan,
                "q75": float(values.quantile(.75)) if len(values) else np.nan,
            })
        for descriptor in DESCRIPTORS:
            subset = group[[descriptor, "working_capacity_mg_g"]].dropna()
            if len(subset) < 3 or subset[descriptor].nunique() < 2:
                rho, p_value = np.nan, np.nan
            else:
                rho, p_value = spearman_values(subset[descriptor], subset.working_capacity_mg_g)
            corr_rows.append({
                "temperature_K": temperature,
                "descriptor": descriptor,
                "spearman_rho": rho,
                "p_value": p_value,
                "n": int(len(subset)),
            })
    pd.DataFrame(summary_rows).to_csv(OUT / f"summary_{name}.csv", index=False)
    pd.DataFrame(corr_rows).to_csv(OUT / f"spearman_{name}.csv", index=False)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--cohort", choices=["available", "full16", "matched3metal", "all"], default="available",
        help="available=all completed pairs; full16=395 structures complete at all conditions; "
             "matched3metal=complete three-metal combinations; all=432 design members (fails if incomplete).",
    )
    args = parser.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    raw = pd.read_csv(SUMMARY)
    descriptors = pd.read_csv(ZEO)
    if descriptors.basename.nunique() < 432:
        raise ValueError("The descriptor table does not cover the 432-structure design library")
    pairs = build_pairs(raw, descriptors)
    sets = cohort_sets(raw)
    if args.cohort == "available":
        selected = pairs
    else:
        names = sets["design_432" if args.cohort == "all" else args.cohort]
        selected = pairs[pairs.current_basename.isin(names)].copy()
        if args.cohort == "all" and len(selected) != 432 * len(TEMPERATURES):
            raise RuntimeError(
                "The 432 design library is not complete at every temperature and pressure; "
                "use --cohort available/full16/matched3metal or finish the missing simulations."
            )
    output_name = {
        "available": "working_capacity_available_v5.csv",
        "full16": "working_capacity_full16_v5.csv",
        "matched3metal": "working_capacity_matched3metal_v5.csv",
        "all": "working_capacity_design432_v5.csv",
    }[args.cohort]
    selected.to_csv(OUT / output_name, index=False)
    write_statistics(args.cohort, selected)
    status = {
        "cohort": args.cohort,
        "structures_in_design_library": len(sets["design_432"]),
        "structures_selected": int(selected.current_basename.nunique()),
        "working_capacity_rows": int(len(selected)),
        "rows_by_temperature": {
            str(t): int(selected[selected.temperature_K.eq(t)].current_basename.nunique())
            for t in TEMPERATURES
        },
        "source": SUMMARY.name,
        "descriptor_source": ZEO.name,
        "working_capacity_definition": "absolute loading at 100 bar minus absolute loading at 5 bar",
    }
    (OUT / f"analysis_status_{args.cohort}.json").write_text(
        json.dumps(status, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(status, indent=2))


if __name__ == "__main__":
    main()
