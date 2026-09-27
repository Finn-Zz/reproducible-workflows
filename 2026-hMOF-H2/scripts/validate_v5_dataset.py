#!/usr/bin/env python3
"""Validate the v5 adsorption table and write an auditable cohort manifest.

The v5 table describes the generated 432-framework design library.  A row is
marked ``complete`` only when the source parser found a finished RASPA output;
incomplete rows are retained in the manifest and are never imputed.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd


HERE = Path(__file__).resolve().parents[1]
DATA = HERE / "data"
OUT = HERE / "outputs"
SUMMARY = DATA / "KISTI_adsorption_overall_summary_v5.csv"
ZEO = DATA / "zeopp_H2_complete_584_new_names.csv"
STRUCTURES = HERE / "structures"
STRUCTURE_MANIFEST = STRUCTURES / "structure_manifest_432.csv"

TEMPERATURES = [77, 120, 160, 200, 233, 253, 273, 298]
PRESSURES = [5, 100]
METALS = {"TiCo", "TiMg", "TiNi"}


def _as_bool(series: pd.Series) -> pd.Series:
    return series.astype(str).str.lower().eq("yes")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--verify-structure-hashes", action="store_true",
        help="recompute the SHA-256 hash of each bundled CIF",
    )
    args = parser.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(SUMMARY)
    descriptors = pd.read_csv(ZEO, usecols=[
        "basename", "relative_path", "result_source", "source_basename",
        "density_g_cm3", "H2_ASA_m2_g", "H2_POAV_cm3_g",
        "PLD_A", "H2_probe_radius_A", "all_zeopp_outputs_complete",
    ])

    required = {
        "current_basename", "basename", "topology", "metal_node", "organic_linker",
        "temperature_K", "pressure_bar", "completed", "result_parse_status",
        "absolute_loading_mg_g_framework",
    }
    missing = required.difference(df.columns)
    if missing:
        raise ValueError(f"Missing v5 columns: {sorted(missing)}")

    expected_rows = 432 * len(TEMPERATURES) * len(PRESSURES)
    if len(df) != expected_rows:
        raise ValueError(f"Expected {expected_rows} v5 rows, found {len(df)}")
    if df["current_basename"].nunique() != 432:
        raise ValueError("The v5 table does not contain exactly 432 structures")
    if df[["current_basename", "temperature_K", "pressure_bar"]].duplicated().any():
        raise ValueError("Duplicate structure/temperature/pressure rows found")
    if set(df["temperature_K"]) != set(TEMPERATURES) or set(df["pressure_bar"]) != set(PRESSURES):
        raise ValueError("Unexpected temperature or pressure values")
    if set(df["metal_node"]) != METALS:
        raise ValueError(f"Unexpected metal nodes: {sorted(set(df['metal_node']))}")

    # The package includes a normalized CIF for every current v5 identifier.
    # This check catches accidental omission or renaming of a rerun input.
    cif_paths = sorted(STRUCTURES.rglob("*.cif"))
    cif_ids = {path.stem for path in cif_paths}
    design_ids = set(df["current_basename"])
    if len(cif_paths) != 432 or cif_ids != design_ids:
        raise ValueError(
            f"Expected one normalized CIF for each of 432 v5 structures; "
            f"found {len(cif_paths)} files, {len(cif_ids & design_ids)} matching IDs"
        )
    if not STRUCTURE_MANIFEST.is_file():
        raise ValueError(f"Missing structure manifest: {STRUCTURE_MANIFEST}")
    structure_manifest = pd.read_csv(STRUCTURE_MANIFEST)
    if set(structure_manifest["current_basename"]) != design_ids or len(structure_manifest) != 432:
        raise ValueError("The 432-CIF structure manifest does not match the v5 design IDs")
    if args.verify_structure_hashes:
        for row in structure_manifest.itertuples(index=False):
            path = STRUCTURES / row.relative_path
            if not path.is_file():
                raise ValueError(f"Missing CIF listed in structure manifest: {path}")
            if path.stat().st_size != int(row.bytes):
                raise ValueError(f"CIF size changed: {path}")
            digest = hashlib.sha256(path.read_bytes()).hexdigest()
            if digest != row.sha256:
                raise ValueError(f"CIF hash changed: {path}")

    structure_meta = df.drop_duplicates("current_basename").copy()
    combo_counts = structure_meta.groupby(
        ["topology", "organic_linker", "connectivity"], dropna=False
    )["metal_node"].nunique()
    if len(combo_counts) != 144 or not (combo_counts == 3).all():
        raise ValueError(
            f"Expected 144 topology-linker combinations with three metals; "
            f"found {len(combo_counts)} combinations and counts {combo_counts.value_counts().to_dict()}"
        )

    complete = _as_bool(df["completed"]) & df["result_parse_status"].eq("parsed")
    df["_complete"] = complete
    by_structure = df.groupby("current_basename", sort=True)
    completed_counts = by_structure["_complete"].sum().astype(int)
    complete_temperatures = (
        df.loc[complete]
        .groupby("current_basename")["temperature_K"]
        .agg(lambda values: ",".join(str(int(x)) for x in sorted(set(values))))
    )
    missing_conditions = (
        df.loc[~complete]
        .assign(condition=lambda x: x.temperature_K.astype(str) + "K/" + x.pressure_bar.astype(str) + "bar")
        .groupby("current_basename")["condition"]
        .agg(lambda values: ";".join(values))
    )

    manifest = structure_meta[[
        "current_basename", "basename", "connectivity", "topology",
        "organic_linker", "metal_node",
    ]].copy()
    manifest["design_library_member"] = True
    manifest["completed_conditions"] = manifest["current_basename"].map(completed_counts).fillna(0).astype(int)
    manifest["complete_all_16_conditions"] = manifest["completed_conditions"].eq(16)
    manifest["complete_temperatures"] = manifest["current_basename"].map(complete_temperatures).fillna("")
    manifest["missing_conditions"] = manifest["current_basename"].map(missing_conditions).fillna("")
    manifest["all_zeopp_outputs_complete"] = manifest["current_basename"].isin(set(descriptors["basename"]))

    # A matched cohort is defined separately from the design library: every
    # member of the topology-linker combination must have all three metals and
    # all 16 adsorption conditions complete.
    full_names = set(manifest.loc[manifest.complete_all_16_conditions, "current_basename"])
    full_meta = manifest[manifest.current_basename.isin(full_names)]
    full_combo_counts = full_meta.groupby(
        ["topology", "organic_linker", "connectivity"], dropna=False
    )["metal_node"].nunique()
    complete_matched_combos = set(full_combo_counts.index[full_combo_counts.eq(3)])
    manifest["complete_three_metal_combo"] = manifest.apply(
        lambda row: row.complete_all_16_conditions
        and (row.topology, row.organic_linker, row.connectivity) in complete_matched_combos,
        axis=1,
    )
    manifest.to_csv(OUT / "cohort_manifest_v5.csv", index=False)

    # Keep the two naming layers explicit: current_basename is the v5/RASPA
    # identifier, while basename/source_basename are the descriptor-table
    # identifiers. This mapping prevents accidental joins to the legacy name.
    descriptor_map = descriptors[["basename", "relative_path", "result_source", "source_basename"]].copy()
    descriptor_map = descriptor_map.rename(columns={"basename": "descriptor_basename", "relative_path": "descriptor_relative_path"})
    cif_map = structure_manifest[["current_basename", "relative_path", "bytes", "sha256"]].rename(
        columns={"relative_path": "cif_relative_path", "bytes": "cif_bytes", "sha256": "cif_sha256"}
    )
    structure_descriptor_manifest = manifest.merge(
        descriptor_map, left_on="current_basename", right_on="descriptor_basename",
        how="left", validate="one_to_one"
    ).merge(cif_map, on="current_basename", how="left", validate="one_to_one")
    if structure_descriptor_manifest["descriptor_relative_path"].isna().any() or structure_descriptor_manifest["cif_relative_path"].isna().any():
        raise ValueError("Structure/descriptor mapping is incomplete")
    structure_descriptor_manifest.to_csv(OUT / "structure_descriptor_manifest_v5.csv", index=False)

    paired_complete_by_temperature = {}
    for temperature in TEMPERATURES:
        at_temperature = df.loc[complete & df.temperature_K.eq(temperature)]
        pressure_counts = at_temperature.groupby("current_basename")["pressure_bar"].nunique()
        paired_complete_by_temperature[str(temperature)] = int((pressure_counts == len(PRESSURES)).sum())

    status = {
        "source": SUMMARY.name,
        "design_structures": int(df["current_basename"].nunique()),
        "design_topology_linker_combinations": int(len(combo_counts)),
        "design_rows": int(len(df)),
        "temperatures_K": TEMPERATURES,
        "pressures_bar": PRESSURES,
        "complete_rows": int(complete.sum()),
        "incomplete_rows": int((~complete).sum()),
        "complete_structures_all_16_conditions": int(manifest.complete_all_16_conditions.sum()),
        "complete_three_metal_combinations": int(len(complete_matched_combos)),
        "complete_three_metal_structures": int(manifest.complete_three_metal_combo.sum()),
        "paired_complete_structures_by_temperature": paired_complete_by_temperature,
        "incomplete_by_status": df.loc[~complete, "result_parse_status"].value_counts().to_dict(),
        "descriptor_rows_available": int(len(descriptors)),
        "descriptor_matches_to_design_library": int(
            set(df.current_basename).intersection(descriptors.basename).__len__()
        ),
        "normalized_cif_files": int(len(cif_paths)),
        "normalized_cif_matches_to_design_library": int(len(cif_ids & design_ids)),
        "structure_descriptor_manifest_rows": int(len(structure_descriptor_manifest)),
        "note": (
            "The 432 count is the generated design-library size. Complete adsorption "
            "counts are reported separately; missing states are not imputed."
        ),
    }
    (OUT / "v5_dataset_status.json").write_text(
        json.dumps(status, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    print(json.dumps(status, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
