# Geometry-Aware Construction of Hypothetical Metal–Organic Frameworks through Coordination-Figure Matching: Hydrogen Adsorption as an Application

**Authors:** Pengyu Zhao, Taekgi Lee, Myoung Soo Lah, and Yongchul G. Chung  
**Repository scope:** computational data and analysis accompanying the hMOF–H₂ manuscript

This repository contains the current, professor-facing hMOF–H₂ adsorption analysis. It was exported from a larger local working archive, so old scripts, temporary files, and cluster-specific paths are intentionally excluded.

## What is included

- **432 designed heterometallic frameworks**: 144 topology–linker combinations × TiCo, TiMg, and TiNi.
- **6,912 adsorption states**: 8 temperatures × 2 pressures.
- v5 status at this snapshot: **6,822 complete rows and 90 pending rows**.
- **432 normalized CIF inputs**, with file hashes in `structures/structure_manifest_432.csv`.
- [`structures/README.md`](structures/README.md): the normalized CIF inputs and the `structures/BBs/` connection-oriented building-block/node coordinate inputs used for structure generation.
- The 584-row Zeo++ master table, including all 432 current v5 structure IDs.
- Validation, working-capacity, correlation, and plotting scripts, together with the generated analysis tables.
- `gcmc/TiCo`, `gcmc/TiMg`, and `gcmc/TiNi`: the metal-specific RASPA3 input bundles (`h2.json`, `force_field.json`, `simulation.json`, and a matching example CIF).
- `screening/qmof_filter/`: the calibrated QMOF local-geometry model, scoring code, parameters, and per-metal score/manifest tables.
- `screening/topology/`: the coordination-figure, building-block/topology compatibility tables, and 14 compact net definitions used upstream.
- `validation/`: the three conventional structural-plausibility result sets for the 432 adsorption structures.
- [`relaxation/`](relaxation/README.md): the Materials Studio Forcite/UFF macro and the LAMMPS/UFF4MOF preparation, relayed minimization, and bonded-topology restoration scripts used before GCMC.

The manuscript-level generated library contains 576 hMOFs; this package focuses on the 432 heterometallic frameworks used for the adsorption analysis. The 432 count is the design-library size for that subset, while the pending adsorption states are identified explicitly and are never imputed.

The checked-in package is portable for the supplied tables, normalized CIFs, screening provenance, and representative GCMC inputs when placed as a project folder inside the group's `reproducible-workflows` repository. The multi-gigabyte raw RASPA archives and the external QMOF reference-CIF archive remain separately hosted. The commands below should be run from this directory.

## Run the current analysis

Install Python 3.10+ and the listed dependencies:

```bash
python -m pip install -r requirements.txt
```

Run from this repository root:

```bash
python scripts/validate_v5_dataset.py --verify-structure-hashes
python scripts/build_working_capacity_v5.py --cohort available
python scripts/build_working_capacity_v5.py --cohort full16
python scripts/build_working_capacity_v5.py --cohort matched3metal
python scripts/plot_v5_figures.py --cohort full16
```

The three cohort choices are deliberate:

- `available`: every completed 5/100-bar pair, with the sample count reported separately at each temperature;
- `full16`: structures complete at all 16 temperature/pressure conditions (395 structures in this snapshot);
- `matched3metal`: complete three-metal topology–linker combinations (384 structures in this snapshot).

The `all` option stops until all 432 structures have complete adsorption data. No missing state is imputed.

## Raw RASPA archives

The multi-gigabyte RASPA archives are not committed to this Git repository. Their download location should be added to [`raw_manifest/README.md`](raw_manifest/README.md) before sending the repository link. The expected archive sizes and SHA-256 values are in [`raw_manifest/raw_archives_manifest.tsv`](raw_manifest/raw_archives_manifest.tsv).

The optional `scripts/update_summary_v5_from_2round_1.py` utility can rebuild v5 from the raw second-round ZIP and the v4 table. The v4 table is included in `data/`; the raw ZIP must be downloaded separately and passed with `--archive`:

```bash
python scripts/update_summary_v5_from_2round_1.py \
  --archive raw/"2round (1).zip" \
  --source data/KISTI_adsorption_overall_summary_v4.csv \
  --target data/KISTI_adsorption_overall_summary_v5_rebuilt.csv \
  --report data/2round_1_update_v5_rebuilt_report.md
```

The bundled `build_kisti_overall_summary.py` file is a parser dependency for that optional rebuild utility; it is not the normal analysis entry point.

## Checksums and scope

From the repository root, verify the package files with:

```bash
sha256sum -c SHA256SUMS.txt
```

The package covers the v5 adsorption summary, normalized CIF inputs, Zeo++ descriptor linkage, the checked-in RASPA3 force-field/simulation inputs, topology and local-geometry screening provenance, working-capacity analysis, validation, and the corresponding temperature/descriptor plots. The three GCMC directories contain representative production input bundles; the complete multi-gigabyte stdout/restart archive remains external and is listed in `raw_manifest/`.

The compact building-block and net files are included for provenance and inspection; the full PORMAKE environment, optimizer logs, and cluster-specific structure-generation workflow remain outside this curated package.

Before publication, choose one named cohort, regenerate the relevant figures, and record the cohort in the manuscript caption or figure metadata.
