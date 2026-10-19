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
- [`controls/coordination_figure/`](controls/coordination_figure/README.md): 448 completed UFF4MOF-optimized CIFs on coordination-figure control nets, the complete 624-design structure-check CSV, and an explanation of why the control set was generated and how it differs from the preliminary L2 comparison.
- The 584-row Zeo++ master table, including all 432 current v5 structure IDs.
- Validation, working-capacity, correlation, and plotting scripts, together with the generated analysis tables.
- [`gcmc/input/`](gcmc/input/) and [`gcmc/output/`](gcmc/output/): paired RASPA3 input and result folders for **all 432 frameworks at 298 K and 5 bar**, organized into 144 structures each for TiCo, TiMg, and TiNi. The [GCMC README](gcmc/README.md) explains the files, conditions, and source manifests. Three separate 77 K, 5 bar examples remain in `gcmc/TiCo/`, `gcmc/TiMg/`, and `gcmc/TiNi/`.
- `data/KISTI_adsorption_overall_summary_v5.csv`: the full 432-design × 8-temperature × 2-pressure results/status table (6,912 rows at this snapshot: 6,822 completed and 90 pending). The available-condition working-capacity and correlation tables are in `outputs/`; pending conditions are not treated as completed calculations.
- `screening/qmof_filter/`: the calibrated QMOF local-geometry model, scoring code, parameters, and per-metal score/manifest tables.
- `screening/topology/`: the coordination-figure, building-block/topology compatibility tables, and 14 compact net definitions used upstream.
- `validation/`: the three conventional structural-plausibility result sets for the 432 adsorption structures.
- [`relaxation/`](relaxation/README.md): the Materials Studio Forcite/UFF macro and the LAMMPS/UFF4MOF preparation, relayed minimization, and bonded-topology restoration scripts used before GCMC.
- The repository includes the RASPA3 result JSON/text for all 432 frameworks at 298 K and 5 bar, plus three earlier 77 K, 5 bar examples. Binary restart states, cluster scheduler files, and the raw archives for the other conditions remain outside this repository.
