# Geometry-Aware Construction of Hypothetical Metal–Organic Frameworks through Coordination-Figure Matching: Hydrogen Adsorption as an Application

**Authors:** Pengyu Zhao, Taekgi Lee, Myoung Soo Lah, and Yongchul G. Chung  
**Repository scope:** computational data and analysis accompanying the hMOF–H₂ manuscript

This repository contains the current, professor-facing hMOF–H₂ adsorption analysis. It was exported from a larger local working archive, so old scripts, temporary files, and cluster-specific paths are intentionally excluded.

## What is included

- **432 designed heterometallic frameworks**: 144 topology–linker combinations × TiCo, TiMg, and TiNi.
- **6,912 adsorption states**: 8 temperatures × 2 pressures.
- **432 normalized CIF inputs**, with file hashes in `structures/structure_manifest_432.csv`.
- [`structures/README.md`](structures/README.md): the normalized CIF inputs and the `structures/BBs/` connection-oriented building-block/node coordinate inputs used for structure generation.
- The 584-row Zeo++ master table, including all 432 current v5 structure IDs.
- Validation, working-capacity, correlation, and plotting scripts, together with the generated analysis tables.
- `gcmc/TiCo`, `gcmc/TiMg`, and `gcmc/TiNi`: the metal-specific RASPA3 input bundles (`h2.json`, `force_field.json`, `simulation.json`, and a matching example CIF).
- `screening/qmof_filter/`: the calibrated QMOF local-geometry model, scoring code, parameters, and per-metal score/manifest tables.
- `screening/topology/`: the coordination-figure, building-block/topology compatibility tables, and 14 compact net definitions used upstream.
- `validation/`: the three conventional structural-plausibility result sets for the 432 adsorption structures.
- [`relaxation/`](relaxation/README.md): the Materials Studio Forcite/UFF macro and the LAMMPS/UFF4MOF preparation, relayed minimization, and bonded-topology restoration scripts used before GCMC.
