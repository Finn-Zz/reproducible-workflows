# Conventional structural-plausibility checks

The three methods reported in the Supporting Information are independent validation outputs rather than GCMC input files:

- `mofclassifier_results.csv` records the learned crystal-likeness scores;
- `mofchecker_results.csv` records coordination, overlap, connectivity, and terminal-environment checks;
- `chen_manz_results.csv` records the Chen–Manz local-defect flags.

The files are split by TiCo, TiMg, and TiNi and retain only the 432 structures used in the adsorption package. Paths are relative to the repository and contain no machine-specific absolute paths. The QMOF local-geometry filter is documented separately in `screening/qmof_filter/` because it uses a continuous reference-distribution score rather than the binary validator criteria.
