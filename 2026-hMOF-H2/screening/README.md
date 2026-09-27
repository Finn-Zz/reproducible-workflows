# Screening and validation provenance

The screening files are kept separate from the RASPA3 input bundles because they answer a different reproducibility question:

- `topology/` contains the coordination-figure and topology/building-block compatibility records used before structure generation.
- `qmof_filter/` contains the calibrated local-geometry model, portable scoring code, parameters, and score/manifest tables for the generated library.
- `../validation/` contains the independent MOFClassifier, MOFChecker, and Chen–Manz result tables.

The `gcmc/TiCo`, `gcmc/TiMg`, and `gcmc/TiNi` folders are the files needed to run the adsorption calculation. The screening folders document how the structures were assessed and should not be interpreted as alternate force-field or simulation inputs.
