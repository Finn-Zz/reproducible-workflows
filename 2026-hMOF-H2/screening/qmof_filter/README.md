# QMOF local-geometry filter

This folder contains the common model and scoring code for the local-geometry filter described in Section 2.3/S5 of the manuscript. It is an upstream structural-plausibility analysis and is kept separate from the three metal-specific GCMC input bundles.

The filter compares bond lengths, bond angles, and organic dihedrals in a relaxed CIF with empirical QMOF reference distributions. The calibrated structure-level cutoff is **3.01256924351463** (99th percentile). A score is a geometric-rarity measure; it is not an energy, a thermodynamic-stability prediction, or a synthetic-feasibility prediction.

## Included files

- `model/qmof_geometry_v1.json.gz`: the exact calibrated model used for the reported scores.
- `model/qmof_geometry_v1.json.diagnostics.json`: model diagnostics.
- `src/mofrelax/` and `pyproject.toml`: the scoring implementation and its NumPy/ASE dependencies.
- `manifests/`: one 144-row candidate manifest and one 144-row score table for each of TiCo, TiMg, and TiNi. The three files together correspond to the 432 CIFs in `structures/`.
- `results/candidate_manifest_full_576.csv` and `results/qmof_scores_full_576.csv`: the complete 576-member filter metadata and score tables. The 144 Cr3F rows are retained for library-level provenance; their CIFs are outside this adsorption-focused package.
- `results/qmof_summary_by_connectivity_metal.csv`: the aggregate filter summary.
- `parameters.json`: the model and calibration settings used for interpretation.

The QMOF reference CIF archive is not redistributed here. The model already contains the fitted reference distributions; the original QMOF source/version and attribution are recorded in `parameters.json` and `CITATION.cff`.

The filter was used to characterize local-geometry plausibility across the generated library. Its outlier label should not be read as an automatic deletion rule for the adsorption dataset; the manuscript reports the structural and adsorption analyses together.

## Re-score the packaged structures

From the repository root, install the local filter package and run the wrapper:

```bash
python -m pip install -e screening/qmof_filter
python screening/qmof_filter/score_package.py \
  --metal all \
  --output /tmp/hmof_qmof_scores_recomputed.csv
```

The output contains both the normalized structure identifier and the RASPA framework basename, so it can be joined to the adsorption table without relying on local machine paths. The checked-in per-metal score tables are reference outputs for this calculation.

For portability, the model metadata stores the original reference-manifest location as the relative label `external/QMOF-v2021-12/manifest.csv`; the fitted reference distributions and calibration values are unchanged.
