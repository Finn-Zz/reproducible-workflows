# RASPA3 H₂ input and output files

The [`input/`](input/) and [`output/`](output/) directories contain the **298 K, 5 bar** run for every one of the 432 designed TiCo, TiMg, and TiNi frameworks. Each directory is divided into `TiCo/`, `TiMg/`, and `TiNi/`, with **144 named framework folders per metal**. The matching folders in `input/` and `output/` refer to the same run and retain its original framework basename.

```text
gcmc/
├── input/TiCo/<framework>/   (and TiMg, TiNi)
│   ├── <framework>.cif
│   ├── h2.json
│   ├── force_field.json
│   └── simulation.json
└── output/TiCo/<framework>/  (and TiMg, TiNi)
    ├── the same four input files
    ├── output_298/            (RASPA3 result JSON and text, plus restart JSON)
    ├── bias_factors_298/
    └── .raspa_success_298
```

[`cohort_298K_5bar_manifest.csv`](cohort_298K_5bar_manifest.csv) maps each of the 432 framework folders to its original run in one of four source archives: `results/298.zip`, `SI.zip`, `stp.zip`, or `need-rerun.zip`. [`cohort_298K_5bar_file_manifest.csv`](cohort_298K_5bar_file_manifest.csv) provides a SHA-256 checksum for every copied file. All 432 copied runs have `ExternalTemperature = 298`, `ExternalPressure = 500000` Pa, a framework `Name` matching the CIF filename, a RASPA result JSON and text file, and an `exit_code=0` success marker.

The output folders contain the input files as well so that each result can be inspected with the exact settings and force field used to produce it. The large binary `restart_data_298.bin` files, PBS scheduler scripts, and scheduler logs were not copied into this GitHub package. The RASPA result JSON/text, restart JSON, bias factors, and success markers **are** included. This package covers one state per framework; it is not the raw output archive for all temperatures and pressures.

## Using another temperature or pressure

Copy the **input** folder for the appropriate metal and framework, then edit `Systems[0].ExternalTemperature` (K) and `Systems[0].ExternalPressure` (Pa) in `simulation.json`. For example, 5 bar is `500000` Pa and 100 bar is `10000000` Pa. Keep the matching CIF, `h2.json`, and metal-specific `force_field.json` with it. If you replace the framework CIF, also change `Systems[0].Name` to its filename stem. The [`generate_simulation_json.py`](generate_simulation_json.py) helper can make this edit from an existing template. **Changing an input does not change an existing result:** RASPA3 must be run again for the new state.

The original calculation grid was 77, 120, 160, 200, 233, 253, 273, and 298 K at 5 and 100 bar. The full state-by-state [`KISTI_adsorption_overall_summary_v5.csv`](../data/KISTI_adsorption_overall_summary_v5.csv) contains 6,912 planned rows (6,822 completed, 90 pending in this snapshot). Its `completed` field identifies rows with available results.

## Earlier example bundles and cell size

The separate top-level [`TiCo/`](TiCo/), [`TiMg/`](TiMg/), and [`TiNi/`](TiNi/) folders are the previously supplied **77 K, 5 bar** examples. Their source files and checksums are listed in [`input_manifest.csv`](input_manifest.csv) and [`result_manifest.csv`](result_manifest.csv); they are distinct from the 432-framework set above.

`NumberOfUnitCells = [1, 1, 1]` in these RASPA3 inputs means one copy of an **already expanded and relaxed CIF**, rather than one copy of the original generated cell. For the three older examples, [`cell_expansion_evidence.csv`](cell_expansion_evidence.csv) and the `cell_reference/` CIFs show a 2 × 2 × 2 expansion of the generated cells before the final GCMC structures were produced. Their final cell lengths are about 31–34 Å, and the reported RASPA3 framework–H₂ van der Waals cutoff is 12 Å. These particular examples do not support a 3 × 3 × 3 expansion claim.
