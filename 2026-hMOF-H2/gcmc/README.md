# RASPA3 H₂ input bundles

This directory contains the input files used by the RASPA3 hydrogen grand-canonical Monte Carlo calculations. The three metal-node families are kept in separate directories because their framework force fields are different:

```text
gcmc/
├── TiCo/
├── TiMg/
└── TiNi/
```

Each directory contains one complete, directly inspectable example bundle from the production workflow:

- the matching framework CIF;
- `h2.json`, the rigid TraPPE H₂ definition;
- `force_field.json`, including the UFF/DREIDING framework terms, TraPPE H₂ terms, Lorentz–Berthelot mixing, and the fitted Morse cross interactions for that metal pair;
- `simulation.json`, containing the Monte Carlo cycle counts, temperature, pressure, charge treatment, and move probabilities.
- `results/77K_5bar/`, the matching RASPA3 run output: human-readable and JSON output, a JSON restart record, the binary restart state, bias factors, and the run-success marker. `result_manifest.csv` records the source run path and SHA-256 checksum of each result file.

The example is the `ceq–N2` framework at 77 K and 5 bar. The CIF and the three JSON files are copied from the same production run directory, so the `Name`, framework filename, and local `ForceField`/`MoleculeDefinition` references agree.

The `NumberOfUnitCells = [1, 1, 1]` entry in these RASPA3 inputs means one copy of the **already expanded and relaxed CIF**, not one copy of the original generated cell. For these three examples, the original generated cell lengths are about 15–16 Å, the selected expanded cell lengths are approximately twice those values, and the GCMC CIF cell lengths are about 31–34 Å after relaxation. Thus a further replication in RASPA3 was not used. The source files for these examples support a **2 × 2 × 2** expansion; they do not support describing these particular runs as 3 × 3 × 3. The reported RASPA3 framework–H₂ van der Waals cutoff is 12 Å. `cell_expansion_evidence.csv` gives the cell lengths, and each metal directory contains the original and selected expanded CIFs under `cell_reference/` for direct inspection.

The complete state-by-state CSV is [`../data/KISTI_adsorption_overall_summary_v5.csv`](../data/KISTI_adsorption_overall_summary_v5.csv). It contains all 6,912 intended rows with an explicit `completed` status; 6,822 rows are complete and 90 remain pending in this snapshot. The three raw examples below are a small inspectable subset of that table.

The v5 adsorption table uses the following state grid:

| Temperature (K) | Pressure (bar) |
|---|---|
| 77, 120, 160, 200, 233, 253, 273, 298 | 5, 100 |

In RASPA3, 5 and 100 bar are represented by `500000` and `10000000` Pa. The production settings are 10,000 initialization cycles and 20,000 production cycles, with `PrintEvery = 100`, `UseChargesFrom = "CIF_File"`, `ChargeMethod = "Ewald"`, and `FugacityCoefficient = 1.0`. The pressure/temperature values in the example `simulation.json` are retained exactly; they are not estimates or post-processing values.

To prepare another state, copy the appropriate metal directory, keep the three files together, replace the CIF and framework `Name`, and update `ExternalTemperature` and `ExternalPressure` in `simulation.json`. The small helper script generates a checked JSON from an existing example:

```bash
python gcmc/generate_simulation_json.py \
  --template gcmc/TiCo/simulation.json \
  --framework-name ceq-TiCo-N19-N2-3c-D3h \
  --temperature 120 \
  --pressure-bar 100 \
  --output /path/to/job/simulation.json
```

Cluster-specific PBS files, account names, executable paths, and zero-byte scheduler logs are omitted. The three example result directories do include generated restart and bias files so that the actual completed runs can be inspected. They are examples, not the full set of 6,912 adsorption states; the full-run archive is recorded separately in `raw_manifest/`.
