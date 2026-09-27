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

The example is the `ceq–N2` framework at 77 K and 5 bar. The CIF and the three JSON files are copied from the same production run directory, so the `Name`, framework filename, and local `ForceField`/`MoleculeDefinition` references agree.

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

Cluster-specific PBS files, account names, executable paths, output directories, restart binaries, and bias-factor files are deliberately omitted from this public input bundle. They are execution-site details or generated state, not force-field definitions. The full run archive is recorded separately in `raw_manifest/`.
