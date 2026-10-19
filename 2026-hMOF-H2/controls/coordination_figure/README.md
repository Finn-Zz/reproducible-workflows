# Coordination-figure control-net structures

## Why this control set was generated

The main construction workflow selects nets by comparing the coordination figure of each metal-assigned 6-connected (6-c) vertex with the trigonal-prismatic arrangement of the six connection directions in the MIL-100-derived metal building block. A control set was designed to test the significance of that geometry-matching choice. It uses the same metal-node families and the same 13 tritopic or 13 tetratopic organic building blocks, but assigns them to bipartite (3,6)-c or (4,6)-c RCSR nets whose metal-assigned 6-c vertices do not consistently have the required trigonal-prismatic coordination figure. The design therefore probes what happens when vertex connectivity is suitable but the local connection geometry is mismatched. It does not by itself establish that any resulting framework can or cannot be synthesized.

The eight (3,6)-c control nets are `pyr`, `spn`, `ant`, `rtl`, `tsx`, `cer`, `laf`, and `hxh`; the four (4,6)-c nets are `soc`, `iac`, `alw`, and `she`. Each metal-node family has 8 × 13 + 4 × 13 = 156 planned net–linker combinations. The four families (`Cr3F`, `TiCo2`, `TiMg2`, and `TiNi2`) give 624 planned structures in this data export.

## Files and status

- `optimized_CIFs/` contains **448 completed UFF4MOF-optimized CIFs** from the `cifs-1014` workflow, grouped by metal node: Cr3F 109, TiCo2 109, TiMg2 119, and TiNi2 111. Original CIF basenames are retained.
- `MOF_structure_checks.csv` is the full, unmodified 624-row source check table. It records optimization status and the available MOFClassifier, MOFChecker, Chen–Manz, and QMOF-filter results. Its `final_cif_path` and report paths reflect the original local Windows workspace and are provenance fields, not portable links.
- `MOF_structure_checks_simplified.csv` is the unmodified seven-column view of the same 624 designs.
- `optimized_CIF_manifest.csv` lists the 448 included CIFs, their relative paths, byte counts, and SHA-256 checksums.

The other **176 designs do not have completed UFF4MOF CIFs** in the source folder and are represented by status rows in the CSV, not by structure files. The full design includes 12 control nets; the completed CIF collection contains 11 because no `hxh` optimization completed. A completed optimization is not equivalent to passing every structural-plausibility check; consult the CSV for each structure's results.

## Relation to the manuscript comparison

The L2 manuscript/SI describes a preliminary **Cr3F-only** comparison of 156 control-net combinations against 144 separately reconstructed structures on geometry-matched nets. Those L2 comparison structures were described as relaxed directly with UFF4MOF, without Forcite pre-relaxation. This uploaded `cifs-1014` collection is a later, four-metal-node workflow: its `*_UFF4MOF` CIFs were generated from `*_UFF` (Materials Studio Forcite/UFF) intermediates. It should therefore not be treated as the raw input/output of the L2 head-to-head comparison, or as proof of the numerical outcomes reported there. The 144 matched Cr3F comparators are not included in this directory.

No H2-adsorption results are provided for these control-net structures. The adsorption analysis in the main repository concerns the separately identified 432 heterometallic structures on the geometry-matched nets.
