# Structure inputs

This directory contains the structure files used by the adsorption analysis and the coordinate inputs used to build the designed frameworks.

- `TiCo/`, `TiMg/`, and `TiNi/` contain the 432 normalized framework CIFs used in the v5 adsorption package, with 144 structures for each heterometallic node family.
- `BBs/` contains the connection-oriented XYZ building blocks and metal-node coordinate inputs used during structure construction. See [`BBs/README.md`](BBs/README.md) for the subdirectory contents and the point-group labels used for the organic building blocks.
- `structure_manifest_432.csv` records the normalized CIF names, relative paths, sizes, and SHA-256 hashes.

The CIFs are the inputs used directly by the GCMC and descriptor workflows. The BBs files are construction provenance and coordinate inputs; they are kept separate from the normalized CIFs so that the origin of each type of structure file remains clear.
