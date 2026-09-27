# Building-block coordinate inputs

This directory contains the connection-oriented XYZ inputs used to represent the linker building blocks and metal nodes in the design library:

- `3-c/` contains the 13 three-connected organic building blocks, labelled by the linker number and connectivity-relevant local point group.
- `4-c/` contains the 13 four-connected organic building blocks, including the rectangular and tetrahedral cases used for the stp-derived and other 4-c topologies.
- `metal/` contains the TiCo, TiMg, TiNi, and Cr3F node coordinate inputs.

These are coordinate inputs/provenance for the generated structures. The normalized CIFs used by the adsorption analysis are in the sibling `../TiCo`, `../TiMg`, and `../TiNi` directories. Cluster job scripts and generated output files are not included here.
