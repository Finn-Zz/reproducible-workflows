# Topology and coordination-figure screening inputs

These files document the upstream compatibility screen used before structure generation. They are separate from the GCMC force fields and from the QMOF local-geometry filter.

- `3,6-Coordination.xlsx` and `4,6-Coordination.xlsx` contain the retained coordination-figure records.
- `bb_symmetry_topology_compatibility.csv` maps building-block connectivity/local point-group records to compatible topology representations.
- `lah_point_group_summary.csv` records the connection-relevant building-block symmetry checks.
- `nets/` contains the 14 compact `.cgd` net definitions used for the reported topology families.

The screening criterion is the coordination figure at the metal-assigned vertex, with topological site symmetry used as a secondary compatibility descriptor. These tables are provenance inputs; they do not alter the GCMC force-field parameters.
