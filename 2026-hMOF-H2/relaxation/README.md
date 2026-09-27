# Structure-relaxation code

This directory contains the two structure-relaxation stages used before the
adsorption calculations:

1. `scripts/forcite_uff_opt.pl` runs Materials Studio Forcite geometry
   optimization with the Universal force field (UFF), ultra-fine quality, and
   cell optimization enabled. It is a MaterialsScript macro and therefore
   must be run inside Materials Studio. Before running it, edit the `$cwd`,
   input, and output directories at the top of the script to match the local
   Windows project location. The comment referring to DREIDING is historical;
   the actual setting in the script is `CurrentForcefield => "Universal"`.
2. The Python scripts in `scripts/` prepare and run the relayed LAMMPS/UFF4MOF
   workflow and restore the bonded topology after cell relaxation:
   `prepare_lammps_inputs_from_ms.py`, `run_lammps_ms_uff4mof.py`,
   `run_lammps_ms_uff4mof_relay_cell.py`, and
   `restore_cell_bonded_cifs.py`.

The LAMMPS workflow first minimizes atomic positions at fixed cell, then
relaxes the cell from the atom-relaxed coordinates with anisotropic
`fix box/relax`, followed by a fixed-cell cleanup. The scripts use the
minimization settings reported in the manuscript/SI.

The Python scripts expect the original working-archive layout under
`Structures_MS`. Set `HMOF_WORKSPACE_ROOT` to the archive root when it is not
the current directory. Set `LAMMPS_INTERFACE_PYTHON` to the Python executable
that provides `lammps_interface`, `LMP_BIN` to the LAMMPS executable, and
optionally `LAMMPS_OMP_THREADS` for the OpenMP thread count. For example:

```bash
export HMOF_WORKSPACE_ROOT=/path/to/hydrogen
export LAMMPS_INTERFACE_PYTHON=/path/to/lammps-interface-python
export LMP_BIN=lmp
python relaxation/scripts/prepare_lammps_inputs_from_ms.py
python relaxation/scripts/run_lammps_ms_uff4mof.py
python relaxation/scripts/run_lammps_ms_uff4mof_relay_cell.py
```

Only the workflow code is included here. The full Forcite/LAMMPS input,
trajectory, log, and restart archives remain in the raw working archive rather
than being duplicated in this repository.
