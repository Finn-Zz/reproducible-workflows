# v5 merge and validation report

This report records the merge of the verified second-round archive `2round (1).zip` into the v4 adsorption summary.

- The archive registry contains 153 second-round jobs.
- 126 jobs have a finished RASPA output with all final loading fields parseable; 27 remain incomplete.
- Of the 126 complete jobs, 59 were already present in v4 with identical values and 67 were newly applied.
- The resulting v5 table has 6,912 rows for 432 designed structures: 6,822 complete rows and 90 incomplete rows.
- The v5 schema is identical to v4; no values were imputed and no uncertainty column was added.

A job is classified as complete only when the final output contains `Simulation finished!`, the expected temperature and pressure, every required loading field, and internally consistent loading units. Status files and success markers are retained as cross-checks, not used as substitutes for a parseable final output.

The v5 table is the packaged analysis input. The rebuild utility is `scripts/update_summary_v5_from_2round_1.py`; it requires the raw ZIP and the v4 table, which must be downloaded from the external raw-data location listed in `raw_manifest/README.md`.
