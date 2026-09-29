"""Cohort-side preprocessing: everything built from the graph and the cohort alone.

⚠ NOTHING IN THIS SUBPACKAGE MAY TAKE, IMPORT, OR NAME A TARGET.

Standing constraint 1 says the baseline model and the output support are
functions of G and D only, never the target. `docs/plan.md` open question 4
recorded that this was "enforced by discipline; it should be enforced by
construction". This subpackage is that construction: no function here has a
target parameter, and no module here imports `target`. A test asserts the import
ban directly, so the rule fails loudly at CI time rather than silently at
release time.
"""
