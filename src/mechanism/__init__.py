"""The mechanism: the cohort baseline, the bounded utility, and the exact sampler.

This is the half of the codebase the guarantee lives in. `cohort/` and `target/`
prepare artifacts; nothing in them knows what `tau` is. Here `tau` is the axis
everything is stamped with.

⚠ THE ASYMMETRY THAT MATTERS. `hmm.py` may not see a target -- `baseline_model`
is a function of `G` and `D` alone (standing constraint 1), and a test asserts
that no function in it takes a target-shaped argument. `utility.py` and `ffbs.py`
MUST see the target: the tilt is the only place it enters. So the import ban that
`cohort/__init__.py` states for a whole subpackage is here a ban on ONE module,
and it is enforced the same way -- by a test that reads the source.

⚠ `toy.py` AND `exact.py` READ NOTHING FROM DISK. They are the falsification
target: a hand-written cohort small enough to enumerate every output path, and a
deliberately naive O(K^T) reference that re-derives the transition formula from
the published Methods rather than importing `hmm.py`. Two implementations that
share code cannot disagree, and a check that cannot fail is decoration -- the
same defect this project documented in PanMixer's own shape assertion.
"""
