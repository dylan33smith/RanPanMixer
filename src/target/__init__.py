"""Target-side preprocessing: the private input, kept strictly downstream.

Runs per target and consumes the cohort artifacts read-only. Split from `cohort`
so that the cohort stage cannot see a target even by accident -- see
`cohort/__init__.py`.
"""
