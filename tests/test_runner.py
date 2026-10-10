"""The parallel runner: memory-bounded sizing, and completeness that is asserted.

⚠ The completeness tests matter more than the sizing ones. A parallel map that
drops a failed job returns a result set that LOOKS complete, and the project has
already been bitten by the same shape of error in a different guise.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import runner as R  # noqa: E402


def _double(x):
    return x * 2


def _boom(x):
    if x == 3:
        raise ValueError("job 3 is broken on purpose")
    return x * 2


def _none(x):
    return None


def test_a_tiny_job_is_core_bound_and_a_huge_one_is_memory_bound():
    """The `reason` string is the point: an over-tight memory estimate should be
    visible as 'MEMORY-bound' rather than showing up as mysterious slowness."""
    small = R.plan_workers(1 << 20, n_jobs=10_000)          # 1 MB a job
    assert small.workers > 1
    assert "CORE-bound" in small.reason, small.reason

    huge = R.plan_workers(200 << 30, n_jobs=10_000)         # 200 GB a job
    assert huge.workers == 1
    assert "MEMORY-bound" in huge.reason, huge.reason


def test_the_worker_count_never_exceeds_the_number_of_jobs():
    """Three jobs must not start forty-eight processes."""
    p = R.plan_workers(1 << 20, n_jobs=3)
    assert p.workers == 3 and "job-count-bound" in p.reason


def test_array_bytes_matches_the_arithmetic_it_claims():
    """215 MB for a chr21 alpha array is quoted all over this repo; check it."""
    assert R.array_bytes(305_886, 88, n_arrays=1) == 305_886 * 88 * 8
    assert abs(R.array_bytes(305_886, 88, n_arrays=1) / 1e6 - 215.3) < 0.5


def test_results_come_back_in_job_order_not_completion_order():
    """Out-of-order results would silently mis-pair every number with its label."""
    jobs = list(range(50))
    assert R.run_jobs(_double, jobs, workers=4) == [j * 2 for j in jobs]


def test_a_single_failed_job_fails_the_whole_run():
    """⚠ THE IMPORTANT ONE. The alternative is a result set that looks complete."""
    with pytest.raises(RuntimeError, match="FAILED"):
        R.run_jobs(_boom, list(range(8)), workers=3)


def test_the_failure_message_names_the_job_so_it_can_be_rerun():
    with pytest.raises(RuntimeError) as e:
        R.run_jobs(_boom, list(range(8)), workers=3)
    msg = str(e.value)
    assert "must not be quoted" in msg
    assert "broken on purpose" in msg, "the underlying error is swallowed"


def test_a_job_returning_nothing_without_raising_is_caught():
    """Defensive, but cheap: a worker that dies oddly leaves a None behind."""
    with pytest.raises(RuntimeError, match="returned nothing"):
        R.run_jobs(_none, list(range(4)), workers=2)


def test_an_empty_job_list_is_refused():
    """Zero jobs completing successfully is the exact ambiguity this guards."""
    with pytest.raises(ValueError, match="no jobs"):
        R.run_jobs(_double, [], workers=2)
