"""Can `scripts/exact_utility.py`'s reporting path report something that is not there?

⚠⚠ WHY THIS EXISTS. The genome-wide utility figure is a FAN-OUT: 44 exact passes, one per
(chromosome, `tau`). CLAUDE.md's rule is that any fan-out must assert its own completeness
and fail loudly, because a script that averages whatever rows it was handed cannot tell
"found nothing at this `tau`" from "ran nothing at this `tau`" -- and `np.mean([])` returns
`nan` with a RuntimeWarning, not an error, so a silently-truncated run would emit a
well-formed table with `nan` cells rather than refusing to report.

⚠ The second gate here is the RELABEL. The ratio row was emitted as `utility_retained` on
2026-10-12, which `docs/terms.md` defines as the SAMPLED `u_path` average -- a level, not a
ratio. `--from-json` exists so a relabel never requires recomputing 44 passes, which makes
it a correctness-preserving path that must be proven to preserve the VALUES.
"""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

_spec = importlib.util.spec_from_file_location(
    "exact_utility", REPO / "scripts" / "exact_utility.py")
eu = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(eu)


def _rows():
    """The measured 2026-10-12 genome-wide result, as the saved JSON holds it."""
    return [
        {"target": "HG00096", "tau": 0.0, "e_q_u": 0.7089676459975968,
         "u_ceiling": 0.8257087407736636, "e_q_u_over_ceiling": 0.858617101877008,
         "chroms": 22},
        {"target": "HG00096", "tau": 0.99, "e_q_u": 0.7089677294498828,
         "u_ceiling": 0.8257087407736636, "e_q_u_over_ceiling": 0.8586172029444692,
         "chroms": 22},
    ]


def test_report_refuses_a_tau_with_no_rows(capsys):
    """THE GATE. A `tau` present in the axis but absent from the results must RAISE.

    This is the mutation: ask for a three-`tau` report from two-`tau` data, exactly what a
    run killed between its second and third pass would leave behind.
    """
    with pytest.raises(AssertionError, match="no saved rows at tau"):
        eu._report(_rows(), [0.0, 0.5, 0.99], "GENOME-WIDE")


def test_report_emits_every_tau_it_was_given(capsys):
    """The CONTROL for the gate above: the same call with matched rows must not raise."""
    eu._report(_rows(), [0.0, 0.99], "GENOME-WIDE")
    out = capsys.readouterr().out
    assert "| `e_q_u` * |" in out
    assert "| `e_q_u_over_ceiling` |" in out
    # the row label `utility_retained` is a DIFFERENT quantity (the sampled u_path) and
    # must never reappear on this table -- see docs/terms.md
    assert "utility_retained" not in out


def test_the_movement_is_reported_not_just_the_level(capsys):
    """The level is reached at `tau` = 0, so the MOVEMENT is the only reportable figure.

    A table that printed only levels would invite exactly the free-lunch reading
    `target_fidelity` exists to prevent, so the movement row is part of the contract.
    """
    eu._report(_rows(), [0.0, 0.99], "GENOME-WIDE")
    out = capsys.readouterr().out
    assert "THE MOVEMENT" in out
    assert "+8.345e-08" in out


def test_from_json_preserves_the_saved_values(tmp_path, capsys):
    """`--from-json` must re-report, never recompute -- and never perturb a value."""
    f = tmp_path / "saved.json"
    f.write_text(json.dumps({"scope": "GENOME-WIDE", "results": _rows()}))
    assert eu.main(["--from-json", str(f)]) == 0
    out = capsys.readouterr().out
    assert "NO marginals computed" in out
    for cell in ("0.708967646", "0.708967729", "0.825708741", "0.858617102"):
        assert cell in out, f"{cell} missing: the re-report changed a saved value"


def test_the_shipped_json_carries_the_relabelled_key():
    """The saved artifact must use the term `docs/terms.md` actually defines."""
    f = REPO / "logs" / "B_exact_utility.json"
    if not f.exists():
        pytest.skip("logs/B_exact_utility.json not present")
    rows = json.loads(f.read_text())["results"]
    assert rows, "the shipped result file is empty"
    for r in rows:
        assert "e_q_u_over_ceiling" in r, "relabel did not reach the saved JSON"
        assert "utility_retained" not in r, (
            "the saved JSON still calls the exact RATIO `utility_retained`, which "
            "docs/terms.md defines as the SAMPLED u_path average")
