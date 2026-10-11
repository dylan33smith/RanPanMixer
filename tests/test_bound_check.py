"""Does the `p_succ_bound` check fire on a REAL excess and stay quiet on noise?

⚠⚠ WHY BOTH DIRECTIONS MATTER. The first version compared a point estimate to the
bound with no uncertainty. At `tau` = 0.05 the bound sits 0.025 above chance while the
balanced accuracy's SE is 0.0223, so a true accuracy of 0.5 exceeds the bound about 13%
of the time PER CELL -- about 2.3 false alarms over 18 cells. It duly fired on the
2026-10-11 re-run (D3_random, `tau` = 0.05, 0.5400 vs 0.5250, z = +0.67) and called it
a falsification of the implementation.

⚠ But a check that merely demands significance could be made quiet by weakening it
until nothing fires. So these tests pin BOTH: a noise-level excess must NOT fail the
run, and a genuine excess MUST.
"""

from __future__ import annotations

import json
import math
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
SCRIPT = REPO / "scripts" / "arm2_attack.py"
sys.path.insert(0, str(REPO / "scripts"))


def _trials(tau: float, n: int, frac_right: float) -> list:
    """n trials at one `tau` whose balanced accuracy is `frac_right` by construction.

    Lambda > 0 under p and < 0 under q is a correct call in each arm, so setting that
    fraction in both arms sets the balanced accuracy exactly.
    """
    k = int(round(frac_right * n))
    out = []
    for i in range(n):
        right = i < k
        cell = {"decoy": "x", "u_q": 0.1, "concordance": None,
                "lambda": (1.0 if right else -1.0),
                "lambda_under_q": (-1.0 if right else 1.0)}
        out.append({"target": f"T{i % 10}", "tau": tau, "seed": i,
                    "log_z_p": 0.15, "u_p": 0.5,
                    "conditions": {"D1_extremal": dict(cell),
                                   "D2_nearest": dict(cell),
                                   "D3_random": dict(cell)}})
    return out


def _run(tmp_path: Path, trials: list) -> tuple[int, str]:
    f = tmp_path / "trials.json"
    f.write_text(json.dumps({"trials": trials}))
    r = subprocess.run([sys.executable, str(SCRIPT), "--from-json", str(f)],
                       capture_output=True, text=True, cwd=REPO)
    return r.returncode, r.stdout + r.stderr


def test_a_genuine_excess_fails_the_run(tmp_path: Path) -> None:
    """⚠ THE POWER TEST. 95% accuracy at `tau` = 0.05, where the bound is 0.525, is an
    excess of 0.425 at an SE of about 0.015 -- around 28 SE. If this does not fail, the
    check has been weakened into decoration."""
    rc, out = _run(tmp_path, _trials(0.05, 500, 0.95))
    assert rc == 1, f"a 28-SE excess did not fail the run.\n{out[-2000:]}"
    assert "falsifies the IMPLEMENTATION" in out


def test_a_noise_level_excess_does_not_fail_the_run(tmp_path: Path) -> None:
    """⚠ THE FALSE-ALARM TEST, which is the one that bit. 54% at `tau` = 0.05 is +0.015
    over the bound at an SE of ~0.022 -- under 1 SE. It must be REPORTED and must not
    fail."""
    rc, out = _run(tmp_path, _trials(0.05, 500, 0.54))
    assert rc == 0, f"a sub-1-SE excess failed the run.\n{out[-2000:]}"
    assert "POINT ESTIMATE" in out and "not claimed" in out.lower()


def test_a_cell_at_chance_is_silent(tmp_path: Path) -> None:
    """The control: nothing above the bound means neither a violation nor a note."""
    rc, out = _run(tmp_path, _trials(0.25, 500, 0.50))
    assert rc == 0
    assert "POINT ESTIMATE" not in out


def test_the_normal_quantile_round_trips_against_erf() -> None:
    """⚠ The threshold is computed without scipy, so the implementation is checked
    rather than trusted -- a wrong quantile silently changes how loud the gate is."""
    import arm2_attack as M
    for alpha, expect in ((0.05, 1.6449), (0.025, 1.9600), (0.01, 2.3263),
                          (0.001, 3.0902)):
        z = M._z_for_one_sided(alpha)
        assert abs(z - expect) < 1e-3, f"alpha={alpha}: got {z}, expected ~{expect}"
        # and the round trip: P(Z > z) must be alpha
        back = 0.5 * (1.0 - math.erf(z / math.sqrt(2.0)))
        assert abs(back - alpha) < 1e-9


def test_the_threshold_tightens_with_more_cells() -> None:
    """Bonferroni: more cells examined must mean a stricter per-cell threshold, or the
    multiple-comparison correction is not happening."""
    import arm2_attack as M
    assert M._z_for_one_sided(0.05 / 1) < M._z_for_one_sided(0.05 / 21)
