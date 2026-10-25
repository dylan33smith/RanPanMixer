"""GATES on the 2026-10-25 failure: `exact_kl_stable.py` ran 22 chromosomes for 3.4
hours in `--etas` mode, then destroyed all of it in the aggregation step.

WHAT WENT WRONG, in one sentence: the knob grid had two names and the knob key had two
spellings, so the aggregator walked a seven-point `tau` grid the workers had never
evaluated and asked for it under a spelling the workers had never written.

  main:    kl = math.fsum(per[c]["kl"][t][str(tau)]["kl"] for c in AUTOSOMES)
  worker:  out[t] = {f"{v:g}": kl_stable(...)}
  -> KeyError: '0.0'          (str(0.0) == '0.0', f"{0.0:g}" == '0')

A parallel audit found 21 further confirmed defects behind that one. The ones that
would have destroyed the RE-run, in order of when they would have fired:
  1. `p_succ_bound(eta)` raises ValueError unless 0 <= tau < 1  -> crash at row 2
  2. `_report(rows, a.taus, ...)` -> `max()` over an empty generator -> crash
  3. the JSON write sat AFTER `_report`, so either crash loses everything again
Each of those is gated below. ⚠ The point of these tests is that they fail on the
ORIGINAL code; a gate written beside its own fix proves nothing, so the ones that can
be mutation-tested say so in their docstring.
"""

from __future__ import annotations

import ast
import json
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
SCRIPT = REPO / "scripts" / "exact_kl_stable.py"
sys.path.insert(0, str(REPO / "scripts"))
sys.path.insert(0, str(REPO / "src"))


# --------------------------------------------------------------------------- key
def test_knob_key_round_trips_and_fixes_the_exact_regression():
    """THE defect. `eta` = 0.0 must key identically from both sides."""
    from exact_kl_stable import knob_key

    # the exact pair that collided
    assert knob_key(0.0) == "0.0", "eta 0.0 must not key as '0'"
    assert f"{0.0:g}" == "0" != knob_key(0.0), (
        "the OLD producer spelling and the new one must differ -- if they agree, this "
        "test is not pinning the regression")

    # every grid this project actually runs, both knobs
    for v in (0.0, 0.05, 0.10, 0.25, 0.50, 0.90, 0.99,
              0.999, 0.9999, 0.99999,
              0.0, 34.33, 68.66, 173.3, 387.39, 774.78,
              250.0, 500.0, 1262.0, 2821.0, 5642.0, 6.6566e7):
        assert float(knob_key(v)) == float(v), f"{v} does not round-trip"


def test_knob_key_is_not_lossy_where_g_would_be():
    """`f"{v:g}"` keeps 6 significant figures, so two near grid points COLLIDE.

    A collision silently overwrites one cell with another's value -- a wrong number,
    not an error. This is the defect class, not just the 0.0 instance.
    """
    from exact_kl_stable import knob_key

    a, b = 1234567.0, 1234568.0          # agree to 6 s.f.
    assert f"{a:g}" == f"{b:g}", "pick a closer pair; these do not collide under :g"
    assert knob_key(a) != knob_key(b), (
        "knob_key collapsed two distinct grid points into one key; one cell would "
        "silently overwrite the other")


def test_both_producers_and_the_consumer_use_the_one_helper():
    """Static drift check: a second spelling ANYWHERE is the bug coming back.

    ⚠ The runtime tests above cannot see a future call site that re-introduces
    `f"{v:g}"`. This walks the source instead.
    """
    tree = ast.parse(SCRIPT.read_text())
    bad = []
    for node in ast.walk(tree):
        # a dict comprehension/literal key built by formatting the knob
        if isinstance(node, ast.JoinedStr):
            txt = ast.unparse(node)
            if ":g}" in txt and "knob" not in txt and "=" not in txt:
                bad.append(f"line {node.lineno}: {txt}")
    assert not bad, (
        f"f-string knob formatting used where knob_key() should be: {bad}. Result "
        f"keys must come from the ONE helper so producer and consumer cannot drift.")

    src = SCRIPT.read_text()
    assert src.count("knob_key(v)") >= 2, (
        "both one_chrom branches (arm1 and arm2) must key with knob_key")
    assert 'per[c]["kl"][t][knob_key(v)]' in src, (
        "the aggregator must look up with knob_key, not str() or an f-string")


# ------------------------------------------------------------------- ordering
def test_the_json_is_written_before_the_report():
    """⚠⚠ THE LOSS MODE. A formatter must never be able to destroy 3.4h of compute.

    On 2026-10-25 all 22 chromosomes finished, `_report` raised, and because the write
    was downstream of it the entire run was lost.
    """
    src = SCRIPT.read_text()
    write = src.index("a.json.write_text")
    report = src.index("_report(rows")
    assert write < report, (
        "the JSON write must come BEFORE the _report call; otherwise a reporting bug "
        "destroys a completed computation, which is exactly what happened")


def test_the_json_carries_the_per_chromosome_results():
    """The expensive part must be recoverable from disk whatever happens next."""
    src = SCRIPT.read_text()
    assert '"per_chrom"' in src, (
        "the JSON must include per_chrom; 22 chromosome-hours held only in memory is "
        "what made a one-line bug cost 3.4 hours")


# ------------------------------------------------------------------ provenance
def test_the_json_stamps_arm_weights_and_knob_rather_than_hardcoding_arm2():
    src = SCRIPT.read_text()
    assert '"arm": "arm2"' not in src, (
        'the JSON hardcoded "arm": "arm2" even under --arm arm1 (which is K = 86), so '
        'an Arm 1 run wrote a file labelled Arm 2')
    assert '"k_states": 88,' not in src, "k_states was hardcoded to 88"
    for key in ('"arm": a.arm', '"weights": a.weights', '"knob": knob'):
        assert key in src, f"the JSON must record {key}"


def test_report_no_longer_hardcodes_arm_2_and_k_88():
    src = SCRIPT.read_text()
    assert "**ARM 2**, K = 88" not in src, (
        "_report printed a fixed 'ARM 2, K = 88' provenance line regardless of --arm")


# ----------------------------------------------------------- crash sites gated
def test_report_names_a_grid_mismatch_instead_of_raising_empty_max():
    """Crash site 2, with the ORIGINAL symptom asserted against.

    `max()` over an empty generator raises a bare `ValueError: max() iterable argument
    is empty`, which says nothing about why. The whole run was already gone by then.
    """
    from exact_kl_stable import _report

    rows = [{"target": "T", "knob": "eta_tau", "knob_value": 0.0, "tau": 0.0,
             "tau_equiv": 0.0, "kl_exact": 0.0, "kl_raw": 0.0, "tv_pinsker": 0.0,
             "acc_ceiling": 0.5, "p_succ_bound": 0.5, "bound_vacuous": False,
             "was_negative": False, "negative_chroms": []}]
    # grid contains a point no row has -- precisely the tau-grid-vs-eta-rows mismatch
    with pytest.raises(AssertionError, match="does not match the grid"):
        _report(rows, [0.0, 0.05], 1, 100, knob="eta_tau", weights="uniform",
                arm="arm1")


def test_duplicate_grid_points_are_refused_up_front():
    """Two grid points sharing one key would overwrite each other silently."""
    r = subprocess.run(
        [sys.executable, str(SCRIPT), "--arm", "arm1", "--targets", "1",
         "--etas", "250", "250", "--workers", "2"],
        capture_output=True, text=True, timeout=600, cwd=REPO)
    assert r.returncode != 0
    assert "duplicate knob values" in (r.stdout + r.stderr), (
        f"a repeated grid point was accepted: {(r.stdout + r.stderr)[-600:]}")


def test_validate_refuses_the_configurations_it_cannot_anchor():
    """Crash-free but worse: a gate that passes on a config it never ran.

    `validate()` reproduces logged chr21 / ARM 2 / inv_support figures and then prints
    'its genome-wide output may be reported'. Under --arm arm1 or --etas it was
    vouching for a configuration it had not touched.
    """
    for extra in (["--arm", "arm1"], ["--etas", "250"],
                  ["--weights", "uniform"]):
        r = subprocess.run([sys.executable, str(SCRIPT), "--validate"] + extra,
                           capture_output=True, text=True, timeout=600, cwd=REPO)
        out = r.stdout + r.stderr
        assert r.returncode != 0, f"--validate {extra} should refuse, got 0:\n{out[-800:]}"
        assert "no anchor on disk" in out, f"wrong refusal for {extra}: {out[-500:]}"


def test_p_succ_bound_is_not_handed_an_eta():
    """Crash site 1. `p_succ_bound` raises unless 0 <= tau < 1; an eta is not a tau."""
    from mechanism import calibration as C

    with pytest.raises(ValueError):
        C.p_succ_bound(250.0)           # the call the old code made

    src = SCRIPT.read_text()
    assert "C.p_succ_bound(tau_eq)" in src, (
        "the bound must be computed from tanh(eta), not from the eta itself")
    assert "tau_eq >= 1.0" in src, (
        "tanh(eta) saturates to 1.0 past eta ~ 19, where p_succ_bound would still "
        "raise -- the vacuous case must be reported, not crashed on")


def test_negative_kl_is_checked_per_chromosome_not_only_on_the_sum():
    """Clamping the 22-term sum lets one negative chromosome hide inside a positive
    total -- a wrong number that prints as a clean one."""
    src = SCRIPT.read_text()
    assert "bad_parts" in src and "negative_chroms" in src, (
        "the sign check must see the per-chromosome parts and name the offenders")


# ------------------------------------------------------------------- resume
def test_checkpoint_refuses_a_foreign_configuration():
    """A cached chromosome from another arm/weighting/grid is a silent wrong answer,
    which is worse than recomputing it. The guard must compare all four fields."""
    src = SCRIPT.read_text()
    for field in ('"_arm"', '"_weights"', '"_grid"', '"_targets"'):
        assert field in src, (
            f"the checkpoint record must carry {field} so --resume cannot reuse a "
            f"chromosome computed under a different configuration")
    assert 'r.get("_arm") == a.arm' in src
