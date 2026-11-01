"""ARM 2's chr21 four-condition envelope driver: the gates that keep it honest.

⚠ WHY THESE EXIST. `scripts/arm2_chr21_ceiling.py` produces the statement that every
Arm 2 chr21 column label bounds the real-genome decoys, which is the claim the attack
tables are read against. Two things could make it return plausible, wrong numbers
without failing: it could DRIFT from the attack's own decoy choice, and it could run
on the wrong `beta_v`. The first is closed by importing rather than copying, the
second by re-measuring `D1_extremal` against the grid its own `eta` was solved for.
These assert both, plus that the instrument check CAN fail -- a gate that cannot is
decoration.

⚠ Nothing here calls `kl_stable`: one call is ~17 s at chr21 scale and the suite is
not the place for it. The one test that needs a divergence substitutes a fake, which
is what lets it assert the GATE's logic rather than the HMM's arithmetic.
"""

from __future__ import annotations

import importlib
import json
import re
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(REPO / "scripts"))

SRC = REPO / "scripts" / "arm2_chr21_ceiling.py"


def _load(name: str):
    """A normal import -- see tests/test_arm2_chr21_kl.py for why not spec_from_file."""
    return importlib.import_module(name)


@pytest.fixture(scope="module")
def mod():
    return _load("arm2_chr21_ceiling")


@pytest.fixture
def verified(mod):
    """⚠ `logs/` IS GITIGNORED ON EVERY BRANCH, so the verified grid is not in git and
    a fresh checkout does not have it. The three refusal tests below read that file to
    build their tampered copies, so they SKIP rather than error when it is absent --
    otherwise this file would fail on the release branch for a reason that has nothing
    to do with the code it gates."""
    if not mod.VERIFIED.exists():
        pytest.skip(f"{mod.VERIFIED.relative_to(REPO)} absent (logs/ is gitignored)")
    return json.loads(mod.VERIFIED.read_text())


@pytest.fixture(scope="module")
def atk():
    return _load("arm2_attack")


def test_decoy_machinery_is_imported_not_copied(mod, atk):
    """The SAME objects, so the measured pairs cannot drift from the attacked pairs."""
    assert mod.pick_decoys is atk.pick_decoys
    assert mod.load_state is atk.load_state
    assert mod.DECOY_SEED is atk.DECOY_SEED
    assert mod.CONDS is atk.CONDS
    assert mod.ceiling_from_kl is _load("arm1_chr21_kl").ceiling_from_kl


def test_no_local_redefinition_of_the_imported_helpers():
    """⚠ `is` above passes even beside a shadowed copy defined LATER in the file.

    A second definition further down would rebind the name at import time and the
    identity test would then fail -- but only if the shadow came before this module
    was read. Asserting on the source closes the case where someone adds a local
    helper under a different name and calls it instead.
    """
    src = SRC.read_text()
    for name in ("pick_decoys", "load_state", "concordance", "ceiling_from_kl",
                 "kl_stable", "extremal"):
        assert not re.search(rf"^def {name}\b", src, re.M), (
            f"{name} is defined locally in {SRC.name}; it must be imported, or this "
            f"measurement and the attack can disagree and neither will fail.")


def test_grid_is_read_from_the_verified_file_not_resolved(mod):
    """The `eta` must come from the file the attack reads its own bound from."""
    assert mod.VERIFIED.name == "arm2_chr21_kl_verified.json"
    src = SRC.read_text()
    assert "eta_for_ceiling" not in src and "refine_eta" not in src, (
        "this driver must not re-solve the grid: a grid re-solved here could not be "
        "compared to the attack cells it is meant to bound.")


def test_load_grid_refuses_a_mismatched_weighting(mod, verified, tmp_path, monkeypatch):
    """Two weightings' `eta` are different quantities wearing one name."""
    d = dict(verified)
    d["weights"] = "inv_support" if mod.WEIGHTS != "inv_support" else "uniform"
    f = tmp_path / "bad_weights.json"
    f.write_text(json.dumps(d))
    monkeypatch.setattr(mod, "VERIFIED", f)
    with pytest.raises(SystemExit, match="weights"):
        mod.load_grid()


def test_load_grid_refuses_a_grid_that_is_not_chr21_scope(mod, verified, tmp_path, monkeypatch):
    """`load_state` builds `beta_v` at scope 'chain', share 1. A grid solved at any
    other share is a different quantity and must not silently label these cells."""
    d = dict(verified)
    d["beta_share"] = 0.0150227            # chr21's genome-wide share
    f = tmp_path / "bad_share.json"
    f.write_text(json.dumps(d))
    monkeypatch.setattr(mod, "VERIFIED", f)
    with pytest.raises(SystemExit, match="beta_share"):
        mod.load_grid()


def test_load_grid_refuses_a_fit_with_no_refined_block(mod, verified, tmp_path, monkeypatch):
    """A grid quoted from a fit is a derived number wearing a measured number's clothes."""
    d = dict(verified)
    for rec in d["per_target"].values():
        rec.pop("refined", None)
    f = tmp_path / "fit_only.json"
    f.write_text(json.dumps(d))
    monkeypatch.setattr(mod, "VERIFIED", f)
    with pytest.raises(SystemExit, match="FIT"):
        mod.load_grid()


def test_instrument_check_fires_when_d1_disagrees(mod, monkeypatch):
    """⚠⚠ THE GATE THAT MAKES THE REST TRUSTWORTHY, AND IT MUST BE ABLE TO FAIL.

    MEASURED 2026-10-31 on the real instrument: a `beta_v` built at chromosome-share
    0.5 instead of 1 returns D1 = 0.179957 where the verified grid says 0.720000 --
    and the SAME wrong run returns D4 = 0.066949, an entirely plausible number under
    an envelope that still appears to hold. So the wrong wiring is not visible in the
    numbers the run exists to report, only in this comparison. Here the divergence is
    faked so the test is fast; what is asserted is that `one_cell` compares D1 against
    the expected value and reports the verdict, not that the HMM is right.
    """
    paths = {"T": 0, **{c: 0 for c in mod.REAL}}
    monkeypatch.setattr(mod, "load_state",
                        lambda: {"paths": paths, "model": None, "hap": None,
                                 "util": type("U", (), {"beta": None})()})
    monkeypatch.setattr(mod, "extremal", lambda p: p)
    monkeypatch.setattr(mod, "kl_stable", lambda *a, **k: {"kl": 0.5})
    ids = {c: c for c in mod.REAL}

    ok = mod.one_cell(("T", "0.80", 400.0, 0.5, ids))
    assert "error" not in ok and ok["d1_ok"] is True, ok

    bad = mod.one_cell(("T", "0.80", 400.0, 0.72, ids))
    assert "error" not in bad, bad
    assert bad["d1_ok"] is False, (
        "a D1 divergence 0.22 away from the value its own eta was solved for was "
        "accepted: the instrument check is decoration and every real-decoy number "
        "this driver reports is unverified.")
    assert bad["d1_err"] == pytest.approx(0.5 - 0.72)


def test_the_atol_is_tight_enough_to_catch_a_real_miswiring(mod):
    """0.54 was the measured error from a wrong share; the tolerance must be far below
    it, and far below any plausible accumulation-order difference."""
    assert mod.D1_ATOL <= 1e-6
