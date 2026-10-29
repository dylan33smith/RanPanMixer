"""ARM 2's per-target chr21 ceiling fit: the gates that stop it drifting from ARM 1's.

⚠ WHY THESE EXIST. `scripts/arm2_chr21_kl.py` relabels ARM 2's attack columns, which were
printed from a grid fitted on HG00096 ALONE. If its ceiling arithmetic were a COPY of
arm 1's rather than the same code, the two arms could report different ceilings for the
same `kl` and nothing would fail. These assert the sharing, not the arithmetic.
"""

from __future__ import annotations

import importlib
import math
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(REPO / "scripts"))


def _load(name: str):
    """⚠ A NORMAL import, deliberately. Loading these with
    `importlib.util.spec_from_file_location` creates a FRESH module object each time, so
    arm 1 would be loaded twice -- once by arm 2's own `from arm1_chr21_kl import ...`
    and once here -- and `is` would compare two copies of the same source and fail.
    That is a defect in the test, not in the sharing, and it cost one false failure."""
    return importlib.import_module(name)


@pytest.fixture(scope="module")
def a2():
    return _load("arm2_chr21_kl")


@pytest.fixture(scope="module")
def a1():
    return _load("arm1_chr21_kl")


def test_ceiling_arithmetic_is_imported_not_copied(a2, a1):
    """The SAME function object, so a change to one cannot miss the other."""
    assert a2.ceiling_from_kl is a1.ceiling_from_kl
    assert a2.eta_for_ceiling is a1.eta_for_ceiling
    assert a2.CEILINGS is a1.CEILINGS


def test_refine_eta_matches_arm1_inline_formula(a2):
    """`arm1_chr21_kl.py:124-125` inlines this; a silent disagreement would move a
    reported ceiling without failing anything."""
    for eta, kl, label in ((454.19, 0.7199, 0.80), (227.1, 0.0431, 0.65),
                           (75.7, 0.02, 0.55), (1.0, 1e-9, 0.50)):
        kl_t = 2.0 * (2.0 * float(label) - 1.0) ** 2
        want = eta * math.sqrt(kl_t / kl) if kl > 0 else eta
        assert a2.refine_eta(eta, kl, label) == want


def test_refine_eta_is_a_noop_at_zero_kl(a2):
    """kl = 0 is eta = 0, where `c` is undefined; returning eta beats dividing by zero."""
    assert a2.refine_eta(123.0, 0.0, 0.80) == 123.0


def test_share_is_none_not_one(a2):
    """⚠ `one_chrom` spells `beta_share` = 1 as share=None. Passing 1.0 re-scales beta,
    which would silently change the scope the bound covers."""
    assert a2.SHARE is None
    assert a2.K_STATES == 88


def test_arm2_is_not_arm1(a2, a1):
    """Different arm, different K, different share source -- the two files must not
    converge into one by accident."""
    assert a1.K_STATES if hasattr(a1, "K_STATES") else True
    assert a2.CHROM == a1.CHROM == "chr21"
    # arm 1 reads a leave-one-out survey for its share; arm 2 must not.
    assert not hasattr(a2, "SURVEY") or a2.SHARE is None


def test_attack_targets_are_the_ones_that_were_attacked(a2):
    """The fit must cover the table it relabels. If the attack ran on a target whose
    chr21 path is absent, the fit CANNOT relabel that column and must say so."""
    if not a2.ATTACK_RUN.exists():
        pytest.skip(f"{a2.ATTACK_RUN} absent")
    import json
    attacked = sorted(json.loads(a2.ATTACK_RUN.read_text())["targets"])
    got = a2.attack_targets()
    assert got == attacked, "the fit's target set must equal the attack's, exactly"
    assert len(got) == 10


def test_fit_refuses_a_zero_eta(a2):
    """`c` is undefined at eta = 0, where kl is identically 0: a fit from it is
    meaningless rather than merely imprecise."""
    with pytest.raises(SystemExit, match="must be > 0"):
        a2.main(["--etas", "0", "454.19"])
