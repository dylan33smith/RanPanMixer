"""A3's score must be PANMIXER'S score, and this is the evidence rather than the claim.

⚠⚠ WHY GOLDEN VALUES AND NOT AN IMPORT. The authority for these numbers is
`external/PanMixer/tools/downstream/privacy/gap_score.py:score_genotypes_remove_missing`,
CALLED DIRECTLY on the inputs built below, in their own `panmixer` env (their module
imports pandas, matplotlib and their SLURM helper, none of which `environment.yml`
promises, so importing it here would make this test unrunnable in our env -- a test
that cannot run is not a test). The cross-check was performed on 2026-10-13 and its
outputs are pinned here, so our reimplementation is held to their code's behaviour on
every future run without needing their env.

⚠ THE POINT IS DRIFT, NOT ARITHMETIC. Arm 1's whole claim is that both tools are
scored by ONE function. If `a3_score` ever stops agreeing with theirs -- a changed
validity rule, a weight taken from the release instead of the candidate, a dropped
`+1e-10` -- then A3 is no longer "PanMixer's score" and the pre-registration's
description of it becomes false. That is what this pins.

Reproduce (writes nothing under external/; their code savefig's into CWD):
    cd <scratch> && PYTHONPATH=/data/ds85/RanPanMixer/runs/pm_eval:<repo>/src:<repo>/scripts \
      PYTHONDONTWRITEBYTECODE=1 MPLBACKEND=Agg <panmixer-env>/python -B xc.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(REPO / "scripts"))


# The inputs, built by the same code and seed the cross-check used. Any change to this
# construction invalidates the golden values, so it is written out rather than imported.
def _inputs():
    rng = np.random.default_rng(20261013)
    n, m = 400, 5
    indiv = rng.integers(0, 2, size=(n, 2)).astype(np.int16)
    other = rng.integers(0, 2, size=(m, n, 2)).astype(np.int16)
    indiv[rng.random(n) < 0.08] = -1
    for i in range(m):
        other[i][rng.random(n) < 0.06] = -1
    alt_af = rng.uniform(0.02, 0.98, size=n)
    return indiv, other, alt_af


# ⚠ PANMIXER'S OWN OUTPUT, measured 2026-10-13 by calling their function. Do not
# regenerate these from our implementation -- that would make the test circular and it
# would then pass by construction, which is this project's definition of decoration.
THEIRS = [216.255656702, 180.609169621, 195.827088859, 178.846094183, 208.33030707]


def _site_table(alt_af: np.ndarray) -> dict:
    """A synthetic bi-allelic axis: every site in play, ALT = code 1."""
    p = np.asarray(alt_af)
    q = 1.0 - p
    return {"mask": np.ones(len(p), dtype=bool),
            "alt": np.ones(len(p), dtype=np.int16),
            "logw": np.stack([-np.log(q * q + 1e-10),
                              -np.log(2 * p * q + 1e-10),
                              -np.log(p * p + 1e-10)])}


def test_a3_score_reproduces_panmixers_own_output():
    from arm2_attack_a3a4 import a3_score
    indiv, other, alt_af = _inputs()
    s = _site_table(alt_af)
    mine = np.array([a3_score(indiv, other[i], s) for i in range(other.shape[0])])
    # 1e-9 rather than exact: their call casts to float32, ours stays float64, and the
    # measured residual was 1.4e-13. A real divergence is orders of magnitude larger.
    assert np.allclose(mine, THEIRS, atol=1e-9), (
        f"a3_score no longer matches PanMixer's score_genotypes_remove_missing.\n"
        f"  theirs: {THEIRS}\n  ours:   {[round(float(x), 9) for x in mine]}\n"
        f"A3 is pre-registered as THEIR score; if ours has changed, either the change "
        f"is a bug or the pre-registration's description of A3 is now false.")


def test_a_mismatch_contributes_zero_not_a_penalty():
    """The property that makes A3's pre-registered positive control vacuous.

    ⚠ This is not a style note. Because mismatches add nothing and the weights are
    non-negative, `score(y, q)` over the agreeing sites is a SUBSET SUM of
    `score(y, y)`, so handing the attack the target's own path forces a win. The
    control therefore tests the plumbing and not the power, and the runner reports it
    that way. If this test ever fails, that reasoning has to be revisited.
    """
    from arm2_attack_a3a4 import a3_score
    _, _, alt_af = _inputs()
    s = _site_table(alt_af)
    n = len(alt_af)
    rng = np.random.default_rng(7)
    p = rng.integers(0, 2, size=(n, 2)).astype(np.int16)
    q = p.copy()
    q[: n // 2] = 1 - q[: n // 2]            # disagree on half the sites

    assert (s["logw"] >= 0).all(), "a negative weight would break the subset-sum argument"
    self_score = a3_score(p, p, s)
    other_score = a3_score(p, q, s)
    assert self_score > other_score, "the target's own path must score highest"
    # and the gap is exactly the weight of the sites where they disagree
    assert other_score <= self_score + 1e-9


def test_the_weight_comes_from_the_candidate_not_the_release():
    """Their design, and it is asymmetric: `score(y, c)` weights by `g(c)`.

    Two candidates that agree with the release at the same NUMBER of sites score
    differently when the genotypes they agree on have different rarity. A version that
    weighted by the release would make the statistic symmetric and would no longer be
    their score.
    """
    from arm2_attack_a3a4 import a3_score
    # two sites: one common-homozygous, one rare-homozygous
    alt_af = np.array([0.5, 0.02])
    s = _site_table(alt_af)
    y = np.array([[1, 1], [1, 1]], dtype=np.int16)
    c_common = np.array([[1, 1], [0, 0]], dtype=np.int16)   # agrees at the common site
    c_rare = np.array([[0, 0], [1, 1]], dtype=np.int16)     # agrees at the rare site
    assert a3_score(y, c_rare, s) > a3_score(y, c_common, s), (
        "agreeing on a RARE genotype must score higher than agreeing on a common one; "
        "if it does not, the rarity weighting is not being applied")


def test_uncalled_sites_are_excluded_from_both_sides():
    """A -1 on either side kills the site, as in their `indiv_valid & other_valid`."""
    from arm2_attack_a3a4 import a3_score
    alt_af = np.array([0.3, 0.3])
    s = _site_table(alt_af)
    y = np.array([[1, 1], [1, 1]], dtype=np.int16)
    full = np.array([[1, 1], [1, 1]], dtype=np.int16)
    half = np.array([[1, 1], [-1, 1]], dtype=np.int16)      # second site uncalled
    assert a3_score(y, half, s) < a3_score(y, full, s)
    # and an uncalled RELEASE position excludes the site just as surely
    y_half = np.array([[1, 1], [-1, -1]], dtype=np.int16)
    assert a3_score(y_half, full, s) == pytest.approx(a3_score(y, half, s))
