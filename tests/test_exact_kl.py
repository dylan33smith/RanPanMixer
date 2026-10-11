"""Does the marginals route to KL agree with BRUTE-FORCE ENUMERATION?

⚠⚠ WHY THIS IS THE ONLY GATE THAT MATTERS FOR `scripts/exact_kl.py`. That script is
meant to CARRY the Arm 2 privacy claim, because sampling provably cannot: `sd(Lambda)`
per draw is 6.15e-04 against a true value of 1.784e-07, so resolving `E[Lambda]` by
sampling needs n ~ 1.2e7. A quantity nothing else can check must be checked against the
one thing that can -- the toy, where every supported output path can be enumerated and
the distributions written down in full.

⚠ The toy is K = 4, T = 7: 2,304 supported state paths. On chr21 the support is
astronomically large, which is exactly why the marginals route exists and why its
agreement HERE is what licenses its use THERE.
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(REPO / "scripts"))

from mechanism import calibration as C                    # noqa: E402
from mechanism import exact, ffbs, hmm, toy, utility       # noqa: E402

TAUS = (0.05, 0.25, 0.5, 0.9)


@pytest.fixture(scope="module")
def setup():
    t = toy.toy_cohort()
    p, q = toy.sharp_pair(t)
    base = exact.enumerate_baseline(t.haplotypes, t.genetic_pos)
    beta = utility.beta_v(t.haplotypes.shape[0] - (t.haplotypes == -1).sum(axis=0),
                          in_chain=np.ones(t.haplotypes.shape[1], bool)).beta
    return t, p, q, base, beta


def _kl_by_enumeration(base, beta, p, q, eta) -> float:
    """KL(Q_p || Q_q) by writing down both distributions in full."""
    u_p = exact.utility_vector(base, p, beta)
    u_q = exact.utility_vector(base, q, beta)
    qp, _ = exact.release_distribution(base, u_p, eta)
    qq, _ = exact.release_distribution(base, u_q, eta)
    m = qp > 0
    return float((qp[m] * (np.log(qp[m]) - np.log(qq[m]))).sum())


def _kl_by_marginals(t, beta, p, q, eta) -> float:
    """KL via the identity the script uses: eta(E u_p - E u_q) - (log Z_p - log Z_q),
    with BOTH expectations taken under Q_p's posterior marginals."""
    m = hmm.from_arrays(t.haplotypes, t.genetic_pos,
                        np.arange(t.haplotypes.shape[1]) * 1000, cohort_id="toy")
    b = ffbs.baseline_log_norm(m)
    util = utility.beta_v(t.haplotypes.shape[0] - (t.haplotypes == -1).sum(axis=0),
                          in_chain=np.ones(t.haplotypes.shape[1], bool))

    def terms(tgt):
        psi = utility.tilt_potentials(m.allowed, t.haplotypes, tgt, util.beta, eta)
        fwd = ffbs.forward(m, psi)
        lz = ffbs.log_z_p(fwd, b, eta, budget=eta * float(util.beta.sum()))
        return lz, ffbs.state_marginals(m, fwd, psi)

    def e_u(gamma, tgt):
        known = tgt >= 0
        match = (t.haplotypes == tgt[None, :]) & known[None, :]
        return float((gamma * match.T).sum(axis=1) @ util.beta)

    lz_p, gam_p = terms(p)
    lz_q, _ = terms(q)
    return eta * (e_u(gam_p, p) - e_u(gam_p, q)) - (lz_p - lz_q)


@pytest.mark.parametrize("tau", TAUS)
def test_marginals_kl_matches_brute_force_enumeration(setup, tau) -> None:
    """⚠ THE LICENSING TEST. If these disagree, every chr21 number from the marginals
    route is unsupported, because nothing at chr21 scale can check it."""
    t, p, q, base, beta = setup
    eta = C.eta_tau(tau)
    enum = _kl_by_enumeration(base, beta, p, q, eta)
    marg = _kl_by_marginals(t, beta, p, q, eta)
    assert enum >= -1e-15, f"enumerated KL is negative: {enum!r}"
    assert abs(enum - marg) < 1e-12 + 1e-9 * abs(enum), (
        f"tau={tau}: enumeration gives KL = {enum!r}, the marginals route gives "
        f"{marg!r}, difference {enum - marg!r}. The marginals route is what "
        f"scripts/exact_kl.py uses at chr21 scale where enumeration is impossible, so "
        f"a disagreement here invalidates every number it produces.")


@pytest.mark.parametrize("tau", TAUS)
def test_the_enumerated_kl_is_non_negative_and_pinsker_holds(setup, tau) -> None:
    """KL >= 0 is a theorem, and TV <= sqrt(KL/2) is the bound the script relies on to
    turn a KL into an attacker ceiling. ⚠ Checked against the ENUMERATED TV, which is
    the real thing rather than a bound."""
    t, p, q, base, beta = setup
    eta = C.eta_tau(tau)
    u_p = exact.utility_vector(base, p, beta)
    u_q = exact.utility_vector(base, q, beta)
    qp, _ = exact.release_distribution(base, u_p, eta)
    qq, _ = exact.release_distribution(base, u_q, eta)
    kl = _kl_by_enumeration(base, beta, p, q, eta)
    tv_true = exact.tv(qp, qq)
    assert kl >= 0.0
    assert tv_true <= math.sqrt(kl / 2.0) + 1e-12, (
        f"tau={tau}: Pinsker VIOLATED -- true TV {tv_true!r} exceeds "
        f"sqrt(KL/2) = {math.sqrt(kl/2.0)!r}. The script's attacker ceiling rests on "
        f"this inequality.")
    # and the ceiling the script would report must not understate the true one
    assert 0.5 * (1 + math.sqrt(kl / 2.0)) >= 0.5 * (1 + tv_true) - 1e-12


@pytest.mark.parametrize("tau", TAUS)
def test_the_true_tv_respects_tau(setup, tau) -> None:
    """⚠ Theorem 1 itself, on the toy: TV(Q_p, Q_q) <= tau for the SHARP pair, which is
    the pair constructed to make the bound as tight as it gets."""
    t, p, q, base, beta = setup
    eta = C.eta_tau(tau)
    qp, _ = exact.release_distribution(base, exact.utility_vector(base, p, beta), eta)
    qq, _ = exact.release_distribution(base, exact.utility_vector(base, q, beta), eta)
    assert exact.tv(qp, qq) <= tau + 1e-12


def test_tau_zero_gives_exactly_zero_kl(setup) -> None:
    """The floor. At tau = 0 both releases ARE the baseline, so the KL is exactly 0 and
    the attacker ceiling is exactly a coin."""
    t, p, q, base, beta = setup
    assert _kl_by_enumeration(base, beta, p, q, 0.0) == pytest.approx(0.0, abs=1e-15)
    assert C.eta_tau(0.0) == 0.0
