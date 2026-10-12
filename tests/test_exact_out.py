"""Does asking for the exact endpoint change the release?

⚠⚠ WHY THIS IS THE GATE. Arm 1 needs both a DRAW and the exact `E[differ]`, and both
need the same `psi` and the same `forward` pass. Computing them in separate runs
recomputes both from scratch — measured 2026-10-12, that doubling was the difference
between ~17.4 h and ~10.9 h for arm 1's 308 cells. So `release()` now takes an
`exact_out` dict and fills it from the forward messages the draw already has.

That is only safe if it is INERT with respect to the draw. `state_marginals` consumes no
randomness and is placed before `backward_sample`, so the released path must be
bit-identical whether or not the exact quantities were requested. If it were not, every
arm 1 release would depend on whether we happened to be measuring it — and the two runs
would silently disagree, with nothing to show which was the real release.

⚠ The second thing tested is agreement: the exact values must match an INDEPENDENT
computation from the marginals, not merely be self-consistent.

⚠⚠ THE POWER OF THE FIRST GATE IS CHARACTERISED, NOT ASSUMED. Mutation-tested
2026-10-12. `backward_sample` builds a FRESH `default_rng(seed)`, so nothing before it
can perturb the draw except by corrupting `fwd` itself. Reversing `fwd.alpha` along the
state axis inside the marginals block is caught at 3 of 3 `tau`. Scaling one element of
`fwd.alpha[0]` by 1e-6 and renormalising is NOT caught, because inverse-CDF sampling
absorbs a perturbation that small. So this is a regression guard with power against
STRUCTURAL corruption of the forward messages and none against float-level noise. That
is the honest scope, and it is the scope that matters: a marginals pass that scribbles on
the messages would corrupt them structurally, not by 1e-6.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

from mechanism import calibration as C  # noqa: E402
from mechanism import ffbs as F  # noqa: E402
from mechanism import hmm as H  # noqa: E402
from mechanism import toy as TOY  # noqa: E402
from mechanism import utility as U  # noqa: E402

TAUS = (0.0, 0.25, 0.9)


@pytest.fixture(scope="module")
def toy():
    t = TOY.toy_cohort()
    util = U.beta_v(t.support)
    m = H.from_arrays(t.haplotypes, t.genetic_pos, t.positions, cohort_id="chrT")
    tgt = np.stack([t.haplotypes[0], t.haplotypes[1]], axis=1)
    return t, util, m, F.baseline_log_norm(m), tgt


@pytest.mark.parametrize("tau", TAUS)
def test_exact_out_does_not_perturb_the_draw(toy, tau):
    """THE GATE. Same seed, with and without `exact_out` -> bit-identical release."""
    t, util, m, base, tgt = toy
    kw = dict(seed=777, target_id="T", baseline_norm=base, stream_key="chrT")
    plain = F.release(m, t.haplotypes, tgt, util, tau, **kw)
    box: dict = {}
    withx = F.release(m, t.haplotypes, tgt, util, tau, exact_out=box, **kw)

    assert np.array_equal(plain.state_path, withx.state_path), (
        "the sampled state path changed when the exact endpoint was requested, so the "
        "marginals pass consumed randomness or reordered it")
    assert np.array_equal(plain.released_path, withx.released_path)
    assert plain.log_z_p == withx.log_z_p
    assert plain.u_achieved == withx.u_achieved
    assert box, "exact_out was left empty"


@pytest.mark.parametrize("tau", TAUS)
def test_exact_out_is_populated_for_both_haplotypes(toy, tau):
    t, util, m, base, tgt = toy
    box: dict = {}
    F.release(m, t.haplotypes, tgt, util, tau, seed=5, target_id="T",
              baseline_norm=base, exact_out=box)
    assert set(box) >= {"e_differ", "u_haps", "called"}
    assert len(box["u_haps"]) == 2, "a diploid release must contribute TWO haplotypes"
    assert box["called"] == int((tgt >= 0).sum())
    assert box["e_differ"] >= 0.0


@pytest.mark.parametrize("tau", TAUS)
def test_the_exact_values_match_an_independent_computation(toy, tau):
    """AGREEMENT, not self-consistency: recompute from the marginals separately."""
    t, util, m, base, tgt = toy
    box: dict = {}
    F.release(m, t.haplotypes, tgt, util, tau, seed=11, target_id="T",
              baseline_norm=base, exact_out=box)

    eta = C.eta_tau(tau)
    differ, u_haps = 0.0, []
    for col in (0, 1):
        t1 = tgt[:, col]
        psi = U.tilt_potentials(m.allowed, t.haplotypes, t1, util.beta, eta / 2.0,
                                unit_spans=m.unit_spans)
        gam = F.state_marginals(m, F.forward(m, psi), psi)
        differ += F.expected_mismatches(m, t.haplotypes, t1, gam)
        known = t1 >= 0
        match = (t.haplotypes == t1[None, :]) & known[None, :]
        u_haps.append(float((gam * match.T).sum(axis=1) @ util.beta))

    assert box["e_differ"] == pytest.approx(differ, rel=0, abs=1e-12)
    for a, b in zip(box["u_haps"], u_haps):
        assert a == pytest.approx(b, rel=0, abs=1e-12)


def test_at_tau_zero_the_endpoint_is_the_untilted_baseline(toy):
    """`tau` = 0 means eta = 0, so the potentials are all-ones and the marginals are
    the BASELINE posterior. The floor needs no special case, and this pins that."""
    t, util, m, base, tgt = toy
    box: dict = {}
    F.release(m, t.haplotypes, tgt, util, 0.0, seed=3, target_id="T",
              baseline_norm=base, exact_out=box)
    differ = 0.0
    for col in (0, 1):
        t1 = tgt[:, col]
        psi = U.tilt_potentials(m.allowed, t.haplotypes, t1, np.zeros_like(util.beta),
                                0.0, unit_spans=m.unit_spans)
        gam = F.state_marginals(m, F.forward(m, psi), psi)
        differ += F.expected_mismatches(m, t.haplotypes, t1, gam)
    assert box["e_differ"] == pytest.approx(differ, rel=0, abs=1e-12)


def test_the_tilt_reduces_expected_differ(toy):
    """DIRECTION. `target_fidelity` is a DIFFER rate, so a higher `tau` must not
    increase it: the tilt is toward the target, never away."""
    t, util, m, base, tgt = toy
    got = []
    for tau in (0.0, 0.5, 0.99):
        box: dict = {}
        F.release(m, t.haplotypes, tgt, util, tau, seed=9, target_id="T",
                  baseline_norm=base, exact_out=box)
        got.append(box["e_differ"])
    assert got[0] >= got[1] >= got[2], (
        f"E[differ] is not monotone non-increasing in tau: {got}. The tilt would be "
        f"pushing the release AWAY from the target.")
