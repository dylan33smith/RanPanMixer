#!/usr/bin/env python
"""`var_q_u` -- the EXACT posterior variance of `u` under the release distribution.

⚠⚠ WHY THIS EXISTS. `logs/B_exact_utility.md` records that `e_q_u` moves +8.345e-08
genome-wide against +3.6368e-06 on chr21 at share 1 -- a factor of **43.6** -- and that
TWO first-order predictions of that factor both failed, in opposite directions
(`mean beta_v` predicts 72, `sum(beta_v^2)` predicts 128.6). `docs/plan.md` states the
untested explanation: allele matches are POSITIVELY CORRELATED along the copying path,
so `Var_Q(u)` exceeds `sum(beta_v^2)`. This script measures `Var_Q(u)` and settles it.

THE IDENTITY THIS RESTS ON. `Q_eta(y) ∝ R_D(y) exp(eta * u(p,y))` is an exponential
family in `eta` with sufficient statistic `u`, so

    d E_Qeta[u] / d eta = Var_Qeta(u).

The movement across the grid is therefore `integral_0^eta_tau Var_Qeta(u) d eta`, and
because `Var` barely moves over that range the first-order form `eta_tau * Var_Q0(u)`
is already tight. So this script does not merely EXPLAIN the measured movement -- it
PREDICTS it from an independently computed quantity, and the comparison is a gate.

⚠ THE ROUTE IS NOT PAIR MARGINALS. `docs/plan.md` and `scripts/exact_utility.py` both
say settling this "needs `Var_Q(u)` from PAIR marginals". It does not. Pair marginals
over all position pairs are O(T^2 K^2) and unreachable at T = 22,134,562. `u` is an
ADDITIVE functional of the hidden path, and FFBS's backward chain is a Markov chain
whose kernels are already in hand, so a single O(TK) sweep carrying the first AND
second conditional moments of the running sum returns the variance exactly. That is
this file. ⚠ The pair-marginal wording in those two files is therefore wrong about the
method, not about the quantity -- correct it there rather than here.

⚠⚠ THE RECURSION IS CENTRED, AND THAT IS LOAD-BEARING. `Var = E[u^2] - E[u]^2` is a
cancellation of about eight orders of magnitude here (`E[u]` ~ 0.59, `Var` ~ 1e-06),
which throws away half of float64's digits for no reason. Subtracting the per-position
posterior mean FIRST makes the recursion compute `E[(u - E[u])^2]` directly, with no
cancellation at all, because `sum_t E[f_t(Z_t)] = E[u]` exactly by linearity. The
centred first moment must then come back at 0, which is a free instrument check on a
quantity whose correctness is otherwise hard to see -- and it is computed through the
alpha-kernels while the centring constants come from the gamma-marginals, so the two
routes share no arithmetic.

WHAT IT REPORTS, and `docs/terms.md` carries all four:
  `var_q_u`            the exact posterior variance of `u` (diploid: of `u_dip`)
  `var_q_u_diag`       what `Var` would be if positions were INDEPENDENT
  `cov_share`          `1 - var_q_u_diag / var_q_u`: the covariance's share of it
  `var_q_differ`       the exact posterior variance of the differing-cell COUNT,
                       whose square root is the SD of every "1 draw" table column

    python scripts/exact_variance.py --validate            # toy + mutation test
    python scripts/exact_variance.py --chroms chr21 --share1 --mc 30
    python scripts/exact_variance.py --taus 0 0.99         # genome-wide, 22 autosomes
"""

from __future__ import annotations

import argparse
import json
import math
import sys
import time
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

import numpy as np                                        # noqa: E402

from mechanism import calibration as C                    # noqa: E402
from mechanism import exact, ffbs, hmm, toy, utility      # noqa: E402

AUTOSOMES = [f"chr{c}" for c in range(1, 23)]


# ---------------------------------------------------------------- the recursion

def centred_moments(switch_mass: np.ndarray, K: int, alpha: np.ndarray,
                    f_of_t, centre: np.ndarray) -> tuple[float, float]:
    """`(E[S], E[S^2])` for `S = sum_t (f_t(Z_t) - centre_t)` under the posterior.

    `Z` is distributed as FFBS's backward chain: `Z_{T-1} ~ alpha_{T-1}` and then
    `P(Z_{t-1} = j | Z_t = i) ∝ alpha_{t-1}(j) * A_t(j, i)`. Conditional on
    `Z_t`, the states BEFORE `t` are independent of those after it -- that is the
    Markov property of the backward chain, in the direction it is sampled -- so
    the prefix moments can be carried forward one position at a time:

        P_t(i) = E[S_t | Z_t = i]   = g_t(i) + M_t(i)
        R_t(i) = E[S_t^2 | Z_t = i] = g_t(i)^2 + 2 g_t(i) M_t(i) + N_t(i)

    with `M_t` and `N_t` the kernel averages of `P_{t-1}` and `R_{t-1}`. The
    stay-switch kernel is rank-one plus diagonal, so each average is one scalar
    dot product plus one elementwise term -- O(K) per position, never a matrix:

        D_t(i) = s/K + (1 - s) alpha_{t-1}(i)                   (alpha sums to 1)
        M_t(i) = [ (s/K) <alpha_{t-1}, P_{t-1}> + (1-s) alpha_{t-1}(i) P_{t-1}(i) ] / D_t(i)

    ⚠ `A` IS NEVER MATERIALISED, for the reason `ffbs.py` point 3 gives: a dense
    per-position matrix is about 21.1 GB on chr21 alone.

    ⚠⚠ `D_t(i)` REACHES EXACTLY ZERO ON REAL DATA, and the unguarded division
    returns NaN for the whole chromosome. It needs two things at once: a step with
    `switch_mass == 0` (zero genetic distance, where `unif` vanishes and the
    transition is the identity) AND a donor with no call at the SOURCE position, so
    `alpha_{t-1}(i) == 0` as well. chr21 has 1,058 zero-distance steps and 7 of them
    sit next to an unavailable donor; the first version of this file divided anyway
    and reported `var_q_u` = nan for all 332,521 positions. ⚠ The project's toy
    cohort carries BOTH features -- a zero-distance step at index 2 and two uncalled
    cells -- but never together, so the enumeration gate passed while the real run
    was NaN. `validate()` therefore builds a second, smaller cohort where they
    coincide, and that case is the gate for this line.

    Such a state is genuinely unreachable: with the identity transition `Z_t = i`
    requires `Z_{t-1} = i`, which has zero forward mass. Its conditional moments
    are multiplied by `alpha_t(i) == 0` everywhere they are used, so any finite
    value does, and the assertion below is what says so rather than assuming it.

    ⚠ `f_of_t(t)` returns the (K,) vector of per-state values at position `t`. It
    is a callable rather than a (T, K) array because that array is 1.4 GB per
    chromosome at genome scale, on a box that is already holding `alpha` and the
    backward messages.
    """
    T = alpha.shape[0]
    g = f_of_t(0) - centre[0]
    P = g.copy()
    R = g * g
    for t in range(1, T):
        sw = switch_mass[t - 1]
        unif = sw / K
        stay = 1.0 - sw
        a_prev = alpha[t - 1]
        D = unif + stay * a_prev
        dead = D <= 0.0
        if dead.any():
            live = alpha[t][dead]
            if (live > 0.0).any():
                raise AssertionError(
                    f"chain position {t}: the backward kernel's normaliser is 0 for "
                    f"{int(dead.sum())} state(s), but {int((live > 0).sum())} of them "
                    f"carry forward mass at {t} (max {float(live.max()):.3e}). With "
                    f"switch_mass = {sw!r} the transition is the identity, so such a "
                    f"state is reachable only from itself at {t-1}, where it has "
                    f"none. The forward pass and this recursion disagree about the "
                    f"chain -- do NOT paper over it here.")
            D = np.where(dead, 1.0, D)
        aP = a_prev * P
        aR = a_prev * R
        M = (unif * aP.sum() + stay * aP) / D
        N = (unif * aR.sum() + stay * aR) / D
        if dead.any():
            # Finite, so the dot products stay finite; weighted by alpha_t == 0.
            M = np.where(dead, 0.0, M)
            N = np.where(dead, 0.0, N)
        g = f_of_t(t) - centre[t]
        P = g + M
        R = g * g + 2.0 * g * M + N
    aT = alpha[T - 1]
    return float(aT @ P), float(aT @ R)


def haplotype_variance(m, alpha: np.ndarray, gamma: np.ndarray,
                       indicator: np.ndarray, weights: np.ndarray) -> dict:
    """`Var` and `E` of `sum_v weights_v * indicator_v(donor)` for ONE haplotype.

    `indicator` is (K, T) boolean -- match-with-the-target for `u`, differs-from-it
    for the cell count -- and `weights` is (T,). Both quantities this script reports
    are of that shape, so they share one implementation and one gate.

    ⚠ MODEL C ONLY. Under model B a unit spans several positions and `f_t` is a
    SUM over the unit, which this does not do; `m.unit_spans is None` is asserted
    rather than assumed, because the silent version would return a variance for a
    different quantity than the one named.
    """
    if m.unit_spans is not None:
        raise NotImplementedError(
            "model B is not supported: a unit spans several positions, so f_t is a "
            "sum over the unit and the per-position centring below is not the right "
            "quantity. The live build is model C (docs/plan.md decision C).")
    T, K = m.T, m.K
    if indicator.shape != (K, T):
        raise ValueError(f"indicator must be (K, T) = {(K, T)}, got {indicator.shape}")
    if weights.shape != (T,):
        raise ValueError(f"weights must be (T,) = {(T,)}, got {weights.shape}")

    # P(indicator at t) under the posterior, then the per-position mean. Exact,
    # and the same expression `scripts/exact_utility.py` uses for E_Q[u].
    p_t = (gamma * indicator.T).sum(axis=1)
    centre = weights * p_t
    mean = float(math.fsum(centre))

    ind = indicator            # local, so the closure does not look up a global
    w = weights

    def f_of_t(t: int) -> np.ndarray:
        return w[t] * ind[:, t]

    centred_mean, var = centred_moments(m.switch_mass, K, alpha, f_of_t, centre)

    # If positions were independent the variance would be the diagonal alone.
    # This is the `sum(beta_v^2)`-style prediction, computed on the SAME posterior
    # rather than assumed, so the covariance share below is a measurement.
    var_diag = float(math.fsum(w * w * p_t * (1.0 - p_t)))

    return {"mean": mean, "var": var, "var_diag": var_diag,
            "centred_mean": centred_mean}


# ------------------------------------------------------------------ one chromosome

def chrom_variance_loo(chrom: str, target: str, tau: float, share: float) -> dict:
    """The same quantities for ARM 1: `target` is a cohort member, held out at K = 86.

    ⚠ ARM 1'S TARGET PATH IS ITS OWN TWO HAPLOTYPE ROWS -- it is a cohort member, so
    no target VCF exists. The panel, the chain, `support` and therefore `beta_v` are
    all the LEAVE-ONE-OUT ones, rebuilt per target, because decision K is re-applied
    per hold-out and the exclusion set is NOT monotone in the donor panel
    (`logs/B_loo_chain_survey.md`). Subsetting arm 2's chain here is the obvious
    implementation and it is wrong.

    ⚠ `share` MUST BE THIS TARGET'S OWN, not arm 2's: leave-one-out moves every
    `T_c`, so each of the 44 targets has its own 22 shares (`loo.loo_shares`). The
    caller passes it, and `logs/B_loo_chain_survey.json` has them per target.
    """
    from cohort import loo as loo_mod

    d = REPO / "data" / "cohort" / chrom
    lc = loo_mod.loo_chain(d, target)
    ic = lc.in_chain
    hap_all = np.asarray(np.load(d / "haplotypes.npy", mmap_mode="r"))
    panel = loo_mod.masked_haplotypes(d, target, haplotypes=hap_all)[:, ic]
    rows = loo_mod.sample_index(
        np.load(d / "haplotype_ids.npy", allow_pickle=True))[target]
    p = hap_all[list(rows)][:, ic].T
    m = hmm.from_arrays(panel, np.load(d / "genetic_pos.npy")[ic],
                        np.load(d / "positions.npy")[ic],
                        cohort_id=f"{chrom}_loo_{target}")
    util = utility.beta_v(lc.support, in_chain=ic, share=share, chrom=chrom)
    return _variance_body(m, panel, p, util, chrom, target, tau, share,
                          arm="arm1", k_states=lc.k_states)


def chrom_variance(chrom: str, target: str, tau: float, share: float,
                   *, loo: str | None = None, keep_alpha: bool = False) -> dict:
    """Every reported quantity for one (chromosome, target, tau). Exact.

    Diploid, like the release itself: two INDEPENDENT chains at `eta_tau / 2`
    (`ffbs.release`), `u_dip = (u1 + u2) / 2`, so

        Var(u_dip) = (Var(u1) + Var(u2)) / 4,
        Var(differing cells) = Var(h1) + Var(h2)      -- a SUM, not a mean.

    ⚠ The halves are not interchangeable. Averaging the cell counts would report a
    variance 4x too small for the row it is meant to be the SD of, and summing the
    utilities would report one 4x too large for a quantity the bound covers.
    """
    d = REPO / "data" / "cohort" / chrom
    ic = np.load(d / "in_chain.npy").astype(bool)
    hap = np.load(d / "haplotypes.npy")[:, ic]
    support = np.load(d / "support.npy")
    gpos = np.load(d / "genetic_pos.npy")[ic]
    pos = np.load(d / "positions.npy")[ic]
    p = np.load(REPO / "data" / "targets" / target / chrom / "path.npy")[ic]

    if loo is not None:
        raise ValueError("call chrom_variance_loo for arm 1, not this with --loo")

    m = hmm.from_arrays(hap, gpos, pos, cohort_id=chrom)
    util = utility.beta_v(support, in_chain=ic, share=share, chrom=chrom)
    return _variance_body(m, hap, p, util, chrom, target, tau, share,
                          arm="arm2", k_states=int(m.K))


def _variance_body(m, hap: np.ndarray, p: np.ndarray, util, chrom: str,
                   target: str, tau: float, share: float, *, arm: str,
                   k_states: int) -> dict:
    """The arm-independent half: two haplotypes, both quantities, one `tau`.

    Split out so arm 1 and arm 2 cannot drift apart in the arithmetic while
    differing, as they must, in the panel, the chain and the shares above.
    """
    eta = C.eta_tau(tau)
    out = {"chrom": chrom, "target": target, "tau": tau, "share": share,
           "arm": arm, "T": int(m.T), "K": int(k_states), "eta_tau": eta,
           "beta_sum": float(util.beta.sum()),
           "beta_sq_sum": float(math.fsum(util.beta * util.beta))}

    v_u = v_u_diag = v_dif = e_u = e_dif = 0.0
    centred = []
    for col in (0, 1):
        tgt = p[:, col]
        known = tgt >= 0
        psi = utility.tilt_potentials(m.allowed, hap, tgt, util.beta, eta / 2.0)
        fwd = ffbs.forward(m, psi)
        gam = ffbs.state_marginals(m, fwd, psi)

        match = (hap == tgt[None, :]) & known[None, :]
        mism = (hap != tgt[None, :]) & known[None, :]

        ru = haplotype_variance(m, fwd.alpha, gam, match, util.beta)
        rd = haplotype_variance(m, fwd.alpha, gam, mism,
                                np.ones(m.T, dtype=np.float64))

        # (u1 + u2)/2 -> variances add at 1/4 each; cell counts add outright.
        v_u += 0.25 * ru["var"]
        v_u_diag += 0.25 * ru["var_diag"]
        e_u += 0.5 * ru["mean"]
        v_dif += rd["var"]
        e_dif += rd["mean"]
        centred += [ru["centred_mean"], rd["centred_mean"]]
        out[f"called_h{col}"] = int(known.sum())
        del psi, fwd, gam, match, mism

    out.update({"e_q_u": e_u, "var_q_u": v_u, "var_q_u_diag": v_u_diag,
                "e_differ": e_dif, "var_q_differ": v_dif,
                "cov_share": 1.0 - v_u_diag / v_u if v_u > 0 else float("nan"),
                "max_abs_centred_mean": max(abs(x) for x in centred)})
    return out


# ----------------------------------------------------------------- the gates

def validate_dead_kernel() -> float:
    """The gate for the zero-normaliser line, against enumeration. Returns the worst
    relative difference.

    ⚠ THIS COHORT EXISTS BECAUSE THE PROJECT'S TOY CANNOT CATCH THE BUG. The toy
    has a zero-distance step and it has uncalled cells, but not at the same place,
    so `D_t(i) = 0` never occurs there and the first version of this file passed the
    enumeration gate while returning NaN for every real chromosome. Here donor 3 is
    uncalled at position 1 and the step from 1 to 2 carries zero genetic distance,
    which is the coincidence chr21 hits 7 times.

    ⚠ MEASURED, 2026-10-22: with the unguarded division this construction returns
    `(nan, nan)` while enumeration returns mean 0.571546 and variance 0.056993. So
    the gate fires on the defect it was written for, rather than merely passing.
    """
    hap = np.array([[0, 0, 1, 0, 1],
                    [1, 1, 0, 1, 0],
                    [0, 1, 1, 1, 0],
                    [1, -1, 0, 0, 1]], dtype=np.int16)
    gpos = np.array([0.0, 1e-5, 1e-5, 3e-5, 6e-5])
    m = hmm.from_arrays(hap, gpos, np.arange(5, dtype=np.int64), cohort_id="dead")
    if m.switch_mass[1] != 0.0 or m.allowed[3, 1]:
        raise AssertionError(
            f"this cohort no longer carries the coincidence it exists for: "
            f"switch_mass[1] = {m.switch_mass[1]!r} (want 0.0) and "
            f"allowed[3,1] = {m.allowed[3,1]!r} (want False). Fix the construction, "
            f"not the assertion -- without both, the gate tests nothing.")
    base = exact.enumerate_baseline(hap, gpos)
    beta = np.full(5, 1.0 / 5)
    tgt = np.array([0, 1, 1, 0, 1], dtype=np.int16)
    worst = 0.0
    print(f"\ndead-kernel cohort: K = {m.K}, T = {m.T}, support M = {base.M}, "
          f"switch_mass[1] = 0 with donor 3 uncalled at position 1")
    print("\n| metric | " + " | ".join(f"tau={x:g}" for x in (0.0, 0.5, 0.99)) + " |")
    print("|" + "---|" * 4)
    cells: dict[str, list[str]] = {"var, enumerated": [], "var, recursion": [],
                                   "rel. difference": [], "all finite?": []}
    for tau in (0.0, 0.5, 0.99):
        eta = C.eta_tau(tau)
        u_vec = exact.utility_vector(base, tgt, beta)
        q, _ = exact.release_distribution(base, u_vec, eta)
        e_enum = float(q @ u_vec)
        v_enum = float(q @ (u_vec * u_vec) - e_enum ** 2)
        psi = utility.tilt_potentials(m.allowed, hap, tgt, beta, eta)
        fwd = ffbs.forward(m, psi)
        gam = ffbs.state_marginals(m, fwd, psi)
        match = (hap == tgt[None, :]) & (tgt[None, :] >= 0)
        r = haplotype_variance(m, fwd.alpha, gam, match, beta)
        rel = abs(r["var"] - v_enum) / max(v_enum, 1e-300)
        worst = max(worst, rel)
        cells["var, enumerated"].append(f"{v_enum:.12e}")
        cells["var, recursion"].append(f"{r['var']:.12e}")
        cells["rel. difference"].append(f"{rel:.2e}")
        cells["all finite?"].append(
            "yes" if math.isfinite(r["var"]) and math.isfinite(r["centred_mean"])
            else "**NO**")
    for k, v in cells.items():
        print(f"| {k} | " + " | ".join(v) + " |")
    if worst > 1e-10:
        raise AssertionError(
            f"the recursion disagrees with enumeration by {worst:.3e} on the cohort "
            f"where the backward kernel's normaliser vanishes. That is the real "
            f"chain's case, not a corner: chr21 hits it 7 times.")
    return worst


def validate() -> int:
    """The toy gate and its mutation test. Nothing here reads the real cohort.

    ⚠ A GATE THAT CANNOT FAIL IS DECORATION (`CLAUDE.md`), and for a variance the
    temptation is exact: the independence formula `sum_v beta_v^2 p_v (1 - p_v)`
    agrees with the truth to several digits whenever the chain switches often, so a
    gate that only checked "is it about the right size" would pass on the wrong
    quantity -- which is the very error this script exists to settle. So the gate
    compares against BRUTE-FORCE ENUMERATION of the release distribution, and then
    asserts that the independence formula FAILS that same comparison.
    """
    t = toy.toy_cohort()
    base = exact.enumerate_baseline(t.haplotypes, t.genetic_pos)
    m = hmm.from_arrays(t.haplotypes, t.genetic_pos,
                        np.arange(t.T, dtype=np.int64), cohort_id="toy")
    beta = np.full(t.T, 1.0 / t.T)
    inputs = toy.toy_inputs(t)
    print(f"toy: K = {t.K}, T = {t.T}, support M = {base.M}, "
          f"{len(inputs)} input path(s)\n")

    print("| metric | " + " | ".join(f"input {i}, tau={x:g}"
                                     for i in range(min(2, len(inputs)))
                                     for x in (0.0, 0.5, 0.99)) + " |")
    rows: dict[str, list[str]] = {"var, enumerated": [], "var, recursion": [],
                                  "rel. difference": [], "var, independence": [],
                                  "centred mean (must be 0)": []}
    worst = 0.0
    worst_indep = 0.0
    for i in range(min(2, len(inputs))):
        tgt = np.asarray(inputs[i]).astype(np.int16)
        for tau in (0.0, 0.5, 0.99):
            eta = C.eta_tau(tau)
            u_vec = exact.utility_vector(base, tgt, beta)
            q, _ = exact.release_distribution(base, u_vec, eta)
            e_enum = float(q @ u_vec)
            v_enum = float(q @ (u_vec * u_vec) - e_enum ** 2)

            psi = utility.tilt_potentials(m.allowed, t.haplotypes, tgt, beta, eta)
            fwd = ffbs.forward(m, psi)
            gam = ffbs.state_marginals(m, fwd, psi)
            match = (t.haplotypes == tgt[None, :]) & (tgt[None, :] >= 0)
            r = haplotype_variance(m, fwd.alpha, gam, match, beta)

            rel = abs(r["var"] - v_enum) / max(v_enum, 1e-300)
            rel_indep = abs(r["var_diag"] - v_enum) / max(v_enum, 1e-300)
            worst = max(worst, rel)
            worst_indep = max(worst_indep, rel_indep)
            if abs(r["mean"] - e_enum) > 1e-12 * max(1.0, abs(e_enum)):
                raise AssertionError(
                    f"E[u] disagrees with enumeration at input {i}, tau = {tau}: "
                    f"{r['mean']!r} against {e_enum!r}. The variance cannot be "
                    f"right if the mean is not.")
            rows["var, enumerated"].append(f"{v_enum:.12e}")
            rows["var, recursion"].append(f"{r['var']:.12e}")
            rows["rel. difference"].append(f"{rel:.2e}")
            rows["var, independence"].append(f"{r['var_diag']:.6e}")
            rows["centred mean (must be 0)"].append(f"{r['centred_mean']:+.1e}")

    print("|" + "---|" * (len(rows["var, recursion"]) + 1))
    for k, v in rows.items():
        print(f"| {k} | " + " | ".join(v) + " |")

    print(f"\nworst relative difference against enumeration: {worst:.3e}")
    if worst > 1e-10:
        raise AssertionError(
            f"the recursion disagrees with brute-force enumeration by {worst:.3e}. "
            f"Both compute the same posterior variance; one of them is wrong.")

    print(f"⚠ MUTATION TEST -- the independence formula, which is the hypothesis "
          f"this script tests, is off by {worst_indep:.3e} on the same toy.")
    if worst_indep < 1e-6:
        raise AssertionError(
            f"the independence formula agrees with the truth to {worst_indep:.3e} on "
            f"this toy, so the gate above cannot tell the two apart and proves "
            f"nothing about the real chain. Pick a toy whose chain switches rarely "
            f"enough to carry covariance.")
    print("✅ the gate separates the two quantities, so it has power.")

    worst_dead = validate_dead_kernel()
    print(f"\n✅ zero-normaliser case agrees with enumeration to {worst_dead:.3e}; "
          f"the unguarded division returned nan here (measured 2026-10-22).")
    return 0


def monte_carlo(chrom: str, target: str, tau: float, share: float, n: int,
                seed: int = 20261022) -> int:
    """Exact `var_q_differ` against the empirical variance of `n` real draws.

    ⚠ ONE forward pass, `n` backward draws from it. That is the right design for
    this check and not only the cheap one: the draws must come from the SAME
    posterior the exact number describes, and recomputing `psi` per draw would
    leave a disagreement ambiguous between the recursion and the setup.

    ⚠ This is a HAPLOID check, deliberately. The diploid combination is two lines
    of arithmetic on top (see `chrom_variance`), and a gate that bundled them could
    not say which half was wrong.
    """
    d = REPO / "data" / "cohort" / chrom
    ic = np.load(d / "in_chain.npy").astype(bool)
    hap = np.load(d / "haplotypes.npy")[:, ic]
    m = hmm.from_arrays(hap, np.load(d / "genetic_pos.npy")[ic],
                        np.load(d / "positions.npy")[ic], cohort_id=chrom)
    util = utility.beta_v(np.load(d / "support.npy"), in_chain=ic,
                          share=share, chrom=chrom)
    p = np.load(REPO / "data" / "targets" / target / chrom / "path.npy")[ic]
    tgt = p[:, 0]
    known = tgt >= 0
    eta = C.eta_tau(tau)

    t0 = time.time()
    psi = utility.tilt_potentials(m.allowed, hap, tgt, util.beta, eta / 2.0)
    fwd = ffbs.forward(m, psi)
    gam = ffbs.state_marginals(m, fwd, psi)
    mism = (hap != tgt[None, :]) & known[None, :]
    match = (hap == tgt[None, :]) & known[None, :]
    rd = haplotype_variance(m, fwd.alpha, gam, mism, np.ones(m.T))
    ru = haplotype_variance(m, fwd.alpha, gam, match, util.beta)
    print(f"{chrom} {target} tau={tau:g} share={share:g} haplotype 0, "
          f"T = {m.T:,}, K = {m.K} ({time.time() - t0:.0f}s)\n")

    rng = np.random.default_rng(seed)
    dif = np.empty(n)
    us = np.empty(n)
    for i in range(n):
        z = ffbs.backward_sample(m, fwd, rng)
        rel = ffbs.expand_to_positions(m, hap, z)
        dif[i] = int(((rel != tgt) & known).sum())
        us[i] = float(utility.u_of_output(rel, tgt, util.beta)[0])
        if (i + 1) % 10 == 0:
            print(f"    {i + 1}/{n} draws ({time.time() - t0:.0f}s)", flush=True)

    print(f"\n{chrom}, {target}, tau = {tau:g}, haplotype 0, n = {n} draws from ONE "
          f"posterior.\n")
    print("| metric | exact | from draws | ratio |")
    print("|---|---|---|---|")
    for name, exact_m, exact_v, obs in (
            ("differing cells", rd["mean"], rd["var"], dif),
            ("`u` (haploid)", ru["mean"], ru["var"], us)):
        se_mean = math.sqrt(exact_v / n)
        print(f"| E[{name}] | {exact_m:.6f} | {obs.mean():.6f} "
              f"| {(obs.mean() - exact_m) / se_mean:+.2f} SE |")
        print(f"| SD[{name}] | {math.sqrt(exact_v):.6f} | {obs.std(ddof=1):.6f} "
              f"| {obs.std(ddof=1) / math.sqrt(exact_v):.4f} |")
    # The sampling distribution of a variance ratio at n draws: chi-square with
    # n-1 degrees of freedom has relative SD sqrt(2/(n-1)) on the variance, so
    # sqrt(1/(2(n-1))) on the SD. Three of those is a wide but finite gate.
    tol = 3.0 * math.sqrt(1.0 / (2.0 * (n - 1)))
    for name, exact_v, obs in (("differing cells", rd["var"], dif),
                               ("u", ru["var"], us)):
        ratio = obs.std(ddof=1) / math.sqrt(exact_v)
        print(f"\n{name}: SD ratio {ratio:.4f}, gate |ratio - 1| <= {tol:.4f} "
              f"(3 sigma at n = {n})")
        if abs(ratio - 1.0) > tol:
            raise AssertionError(
                f"the exact SD for {name} disagrees with {n} real draws: ratio "
                f"{ratio:.4f} against a 3-sigma gate of {tol:.4f}. Either the "
                f"recursion is wrong or the draws do not come from the posterior "
                f"it describes.")
    print("\n✅ the exact variance agrees with real draws on both quantities.")
    return 0


# -------------------------------------------------------------------- reporting

def report(rows: list[dict], taus: list[float], scope: str,
           measured: dict | None) -> None:
    """Metrics as ROWS, configurations as COLUMNS (`CLAUDE.md` reporting format)."""
    n = len({r["target"] for r in rows})
    tot_T = sum(r["T"] for r in rows if r["tau"] == taus[0])
    print(f"\n{scope}: exact posterior variance of `u`. n = {n} target(s), "
          f"T = {tot_T:,}, ZERO Monte-Carlo error.\n")
    print("| metric | " + " | ".join(f"tau = {x:g}" for x in taus) + " |")
    print("|" + "---|" * (len(taus) + 1))

    def col(key, fmt=".9e"):
        out = []
        for x in taus:
            vals = [r[key] for r in rows if r["tau"] == x]
            assert vals, f"no rows at tau = {x}; the report would be silent"
            out.append(format(float(np.mean(vals)), fmt))
        return out

    # `var` and `e` are SUMS over chromosomes: one release, one budget. The
    # chromosomes are independent given the panel, so their variances add.
    agg: dict[float, dict] = {}
    for x in taus:
        rs = [r for r in rows if r["tau"] == x]
        agg[x] = {
            "var_q_u": math.fsum(r["var_q_u"] for r in rs),
            "var_q_u_diag": math.fsum(r["var_q_u_diag"] for r in rs),
            "var_q_differ": math.fsum(r["var_q_differ"] for r in rs),
            "e_q_u": math.fsum(r["e_q_u"] for r in rs),
            "e_differ": math.fsum(r["e_differ"] for r in rs),
            "beta_sq_sum": math.fsum(r["beta_sq_sum"] for r in rs),
            "max_centred": max(r["max_abs_centred_mean"] for r in rs),
            "eta_tau": rs[0]["eta_tau"],
        }

    def arow(label, key, fmt=".9e"):
        print(f"| {label} | " + " | ".join(format(agg[x][key], fmt)
                                           for x in taus) + " |")

    arow("`var_q_u` *", "var_q_u")
    arow("`var_q_u_diag` (independence)", "var_q_u_diag")
    print("| `cov_share` | " + " | ".join(
        f"{1.0 - agg[x]['var_q_u_diag'] / agg[x]['var_q_u']:.6f}"
        if agg[x]["var_q_u"] > 0 else "n/a -- Var = 0 at tau = 0 only if u is "
        "degenerate" for x in taus) + " |")
    arow("`sum beta_v^2`", "beta_sq_sum")
    arow("`var_q_differ` *", "var_q_differ")
    print("| `sd_q_differ` (cells) | " + " | ".join(
        f"{math.sqrt(agg[x]['var_q_differ']):,.1f}" for x in taus) + " |")
    arow("`e_q_u`", "e_q_u", ".9f")
    arow("`e_differ`", "e_differ", ".1f")
    print("| `eta_tau`", end="")
    print(" | " + " | ".join(f"{agg[x]['eta_tau']:.6f}" for x in taus) + " |")
    print("| max abs centred mean (gate: ~0) | " + " | ".join(
        f"{agg[x]['max_centred']:.1e}" for x in taus) + " |")

    lo, hi = min(taus), max(taus)
    if lo != hi:
        print(f"\n**THE PREDICTION.** `d E_Q[u] / d eta = Var_Q(u)`, so the movement "
              f"across the grid is the integral of the variance.")
        v_lo, v_hi = agg[lo]["var_q_u"], agg[hi]["var_q_u"]
        de = agg[hi]["eta_tau"] - agg[lo]["eta_tau"]
        first = de * v_lo
        trap = de * 0.5 * (v_lo + v_hi)
        print(f"\n| quantity | value |")
        print(f"|---|---|")
        print(f"| `eta_tau` span | {de:.6f} |")
        print(f"| first order, `d_eta * Var(tau = {lo:g})` | {first:.6e} |")
        print(f"| trapezoid over the two endpoints | {trap:.6e} |")
        print(f"| measured `e_q_u` movement, THIS run | "
              f"{agg[hi]['e_q_u'] - agg[lo]['e_q_u']:.6e} |")
        if measured:
            for k, v in measured.items():
                print(f"| {k} | {v:.6e} |")
                print(f"| trapezoid / that | {trap / v:.4f} |")
        print(f"\n⚠ The two endpoints bracket the integral only if `Var` is monotone "
              f"in `eta` over the range, which is not proven here. The gap between "
              f"the first-order and trapezoid figures is "
              f"{abs(trap - first) / max(trap, 1e-300):.2%} of the trapezoid, so at "
              f"this scope the distinction does not carry the conclusion.")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--validate", action="store_true",
                    help="toy enumeration gate + mutation test; reads no real data")
    ap.add_argument("--mc", type=int, default=0,
                    help="also check the exact SD against this many real draws")
    ap.add_argument("--targets", type=int, default=1)
    ap.add_argument("--only-targets", nargs="+", default=None)
    ap.add_argument("--taus", type=float, nargs="+", default=[0.0, 0.99])
    ap.add_argument("--chroms", nargs="+", default=AUTOSOMES)
    ap.add_argument("--share1", action="store_true",
                    help="give the named chromosome the WHOLE budget (share = 1), "
                         "which is the chr21 scope logs/B_beta_concentration.md used")
    ap.add_argument("--loo", metavar="MEMBER", default=None,
                    help="ARM 1: hold out this cohort member at K = 86 and use its own "
                         "haplotypes as the target. Needs --share, because a held-out "
                         "target has its OWN 22 shares (logs/B_loo_chain_survey.json).")
    ap.add_argument("--share", type=float, default=None,
                    help="explicit beta_share for the single named chromosome, e.g. "
                         "0.015024898412119904 for HG00438's chr21 under arm 1")
    ap.add_argument("--json", type=Path, default=None)
    a = ap.parse_args(argv)

    if a.validate:
        return validate()

    T = {c: int(np.load(REPO / "data" / "cohort" / c / "in_chain.npy").sum())
         for c in AUTOSOMES}
    tot = sum(T.values())
    if a.share is not None:
        if len(a.chroms) != 1:
            ap.error("--share names ONE chromosome's share; pass one --chroms")
        if a.share1:
            ap.error("--share and --share1 both set the share; pass one")
        sh = {a.chroms[0]: a.share}
        scope = f"{a.chroms[0]} at share {a.share:.9g}"
    elif a.share1:
        if len(a.chroms) != 1:
            ap.error("--share1 is a single-chromosome scope; pass one --chroms")
        sh = {a.chroms[0]: 1.0}
        scope = f"{a.chroms[0]} at SHARE 1 (the whole budget on one chromosome)"
    else:
        sh = {c: T[c] / tot for c in a.chroms}
        scope = ("GENOME-WIDE, 22 autosomes" if len(a.chroms) == 22
                 else f"{len(a.chroms)} chromosome(s) at genome-wide share")

    if a.loo:
        if a.share is None:
            ap.error(
                "--loo needs --share: leave-one-out moves every T_c, so a held-out "
                "target has its OWN 22 shares and arm 2's would make sum_v beta_v != 1 "
                "over the scope tau covers (standing constraint 2). "
                "logs/B_loo_chain_survey.json carries them per target.")
        if a.only_targets or a.mc:
            ap.error("--loo IS the target and --mc is arm-2 only; pass neither")
        tids = [a.loo]
        scope = f"ARM 1, {a.loo} held out at K = 86, " + scope
    elif a.only_targets:
        tids = list(a.only_targets)
    else:
        tids = sorted(p.name for p in (REPO / "data" / "targets").iterdir()
                      if p.is_dir())[: a.targets]

    print(f"{scope}: share covered = {math.fsum(sh.values()):.6f}, "
          f"T = {sum(T[c] for c in a.chroms):,} of {tot:,}")
    print(f"targets: {', '.join(tids)}")
    panel_note = ("K = 86, both of the target's haplotypes held out, with its own "
                  "chain and its own shares" if a.loo else "K = 88, nothing held out")
    print(f"⚠ ARM {'1' if a.loo else '2'}: {panel_note}. The two arms carry different "
          f"ceilings and different utility functions, so no figure crosses between "
          f"them (standing constraint 8).\n")

    if a.mc:
        monte_carlo(a.chroms[0], tids[0], max(a.taus), sh[a.chroms[0]], a.mc)
        print()

    rows = []
    t0 = time.time()
    for tid in tids:
        for tau in a.taus:
            for c in a.chroms:
                r = (chrom_variance_loo(c, tid, tau, sh[c]) if a.loo
                     else chrom_variance(c, tid, tau, sh[c]))
                rows.append(r)
                print(f"    {tid} tau={tau:g} {c}: var_q_u={r['var_q_u']:.6e} "
                      f"cov_share={r['cov_share']:.4f} "
                      f"sd_differ={math.sqrt(r['var_q_differ']):,.1f} "
                      f"({time.time() - t0:.0f}s)", flush=True)

    # The figures this run is meant to predict, from the logs named in the
    # docstring. Quoted as constants so the comparison is explicit and auditable.
    measured = ({"measured chr21-at-share-1 movement "
                 "(logs/B_beta_concentration.md)": 3.6368e-06}
                if a.share1 else
                {"measured genome-wide movement (logs/B_exact_utility.md)": 8.345e-08}
                if len(a.chroms) == 22 else None)
    report(rows, a.taus, scope, measured)

    if a.json:
        a.json.write_text(json.dumps(
            {"scope": scope, "arm": "arm2", "share1": bool(a.share1),
             "results": rows}, indent=1) + "\n")
        print(f"\nJSON -> {a.json}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
