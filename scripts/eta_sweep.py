#!/usr/bin/env python
"""The RESPONSE CURVE: `E_Q[u]` against `eta` directly, far past any usable `tau`.

⚠⚠ EXPLORATORY AND NOT PRE-REGISTERED, and it is not a proposal to release anything.
Dylan, 2026-10-23: stop treating `tau` as the knob and sweep `eta` itself, to find out
whether the mechanism has ANY regime in which the tilt moves the release. Every row
here is a CHARACTERISATION of the mechanism family, never an operating point.

WHY `eta` AND NOT `tau`, AND IT IS NOT A DODGE. The two are a bijection,
`eta = arctanh(tau)` and `tau = tanh(eta)`, so sweeping `eta` IS sweeping `tau` --
there is no extra freedom hiding in the reparameterisation. But `tau` cannot CARRY
the high end: `tanh(eta)` rounds to exactly 1.0 in float64 beyond `eta` ~ 19, and
`calibration.eta_tau` rightly refuses `tau` = 1, which is the vacuous guarantee
(TV <= 1 holds between any two distributions). So above `eta` ~ 19 the parameter must
be `eta`, and the guarantee must be reported in LOG SPACE:

    log10(1 - tau) = log10(2) - 2*eta/log(10) + O(e^-4eta)

⚠⚠ THE GUARANTEE IS STILL REPORTED, ALWAYS. `1 - tau` and `epsilon_tau = 2*eta` are
printed on every row precisely because an `eta` table with no `tau` column would read
as a mechanism with a privacy parameter while actually being a mechanism without one.
At `eta` = 20 the theorem already permits an attacker to identify the target with
certainty. Dropping that column is the one thing this script will not do.

WHAT THE CURVE IS. `Q_eta(y) ∝ R_D(y) exp(eta * u(p,y))` is an exponential family in
`eta` with sufficient statistic `u`, so

    E_Qeta[u]  is monotone increasing in `eta`,   d E_Qeta[u] / d eta = Var_Qeta(u).

It runs from `E_R[u]` at `eta` = 0 (the untilted cohort mosaic, carrying zero target
information) up to `max_y u(p,y)` over the output support as `eta` -> infinity, where
the tilt degenerates to argmax. ⚠ That limit is a CHARACTERISATION, not a decoder:
standing constraint 3 forbids argmax as a release and the canonical decoder stays
`ffbs`. The limit is approached here by pushing `eta` high, never by running Viterbi.

So "does the mechanism have an effective regime" becomes "at what `eta` does the
curve leave its floor", which is exact and needs no sampling. ⚠ `Var_Qeta(u)` is
reported beside it because it IS the slope: a flat curve and a tiny variance are the
same statement, and the variance says how far the next unit of `eta` would move it.

THE THRESHOLDS THIS WAS BUILT TO BRACKET, predicted BEFORE the first run from the
repo's own inequalities (`docs/plan.md` R5, the reusable criterion):
  * buy ONE donor switch: needs `eta`/2 * (mass it controls) > the diversion cost,
    median 10.1529 nats on chr21, so `eta` > 20.31 with the WHOLE budget on it.
  * move a TYPICAL position: needs (`eta`/2) * `beta_v` ~ O(1), so `eta` ~ 2/`beta_v`
    -- about 6.65e5 on chr21 at `beta_share` = 1, and 4.43e7 genome-wide.
  * ⚠ The onset is set by `max(beta_v)`, NOT the mean: `max/mean` = 77.0 on this
    cohort, so the most heavily weighted positions start moving about 77x earlier,
    at `eta` ~ 2/max(`beta_v`) = 8.6e3 on chr21 at chromosome scope. Expect a
    gradual onset from ~1e4, not a step at the mean-based threshold.
  * `psi` overflows float64 at (`eta`/2) * max(`beta_v`) > 709 -- the WALL, and it
    is a property of `beta`'s scope rather than of the science: `eta` = 6.12e6 at
    chr21 chromosome scope, 4.08e8 at chr21's genome share. `main` refuses an `eta`
    past the wall UP FRONT rather than after the expensive passes, and `one_eta`
    asserts finiteness anyway.

⚠ THE REGIME VARIABLE IS `eta` * `beta_v`, NOT `eta`. That is what makes the toy gate
meaningful: the toy's `beta_v` = 1/7 is five orders of magnitude LARGER than a real
chromosome's, so toy `eta` = 2000 sits at the same `eta`*`beta_v` = 285.7 as chr21 at
`eta` = 1.23e6. The gate therefore covers the sweep's own regime, and it reports both
axes so that correspondence is checkable rather than asserted.

    python scripts/eta_sweep.py --validate
    python scripts/eta_sweep.py --chrom chr21 --target HG00096 --json out.json
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
sys.path.insert(0, str(REPO / "scripts"))

import numpy as np                                        # noqa: E402

from mechanism import calibration as C                    # noqa: E402
from mechanism import exact, ffbs, hmm, toy, utility       # noqa: E402

from exact_utility import AUTOSOMES                        # noqa: E402
from exact_utility import alt_totals, rarity_totals        # noqa: E402
from exact_utility import rarity_weights                   # noqa: E402
from exact_utility import shares as exact_utility_shares   # noqa: E402
from exact_variance import haplotype_variance              # noqa: E402
from exact_kl_stable import extremal, kl_stable            # noqa: E402

# log-spaced, and deliberately spanning both predicted thresholds and the pinned
# grid's own top (`eta_tau`(0.99) = 2.6467, `eta_tau`(0.99999) = 6.1030)
DEFAULT_ETAS = (0.0, 1.0, 2.6466524124, 6.1030338228, 20.0, 60.0, 200.0,
                1e3, 1e4, 1e5, 1e6, 3e6)


def build_betas(d, ic, p_chain, weights: str, *, share=None, chrom=None,
                alt_tot=None):
    """The per-haplotype `beta` pair for one chromosome, under either weighting.

    `production` -- `w_v` = 1/support(v) over every chain position, normalised to this
    chromosome's decision-J share. Target-INDEPENDENT, so both chains share one array.

    `target_alt` -- uniform over the positions where THAT haplotype carries a NON-REF
    allele, normalised by the GENOME-WIDE count of them (`alt_tot`), so the 22
    chromosomes' weights still sum to one budget (standing constraint 4). The
    expression is lifted from `scripts/alt_beta_sweep.py` so there is one definition.
    ⚠ TARGET-DEPENDENT. Theorem 1 still holds -- `beta` >= 0, sums to 1, `u` stays in
    [0,1] -- but the utility function is now a function of the target, so it is NOT
    interchangeable with the production one and the two must never share a column.
    """
    if weights in ("uniform", "production"):
        # "production" is the legacy name for inv_support; "uniform" is the
        # released tool's weighting as of 2026-10-24 (argued in docs/terms.md)
        w = "uniform" if weights == "uniform" else "inv_support"
        b = utility.beta_v(np.load(d / "support.npy"), in_chain=ic, weights=w,
                           **({} if share is None else {"share": share,
                                                        "chrom": chrom})).beta
        return (b, b)
    if weights == "target_alt":
        if alt_tot is None:
            raise ValueError("target_alt needs alt_tot, the genome-wide non-ref count")
        return tuple((p_chain[:, k] > 0) / float(alt_tot[k]) for k in (0, 1))
    if weights == "target_rarity":
        # `w_v = -log f_v`, f_v the cohort frequency of the TARGET's allele. Downweights
        # common alleles continuously instead of by the ref/non-ref cut, gives ZERO to
        # an allele every donor carries, and gives zero to an allele NO donor carries
        # (unmatchable, so budget there is wasted). See exact_utility.rarity_weights.
        if alt_tot is None:
            raise ValueError("target_rarity needs the genome-wide rarity total")
        hap_all = np.asarray(np.load(d / "haplotypes.npy", mmap_mode="r"))[:, ic]
        return tuple(rarity_weights(hap_all, p_chain[:, k]) / float(alt_tot[k])
                     for k in (0, 1))
    raise ValueError(f"unknown weights {weights!r}")


def genome_chrom(args) -> dict:
    """One chromosome of the GENOME-WIDE curve, at its decision-J share.

    A worker body, so it takes one tuple and never raises into the pool: a fan-out
    that loses a chromosome silently would produce a complete-looking curve over 21.
    """
    chrom, target, share, etas, with_kl, weights, alt_tot = args
    try:
        d = REPO / "data" / "cohort" / chrom
        ic = np.load(d / "in_chain.npy").astype(bool)
        hap = np.load(d / "haplotypes.npy")[:, ic]
        m = hmm.from_arrays(hap, np.load(d / "genetic_pos.npy")[ic],
                            np.load(d / "positions.npy")[ic], cohort_id=chrom)
        p = np.load(REPO / "data" / "targets" / target / chrom / "path.npy")[ic]
        betas = build_betas(d, ic, p, weights, share=share, chrom=chrom,
                            alt_tot=alt_tot)
        wall = psi_wall(betas)
        over = [e for e in etas if e >= wall]
        if over:
            raise ValueError(f"{chrom}: eta {over} past its wall {wall:.3e}")
        ceil = utility.diploid_utility(utility.u_ceiling(p[:, 0], betas[0]),
                                       utility.u_ceiling(p[:, 1], betas[1]))
        base = ffbs.baseline_log_norm(m)
        rows = []
        for eta in etas:
            # ⚠ NAME THE `eta` IN ANY FAILURE. The forward pass can dead-end at a
            # zero-distance step once the posterior concentrates (measured 2026-10-23
            # on chr7 and chr10), and that is eta-DEPENDENT, so a traceback without
            # the eta does not say whether the curve died at its second point or its
            # eighth -- which is the only thing that decides what is reportable.
            try:
                r = one_eta(m, hap, p, betas, eta, with_kl=with_kl, base=base)
            except Exception as e:
                raise RuntimeError(
                    f"{chrom} FAILED AT eta={eta:g} (rows up to "
                    f"{rows[-1]['eta']:g} succeeded): {type(e).__name__}: {e}"
                ) from e
            r["chrom"] = chrom
            rows.append(r)
        return {"chrom": chrom, "share": share, "t_chain": int(ic.sum()),
                "u_ceiling": ceil, "weights": weights,
                "beta_sum": float(sum(b.sum() for b in betas) / 2.0),
                "beta_max": max(float(b.max()) for b in betas),
                "wall": wall, "rows": rows}
    except Exception as e:                  # a worker must report, never vanish
        import traceback
        return {"chrom": chrom, "error": f"{type(e).__name__}: {e}",
                "traceback": traceback.format_exc()[-1200:]}


def psi_wall(betas) -> float:
    """The largest `eta` whose potentials stay finite: (`eta`/2)*max(`beta_v`) < 709.

    A property of `beta`'s SCOPE, not of the mechanism. Returned so the sweep can
    refuse an out-of-range `eta` before spending a forward pass on it. Takes the
    per-haplotype pair and uses the WORST of the two.
    """
    b = betas if isinstance(betas, (tuple, list)) else (betas, betas)
    return 2.0 * 709.0 / max(float(x.max()) for x in b)


def log10_one_minus_tau(eta: float) -> float:
    """`log10(1 - tanh(eta))`, computed in log space because float64 cannot hold it.

    `1 - tanh(x) = 2/(e^{2x} + 1)`, so the log is exact to the `log1p` term and the
    asymptote `log10(2) - 2x/ln10` is what survives at large `x`.
    """
    if eta == 0.0:
        return 0.0
    return math.log10(2.0) - (2.0 * eta) / math.log(10.0) - math.log10(
        1.0 + math.exp(-2.0 * eta)) if eta < 350 else (
        math.log10(2.0) - (2.0 * eta) / math.log(10.0))


def tau_label(eta: float) -> str:
    """The equivalent `tau`, printable at any `eta`."""
    if eta == 0.0:
        return "0"
    g = log10_one_minus_tau(eta)
    if g > -15:
        return f"{math.tanh(eta):.15f}"
    return f"1 - 1e{g:.4g}"


def enumerated_moments(base, u_vec: np.ndarray, eta: float) -> tuple[float, float]:
    """`E_Q[u]` and `Var_Q(u)` by direct enumeration, in LOG SPACE.

    ⚠ WHY NOT `exact.release_distribution`. That computes `prob * exp(eta * u)`
    directly, and `u` reaches 1 by construction, so it overflows float64 at
    `eta` ~ 709 WHATEVER `beta` is -- a far lower wall than the sampler's
    709/max(`beta_v`), and therefore the binding constraint on how far this gate
    could check. Measured: it returned nan at `eta` = 2000 and the first version of
    this gate compared nan to nan and called it agreement.

    Log-sum-exp removes the wall entirely and leaves the SAME quantity, so the
    reference now reaches as far as the sampler does. `exact.py` is deliberately not
    touched: it is load-bearing for the signed pre-registration's toy gates, and a
    local reference that owes nothing to the code under test is the better gate
    anyway.
    """
    lw = np.log(base.prob) + eta * u_vec
    lw -= lw.max()
    q = np.exp(lw)
    q /= q.sum()
    e = float(q @ u_vec)
    return e, float(q @ (u_vec - e) ** 2)


def one_eta(m, hap: np.ndarray, p: np.ndarray, betas, eta: float,
            *, with_kl: bool = False, base: float | None = None) -> dict:
    """Exact `E_Q[u]`, `Var_Q(u)` and `E[differ]` for one `eta`, diploid.

    Two forward + state_marginals passes, one per haplotype chain, each at `eta`/2 as
    the diploid aggregation requires. No sampling anywhere.

    ⚠ `betas` is a PER-HAPLOTYPE pair. The production weighting gives both chains the
    same array, but `target_alt` does NOT: its support is the positions where THAT
    haplotype carries a non-ref allele, and the two haplotypes differ there (2,493,299
    against 2,501,582 for HG00096). Sharing one array across the chains would weight
    each haplotype by the other's non-ref set.
    """
    h = eta / 2.0
    means, vars_, diffs, called = [], [], [], []
    nonref, n_nonref = [], []
    # ⚠⚠ PMI IS PANMIXER'S OWN PRIVACY METRIC and we can report it in CLOSED FORM.
    # PMI(p, y) = log[P(y|p)/P(y)] = log[Q_p(y)/R_D(y)] = eta*u(p,y) - log Z_p,
    # because Z_p = E_{R_D}[exp(eta u)] by construction. So
    #     E_Q[PMI] = eta*E_Q[u] - log Z_p = KL(Q_p || R_D),
    # the divergence from the release distribution to the UNTILTED baseline. They
    # estimate it blockwise (their Supplementary Lemma 1); this is exact from the
    # same forward pass as e_q_u. ⚠ It is a DIFFERENT question from acc_ceiling:
    # acc_ceiling asks whether an attacker can tell WHICH INPUT produced an output;
    # PMI asks how statistically DEPENDENT the output is on the input. Report both.
    log_z = 0.0 if base is not None else None
    for col in (0, 1):
        tgt = p[:, col]
        beta = betas[col]
        psi = utility.tilt_potentials(m.allowed, hap, tgt, beta, h)
        if not np.isfinite(psi).all():
            raise AssertionError(
                f"psi overflowed at eta={eta:g}: max (eta/2)*beta_v = "
                f"{h * float(beta.max()):.1f} nats, and exp() dies near 709. The "
                f"curve cannot be continued past this eta at this beta scope.")
        fwd = ffbs.forward(m, psi)
        gam = ffbs.state_marginals(m, fwd, psi)
        match = (hap == tgt[None, :]) & (tgt[None, :] >= 0)
        if base is not None:
            log_z += ffbs.log_z_p(fwd, base, h, budget=h * float(beta.sum()))
        r = haplotype_variance(m, fwd.alpha, gam, match, beta)
        means.append(r["mean"])
        vars_.append(r["var"])
        # ⚠⚠ THE ONE UTILITY NUMBER THAT IS COMPARABLE ACROSS WEIGHTINGS. `e_q_u` is
        # beta-weighted, so production's and target_alt's are DIFFERENT quantities
        # (constraint 8). `e_differ` counts every called cell, and is dominated by the
        # ~19.6M REF positions a cohort mosaic matches by chance. Neither answers
        # "how many of the target's NON-REF alleles did the release recover", which is
        # the biologically interesting count and is weighting-INDEPENDENT. This is it:
        # the expected number of non-ref positions where the release matches.
        pm = (gam * match.T).sum(axis=1)
        nonref.append(float(pm[tgt > 0].sum()))
        n_nonref.append(int((tgt > 0).sum()))
        # E[differ] over CALLED cells, the same quantity ffbs.expected_mismatches gives
        diffs.append(float(ffbs.expected_mismatches(m, hap, tgt, gam)))
        called.append(int((tgt >= 0).sum()))
        del psi, fwd, gam, match
    # the diploid utility is the MEAN of the two chains, so Var scales by 1/4
    e_q_u = 0.5 * (means[0] + means[1])
    var = 0.25 * (vars_[0] + vars_[1])
    out = {"eta": eta, "epsilon": 2.0 * eta, "tau": tau_label(eta),
           "log10_one_minus_tau": log10_one_minus_tau(eta),
           "e_q_u": e_q_u, "var_q_u": var, "sd_q_u": math.sqrt(max(var, 0.0)),
           "e_differ": diffs[0] + diffs[1], "called_cells": called[0] + called[1],
           "e_match_nonref": nonref[0] + nonref[1],
           "n_nonref": n_nonref[0] + n_nonref[1],
           "u_hap": means, "var_hap": vars_}
    if log_z is not None:
        # the diploid tilt is exp(eta*u_mean) = exp((eta/2)u_0)*exp((eta/2)u_1), so
        # log Z is the SUM of the two chains' and PMI uses the full eta on u_mean
        out["log_z_p"] = log_z
        out["pmi"] = eta * e_q_u - log_z
    if with_kl:
        # the SAME kl_stable the genome-wide ceiling uses, driven by eta
        k = kl_stable(m, hap, betas, p, extremal(p), eta=eta)
        tv = min(1.0, math.sqrt(max(k["kl"], 0.0) / 2.0))
        out |= {"kl_exact": k["kl"], "tv_pinsker": tv,
                "acc_ceiling": 0.5 * (1.0 + tv)}
    return out


# --------------------------------------------------------------------- the gate

def validate() -> int:
    """⚠ THE GATE, AND IT MUST BE ABLE TO FAIL. On the hand-written toy the release
    distribution is ENUMERABLE, so `E_Q[u]` and `Var_Q(u)` have closed forms that owe
    nothing to the forward pass. This checks the sweep's route against them AT THE
    EXTREME `eta` THE SWEEP IS FOR -- the pinned-grid validation in
    `scripts/exact_variance.py` only ever ran at `tau` <= 0.99, i.e. `eta` <= 2.65,
    and a route that is right there can still die at `eta` = 1e6.
    """
    t = toy.toy_cohort()
    m = hmm.from_arrays(t.haplotypes, t.genetic_pos, t.positions, cohort_id="toy")
    base = exact.enumerate_baseline(t.haplotypes, t.genetic_pos)
    beta = np.full(t.T, 1.0 / t.T)
    tgt = t.haplotypes[0].copy()
    u_vec = exact.utility_vector(base, tgt, beta)
    # ⚠ The toy's own wall. `validate` passes the FULL eta to tilt_potentials (not
    # eta/2 -- the toy is haploid here), so its wall is 709/max(beta_v), and the
    # first version of this gate asked for eta = 1e4 and 1e6, sailed past it, and
    # compared nan to nan. A gate that reads nan as agreement is worse than none.
    # the SAMPLER's wall; the log-space reference above has none
    wall = 709.0 / float(beta.max())
    mx = float(beta.max())
    print("TOY GATE. Enumerated support against the forward/marginals route.\n")
    print(f"K = {m.K}, T = {m.T}, support = {len(u_vec)} paths, "
          f"u in [{u_vec.min():.6f}, {u_vec.max():.6f}], max(beta_v) = {mx:.4f}")
    print(f"⚠ THE REGIME VARIABLE IS `eta`*`beta_v`, so the toy's small T covers the "
          f"sweep's range: toy `eta` = 2000 sits where chr21 sits at `eta` ~ 1.2e6.")
    print(f"⚠ The toy's OWN overflow wall is `eta` = {wall:.1f} "
          f"(709/max(beta_v)); nothing is compared past it.\n")
    print(f"| toy `eta` | `eta`*max(`beta_v`) | `E_Q[u]` enumerated | route | "
          f"rel. diff | `Var` enumerated | route | rel. diff |")
    print("|---|---|---|---|---|---|---|---|")
    # ⚠⚠ THE TWO QUANTITIES HAVE DIFFERENT VALID RANGES, MEASURED HERE RATHER THAN
    # ASSUMED. `E_Q[u]` agrees to 1e-16 at every `eta`. `Var_Q(u)` agrees only while
    # it is above a FLOOR near 1e-33: the centred-moment recursion in
    # `haplotype_variance` cancels, so once the posterior collapses onto one path it
    # reports a constant ~7.7e-34 while the true variance keeps shrinking past
    # 1e-270. Checking them under one tolerance would either fail on a harmless
    # floor or pass by loosening the tolerance until nothing was tested.
    VAR_FLOOR_PROBE = 1e-28
    worst_mean = worst_var = 0.0
    floored = []
    for eta in (0.0, 1.0, 6.1030338228, 20.0, 200.0, 2000.0, 4400.0):
        assert eta < wall, f"validation eta {eta} is past the toy's wall {wall}"
        e_ref, v_ref = enumerated_moments(base, u_vec, eta)
        psi = utility.tilt_potentials(m.allowed, t.haplotypes, tgt, beta, eta)
        fwd = ffbs.forward(m, psi)
        gam = ffbs.state_marginals(m, fwd, psi)
        match = (t.haplotypes == tgt[None, :]) & (tgt[None, :] >= 0)
        r = haplotype_variance(m, fwd.alpha, gam, match, beta)
        assert np.isfinite([e_ref, v_ref, r["mean"], r["var"]]).all(), (
            f"eta={eta:g}: a non-finite value entered the comparison, so the gate "
            f"would be comparing nan to nan")
        de = abs(r["mean"] - e_ref) / max(abs(e_ref), 1e-300)
        worst_mean = max(worst_mean, de)
        if v_ref > VAR_FLOOR_PROBE:
            dv = abs(r["var"] - v_ref) / v_ref
            worst_var = max(worst_var, dv)
            dvs = f"{dv:.2e}"
        else:
            floored.append((eta, v_ref, r["var"]))
            dvs = "n/a (below the floor)"
        print(f"| {eta:g} | {eta * mx:.2f} | {e_ref:.12f} | {r['mean']:.12f} | "
              f"{de:.2e} | {v_ref:.6e} | {r['var']:.6e} | {dvs} |")
    print(f"\nworst relative difference, `E_Q[u]`: {worst_mean:.3e}  (all `eta`)")
    print(f"worst relative difference, `Var_Q(u)` above {VAR_FLOOR_PROBE:g}: "
          f"{worst_var:.3e}")
    if floored:
        print(f"\n⚠⚠ **`var_q_u` HAS A FLOOR AND IT ERRS IN THE FLATTERING "
              f"DIRECTION.** Where the true variance falls below ~1e-28 the route "
              f"reports a constant instead:")
        for eta, v_ref, v_got in floored:
            print(f"  `eta` = {eta:g}: true {v_ref:.3e}, route {v_got:.3e} "
                  f"(overstated by 10^{math.log10(v_got / max(v_ref, 1e-300)):.0f})")
        print(f"  `Var` is the CURVE'S SLOPE, so an overstated floor overstates how "
              f"much the next unit of `eta` would buy. **Never quote `var_q_u` "
              f"below ~1e-28**; `e_q_u` itself is unaffected and exact throughout.")
        assert all(v_got >= v_ref for _, v_ref, v_got in floored), (
            "the var floor runs BELOW the true value somewhere, so it is not a "
            "simple floor and the direction of the error is not known")
    # ⚠ AND THE WALL ITSELF MUST BE REAL: past it the potentials MUST go non-finite,
    # so the finiteness guard in `one_eta` has something to catch. Asserting this is
    # what stops the wall being a number nobody checked.
    over = utility.tilt_potentials(m.allowed, t.haplotypes, tgt, beta, wall * 1.5)
    assert not np.isfinite(over).all(), (
        f"psi stayed finite at eta={wall*1.5:.1f}, past the computed wall "
        f"{wall:.1f} -- the wall is wrong and one_eta's guard would never fire")
    print(f"✅ and past the wall (`eta` = {wall*1.5:.0f}) the potentials DO go "
          f"non-finite, so `one_eta`'s finiteness guard has a real defect to catch.")
    # ⚠ The monotonicity and the limit are PROPERTIES OF THE FAMILY, so assert them
    # here rather than trusting the table to be read.
    seq = []
    for eta in (0.0, 1.0, 10.0, 100.0, 1000.0, 4400.0):
        seq.append(enumerated_moments(base, u_vec, eta)[0])
    for a, b in zip(seq, seq[1:]):
        assert b >= a - 1e-12, f"E_Q[u] is not monotone in eta: {a} -> {b}"
    assert abs(seq[-1] - u_vec.max()) < 1e-9, (
        f"E_Q[u] does not approach max_y u as eta grows: {seq[-1]} vs {u_vec.max()}")
    print(f"✅ monotone in `eta`, and the limit reaches max_y u = {u_vec.max():.12f} "
          f"(the argmax limit, a characterisation and NOT a decoder).")
    if max(worst_mean, worst_var) > 1e-8:
        raise SystemExit(
            f"the sweep's route does NOT reproduce the enumerated curve (mean "
            f"{worst_mean:.3e}, var {worst_var:.3e}). No swept number may be "
            f"reported until it does.")
    print(f"\n✅ `e_q_u` reproduces the enumerated curve to {worst_mean:.2e}, and "
          f"`var_q_u` to {worst_var:.2e} above its floor, at every `eta` up to "
          f"`eta`*max(`beta_v`) = {4400.0 * mx:.0f} -- the regime chr21 reaches at "
          f"`eta` ~ {2 * 4400.0 * mx / 2.3156e-04:.2e} (the sweep tilts at `eta`/2).")
    return 0


# --------------------------------------------------------------------- the sweep

def run_genome(a) -> int:
    """The GENOME-WIDE curve: all 22 autosomes at decision-J shares, summed.

    ⚠ THIS IS THE ONE THAT MATTERS, because the release IS genome-wide and the
    single-chromosome curve CANNOT be rescaled into it: the `eta` * `beta_v` scaling
    is exact WITHIN a chromosome (measured to 8-9 digits) but chr1 sits 1.0-1.4
    points of ceiling fraction above chr21 at matched `eta` * mean(`beta_v`), so the
    chromosomes are not one curve and the composition has to be run.
    """
    from concurrent.futures import ProcessPoolExecutor, as_completed
    sh, T, tot = exact_utility_shares()
    print(f"⚠⚠ EXPLORATORY, NOT PRE-REGISTERED, AND NOT AN OPERATING POINT.")
    print(f"THE GENOME-WIDE RESPONSE CURVE. {a.target}, ARM 2, 22 autosomes, "
          f"T = {tot:,}, shares by decision J summing to {sum(sh.values()):.6f}, "
          f"diploid, EXACT via state marginals.")
    print(f"{len(a.etas)} `eta` x 22 chromosomes, {a.workers} workers.\n")
    alt_tot = (alt_totals(a.target) if a.weights == "target_alt"
               else rarity_totals(a.target) if a.weights == "target_rarity"
               else None)
    if alt_tot is not None:
        # ⚠ the normaliser means different things per weighting, so do NOT print one
        # label for both: target_alt's is a COUNT of non-ref positions, rarity's is a
        # SUM of -log f and is not a count of anything.
        if a.weights == "target_alt":
            print(f"`target_alt`: {alt_tot[0]:,.0f} and {alt_tot[1]:,.0f} non-ref "
                  f"chain positions per haplotype, so beta_v = "
                  f"{1.0/alt_tot[0]:.4e} there and 0 elsewhere.\n")
        else:
            print(f"`target_rarity`: sum(-log f) = {alt_tot[0]:,.1f} and "
                  f"{alt_tot[1]:,.1f} per haplotype (a WEIGHT TOTAL, not a count); "
                  f"beta_v = -log(f_v) / that, bounded by log(K)/total = "
                  f"{__import__('math').log(88)/alt_tot[0]:.4e}.\n")
    jobs = [(c, a.target, sh[c], list(a.etas), a.with_kl, a.weights, alt_tot)
            for c in AUTOSOMES]
    out, t0 = {}, time.time()
    with ProcessPoolExecutor(max_workers=a.workers) as ex:
        futs = {ex.submit(genome_chrom, j): j[0] for j in jobs}
        for n, f in enumerate(as_completed(futs), 1):
            r = f.result()
            out[r["chrom"]] = r
            tag = f"❌ {r['error']}" if "error" in r else f"✅ T = {r['t_chain']:,}"
            print(f"  [{n}/22] {r['chrom']}: {tag} "
                  f"[{(time.time()-t0)/60:.0f} min]", flush=True)
    errs = [r for r in out.values() if "error" in r]
    # ⚠ completeness asserted from the RESULT SET, not from the loop's bookkeeping
    missing = [c for c in AUTOSOMES if c not in out or "error" in out[c]]
    if missing:
        for r in errs[:3]:
            print(f"\n❌ {r['chrom']}\n{r['traceback']}")
        raise SystemExit(f"INCOMPLETE: {len(missing)} chromosome(s) missing "
                         f"({missing[:5]}). A curve over 21 autosomes is not the "
                         f"genome-wide curve and must not be reported as one.")

    ceil = math.fsum(out[c]["u_ceiling"] for c in AUTOSOMES)
    agg = []
    for i, eta in enumerate(a.etas):
        cells = [out[c]["rows"][i] for c in AUTOSOMES]
        assert all(abs(x["eta"] - eta) < 1e-9 for x in cells), "eta misalignment"
        row = {"eta": eta, "epsilon": 2.0 * eta, "tau": tau_label(eta),
               "log10_one_minus_tau": log10_one_minus_tau(eta),
               "e_q_u": math.fsum(x["e_q_u"] for x in cells),
               "var_q_u": math.fsum(x["var_q_u"] for x in cells),
               "e_differ": math.fsum(x["e_differ"] for x in cells),
               "e_match_nonref": math.fsum(x["e_match_nonref"] for x in cells),
               "pmi": math.fsum(x.get("pmi", 0.0) for x in cells),
               "log_z_p": math.fsum(x.get("log_z_p", 0.0) for x in cells),
               "n_nonref": sum(x["n_nonref"] for x in cells),
               "called_cells": sum(x["called_cells"] for x in cells),
               "u_ceiling": ceil}
        row["e_q_u_over_ceiling"] = row["e_q_u"] / ceil if ceil else float("nan")
        if a.with_kl:
            row["kl_exact"] = math.fsum(x["kl_exact"] for x in cells)
            tv = min(1.0, math.sqrt(max(row["kl_exact"], 0.0) / 2.0))
            row |= {"tv_pinsker": tv, "acc_ceiling": 0.5 * (1.0 + tv)}
        agg.append(row)

    base = agg[0]["e_q_u"] if agg[0]["eta"] == 0.0 else None
    print(f"\n`u_ceiling` = {ceil:.9f} (summed over 22 autosomes).")
    print(f"\nMetrics are ROWS, `eta` is the COLUMN.\n")
    print("| metric | " + " | ".join(f"{r['eta']:g}" for r in agg) + " |")
    print("|" + "---|" * (len(agg) + 1))
    print("| equivalent `tau` | " + " | ".join(r["tau"] for r in agg) + " |")
    print("| `epsilon_tau` | " + " | ".join(f"{r['epsilon']:.4g}" for r in agg) + " |")
    print("| `e_q_u` * | " + " | ".join(f"{r['e_q_u']:.12f}" for r in agg) + " |")
    print("| `e_q_u_over_ceiling` | " + " | ".join(
        f"{r['e_q_u_over_ceiling']:.9f}" for r in agg) + " |")
    if base is not None:
        gap = ceil - base
        print("| movement from `eta` = 0 | " + " | ".join(
            f"{r['e_q_u'] - base:+.4e}" for r in agg) + " |")
        print("| % of the recoverable gap | " + " | ".join(
            f"{100*(r['e_q_u']-base)/gap:.3f}%" for r in agg) + " |")
    print("| `var_q_u` | " + " | ".join(f"{r['var_q_u']:.4e}" for r in agg) + " |")
    print("| `e_differ` cells | " + " | ".join(
        f"{r['e_differ']:,.1f}" for r in agg) + " |")
    print("| `e_match_nonref` * [COMPARABLE across weightings] | " + " | ".join(
        f"{r['e_match_nonref']:,.1f}" for r in agg) + " |")
    if any(r.get("pmi") for r in agg):
        print("| `pmi` * [PanMixer's OWN privacy metric, exact] | " + " | ".join(
            f"{r['pmi']:.4e}" for r in agg) + " |")
        print("| `log_z_p` | " + " | ".join(
            f"{r['log_z_p']:.4e}" for r in agg) + " |")
    if base is not None:
        n0 = agg[0]["e_match_nonref"]
        print("| non-ref alleles RECOVERED vs `eta_tau` = 0 | " + " | ".join(
            f"{r['e_match_nonref'] - n0:+,.1f}" for r in agg) + " |")
    if a.with_kl:
        print("| `kl_exact` | " + " | ".join(f"{r['kl_exact']:.4e}" for r in agg) + " |")
        print("| `acc_ceiling` * | " + " | ".join(
            f"{r['acc_ceiling']:.9f}" for r in agg) + " |")
    print("| `p_succ_bound` * | " + " | ".join(
        ("0.5" if r["eta"] == 0 else
         (f"{(1+math.tanh(r['eta']))/2:.9f}" if r["log10_one_minus_tau"] > -15
          else "1 - 5e" + f"{r['log10_one_minus_tau']:.4g}")) for r in agg) + " |")
    print(f"\n⚠⚠ **READ THE MOVEMENT, NOT THE LEVEL.** `e_q_u` at `eta` = 0 is the "
          f"untilted cohort mosaic carrying ZERO target information.")
    print(f"⚠ `var_q_u` is summed over chromosomes, which is correct because the 22 "
          f"chains are independent given the shares; it FLOORS near 1e-33 per "
          f"chromosome and must not be quoted below ~1e-28 (docs/terms.md).")
    print(f"⚠ Standing constraint 3: the `eta` -> infinity limit is argmax, NOT a "
          f"release procedure. Constraint 5: one release per genome; this is an "
          f"evaluation device.")
    if a.json:
        a.json.write_text(json.dumps(
            {"target": a.target, "arm": "arm2", "scope": "genome", "T": tot,
             "shares": sh, "u_ceiling": ceil, "exploratory": True,
             "aggregate": agg,
             "per_chrom": {c: {k: v for k, v in out[c].items() if k != "rows"}
                           for c in AUTOSOMES}}, indent=1) + "\n")
        print(f"\nJSON -> {a.json}")
    return 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--validate", action="store_true")
    ap.add_argument("--genome", action="store_true",
                    help="all 22 autosomes at decision-J shares, summed")
    ap.add_argument("--workers", type=int, default=11)
    ap.add_argument("--chrom", default="chr21")
    ap.add_argument("--target", default="HG00096")
    ap.add_argument("--etas", type=float, nargs="+", default=list(DEFAULT_ETAS))
    ap.add_argument("--share", type=float, default=None,
                    help="beta_share; default None = chromosome-scope (sums to 1), "
                         "which is 66.6x MORE generous per position than genome-wide")
    ap.add_argument("--weights",
                    choices=("uniform", "production", "target_alt",
                             "target_rarity"),
                    default="uniform",
                    help="production = w_v 1/support(v) over the whole chain; "
                         "target_alt = uniform over the target haplotype's NON-REF "
                         "positions only (target-dependent; a DIFFERENT utility "
                         "function, never pooled with production)")
    ap.add_argument("--with-kl", action="store_true",
                    help="also compute kl_exact/acc_ceiling against D1_extremal "
                         "(doubles the passes)")
    ap.add_argument("--json", type=Path, default=None)
    a = ap.parse_args(argv)

    if a.validate:
        return validate()
    if a.genome:
        return run_genome(a)

    d = REPO / "data" / "cohort" / a.chrom
    ic = np.load(d / "in_chain.npy").astype(bool)
    hap = np.load(d / "haplotypes.npy")[:, ic]
    m = hmm.from_arrays(hap, np.load(d / "genetic_pos.npy")[ic],
                        np.load(d / "positions.npy")[ic], cohort_id=a.chrom)
    kw = {} if a.share is None else {"share": a.share, "chrom": a.chrom}
    p = np.load(REPO / "data" / "targets" / a.target / a.chrom / "path.npy")[ic]
    betas = build_betas(d, ic, p, a.weights,
                        alt_tot=(alt_totals(a.target) if a.weights == "target_alt"
                                 else rarity_totals(a.target)
                                 if a.weights == "target_rarity" else None), **kw)
    beta = betas[0]
    ceil = utility.diploid_utility(utility.u_ceiling(p[:, 0], betas[0]),
                                   utility.u_ceiling(p[:, 1], betas[1]))

    print(f"⚠⚠ EXPLORATORY, NOT PRE-REGISTERED, AND NOT AN OPERATING POINT. The "
          f"signed grid ends at `tau` = 0.99 (`eta` = 2.6467).")
    print(f"THE RESPONSE CURVE. {a.target}, {a.chrom}, ARM 2, K = {m.K}, "
          f"T = {int(ic.sum()):,}, diploid, "
          f"`beta` scope {'chromosome (sums to 1)' if a.share is None else f'share {a.share:g}'}, "
          f"sum(beta) = {beta.sum():.6f}, mean(beta_v) = {beta.mean():.4e}.")
    print(f"`u_ceiling` = {ceil:.9f}. EXACT via state marginals; no sampling.")
    # ⚠ REFUSE AN OUT-OF-RANGE `eta` UP FRONT. The wall is set by max(beta_v) and is
    # a property of beta's scope; discovering it after an hour of forward passes is
    # the avoidable version of this failure.
    wall = psi_wall(betas)
    print(f"`psi` wall for this scope: `eta` < {wall:.3e} "
          f"(max(beta_v) = {beta.max():.4e}, max/mean = {beta.max()/beta.mean():.1f}).")
    over = [e for e in a.etas if e >= wall]
    if over:
        raise SystemExit(
            f"eta {over} is at or past this beta scope's overflow wall "
            f"{wall:.3e}: exp((eta/2)*max(beta_v)) would exceed float64. Lower the "
            f"eta, or widen the scope (a genome share puts the wall at ~4e8).")
    print()

    rows, t0 = [], time.time()
    for eta in a.etas:
        r = one_eta(m, hap, p, betas, eta, with_kl=a.with_kl,
                    base=ffbs.baseline_log_norm(m))
        r |= {"u_ceiling": ceil,
              "e_q_u_over_ceiling": r["e_q_u"] / ceil if ceil else float("nan")}
        rows.append(r)
        print(f"  eta={eta:<12g} e_q_u={r['e_q_u']:.12f} "
              f"var={r['var_q_u']:.4e} ({time.time()-t0:.0f}s)", flush=True)

    base = rows[0]["e_q_u"] if rows and rows[0]["eta"] == 0.0 else None
    print(f"\nMetrics are ROWS, `eta` is the COLUMN.\n")
    hdr = "| metric | " + " | ".join(f"{r['eta']:g}" for r in rows) + " |"
    print(hdr); print("|" + "---|" * (len(rows) + 1))
    print("| equivalent `tau` | " + " | ".join(r["tau"] for r in rows) + " |")
    print("| `epsilon_tau` | " + " | ".join(f"{r['epsilon']:.4g}" for r in rows) + " |")
    print("| `e_q_u` * | " + " | ".join(f"{r['e_q_u']:.12f}" for r in rows) + " |")
    print("| `e_q_u_over_ceiling` | " + " | ".join(
        f"{r['e_q_u_over_ceiling']:.9f}" for r in rows) + " |")
    if base is not None:
        print("| movement from `eta` = 0 | " + " | ".join(
            f"{r['e_q_u'] - base:+.4e}" for r in rows) + " |")
    print("| `var_q_u` | " + " | ".join(f"{r['var_q_u']:.4e}" for r in rows) + " |")
    print("| `e_differ` cells | " + " | ".join(
        f"{r['e_differ']:,.1f}" for r in rows) + " |")
    if any(r.get("pmi") for r in rows):
        print("| `pmi` * [PanMixer's OWN privacy metric, exact] | " + " | ".join(
            f"{r['pmi']:.4e}" for r in rows) + " |")
    if a.with_kl:
        print("| `kl_exact` | " + " | ".join(
            f"{r['kl_exact']:.4e}" for r in rows) + " |")
        print("| `acc_ceiling` * | " + " | ".join(
            f"{r['acc_ceiling']:.9f}" for r in rows) + " |")
    print("| `p_succ_bound` * | " + " | ".join(
        ("0.5" if r["eta"] == 0 else
         (f"{(1+math.tanh(r['eta']))/2:.9f}" if r["log10_one_minus_tau"] > -15
          else "1 - 5e" + f"{r['log10_one_minus_tau']:.4g}")) for r in rows) + " |")

    if base is not None:
        reach = max(r["e_q_u"] for r in rows)
        print(f"\n⚠⚠ **READ THE MOVEMENT, NOT THE LEVEL.** `e_q_u` = {base:.9f} is "
              f"reached at `eta` = 0, where the tilt is OFF and the release carries "
              f"zero target information. The highest `eta` swept reaches "
              f"{reach:.9f}, a movement of {reach - base:+.4e} over a privacy "
              f"parameter that has gone from `tau` = 0 to "
              f"`tau` = {rows[-1]['tau']}.")
    print(f"⚠ `beta` scope is "
          f"{'CHROMOSOME' if a.share is None else f'share {a.share:g}'}; a "
          f"chromosome-scope curve is far more generous per position than the "
          f"genome-wide release and the two must never be pooled.")
    print(f"⚠ Standing constraint 3: the `eta` -> infinity limit is argmax, which is "
          f"NOT a release procedure. These rows characterise the family.")

    if a.json:
        a.json.write_text(json.dumps(
            {"target": a.target, "chrom": a.chrom, "arm": "arm2",
             "T": int(ic.sum()), "K": int(m.K), "beta_sum": float(beta.sum()),
             "beta_mean": float(beta.mean()), "share": a.share,
             "u_ceiling": ceil, "exploratory": True, "rows": rows}, indent=1) + "\n")
        print(f"\nJSON -> {a.json}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
