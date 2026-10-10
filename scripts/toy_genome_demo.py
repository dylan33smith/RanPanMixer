#!/usr/bin/env python
"""The whole mechanism on a two-chromosome toy genome, enumerated exactly.

⚠ WHY THIS EXISTS AND WHAT IT IS NOT. `tests/test_diploid_gate.py` already gates the
diploid and multi-segment release — it checks that bounds HOLD. This is a different
object: a DEMONSTRATION that shows what the mechanism DOES. Every number here is
exact by enumeration over the full output support, with no sampling and no finite
differences, so the tilt's behaviour can be read off rather than inferred.

⚠ AND IT IS THE ONE PLACE THE TILT IS VISIBLE. On chr21 `beta_v` averages 3.27e-06
per position, so the tilt reweights each position by about one part in a million and
the released path moves while expected utility does not (deferred decision R5). Here
there are 9 chain positions, so `beta_v` averages 0.11 — about 34,000x larger — and
the same mechanism, running the same code, moves the output distribution visibly.
The contrast is the point: it separates "the mechanism does nothing" from "the budget
is spread too thin at genome scale", which are very different claims and only one of
them is true.

THE TOY GENOME. Two chromosomes, because two is enough to exercise every piece of the
genome-wide budget algebra that 22 would:
  chrA = the main toy, K = 4, T = 7, 96 output paths (a zero-distance step, a missing
         donor cell, a multi-allelic site, a donor that is another minus a cell)
  chrB = the small segment, K = 3, T = 2, 4 output paths (a wide step that drives
         switching to near-uniform)
Shares follow DECISION J — by variant count, `T_c / T_total` — so chrA takes 7/9 and
chrB takes 2/9. Nothing here is hand-tuned to make the result look better; the shares
are the rule the pre-registration pins.

    python scripts/toy_genome_demo.py [--draws 20000] [--json out.json]
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

import attacks as A  # noqa: E402
from mechanism import calibration as C  # noqa: E402
from mechanism import exact, ffbs, hmm, toy, utility  # noqa: E402

TAUS = (0.0, 0.05, 0.10, 0.25, 0.50, 0.90, 0.99)


def inv_support_beta(support: np.ndarray, share: float) -> np.ndarray:
    w = 1.0 / support.astype(np.float64)
    return w * (share / w.sum())


def build():
    """The two chromosomes, their baselines, and shares by decision J."""
    a = toy.toy_cohort()
    _, b = toy.toy_segments()
    t_total = a.T + b.T
    shares = (a.T / t_total, b.T / t_total)          # DECISION J: by variant count
    assert abs(sum(shares) - 1.0) < 1e-15
    bases = (exact.enumerate_baseline(a.haplotypes, a.genetic_pos),
             exact.enumerate_baseline(b.haplotypes, b.genetic_pos))
    betas = (inv_support_beta(a.support, shares[0]),
             inv_support_beta(b.support, shares[1]))
    return (a, b), bases, shares, betas


def joint_baseline(bases):
    return np.outer(bases[0].prob, bases[1].prob).ravel()


def release_and_u(bases, betas, targets, eta, ploidy):
    """Exact `Q` over the product support, `log Z`, and the composite `u`."""
    factors = []
    for h in range(ploidy):
        for c in (0, 1):
            factors.append(exact.Factor(bases[c], targets[h][c], betas[c],
                                        f"hap{h}chr{'AB'[c]}"))
    q, log_z, _ = exact.composite_release(factors, eta, ploidy=ploidy)
    u = exact.composite_utility(factors, ploidy=ploidy)
    return q, log_z, u


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--draws", type=int, default=20000)
    ap.add_argument("--json", type=Path, default=None)
    a = ap.parse_args(argv)

    (ca, cb), bases, shares, betas = build()
    print(f"# The mechanism on a two-chromosome toy genome, exactly\n")
    print(f"- **chrA** K = {ca.K}, T = {ca.T}, {bases[0].M} output paths "
          f"(of {bases[0].n_enumerated:,} state paths enumerated)")
    print(f"- **chrB** K = {cb.K}, T = {cb.T}, {bases[1].M} output paths")
    print(f"- shares by **decision J** (variant count): chrA "
          f"{shares[0]:.6f} = {ca.T}/{ca.T + cb.T}, chrB {shares[1]:.6f}")
    print(f"- mean `beta_v` = {1.0/(ca.T+cb.T):.4f} per position, against "
          f"**3.27e-06** on chr21 — about **{(1.0/(ca.T+cb.T))/3.269192e-06:,.0f}x** larger")
    print(f"- haploid product support {bases[0].M * bases[1].M}, "
          f"diploid {(bases[0].M * bases[1].M)**2:,}\n")

    # the saturating input pair, generalised across chromosomes
    p_hap = (ca.haplotypes[0].copy(), cb.haplotypes[0].copy())
    q_hap = (ca.haplotypes[1].copy(), cb.haplotypes[1].copy())
    r_joint = joint_baseline(bases)

    rows = []
    for tau in TAUS:
        eta = C.eta_tau(tau)
        qp, lzp, up = release_and_u(bases, betas, [p_hap], eta, 1)
        qq, lzq, uq = release_and_u(bases, betas, [q_hap], eta, 1)
        # EXACT moments by enumeration -- not finite differences. The derivative
        # route returns a negative variance at genome-wide shares (chunk 4), so the
        # demonstration must not depend on it even where it would work.
        e_u = float(qp @ up)
        var_u = float(qp @ (up - e_u) ** 2)
        rows.append({
            "tau": tau, "eta_tau": eta,
            "tv_from_prior": exact.tv(qp, r_joint),
            "tv_between_inputs": exact.tv(qp, qq),
            "ratio_spread": exact.ratio_spread(qp, qq),
            "exp_eps": math.exp(C.epsilon_tau(tau)),
            "adv": A.exact_advantage(qp, qq),
            "p_succ_bound": C.p_succ_bound(tau),
            "log_z": lzp, "e_u": e_u, "var_u": var_u,
            "u_max": float(up.max()),
        })

    print("## Haploid, exact. Metrics are rows; `tau` is the column.\n")
    hdr = "| metric | " + " | ".join(f"{r['tau']}" for r in rows) + " |"
    print(hdr); print("|" + "---|" * (len(rows) + 1))

    def line(label, key, fmt="{:.6f}", star=False):
        cells = " | ".join(fmt.format(r[key]) for r in rows)
        print(f"| {'**' if star else ''}{label}{'**' if star else ''} | {cells} |")

    line("TV(Q, R_D) — how far the release moves from the prior", "tv_from_prior", star=True)
    line("`two_target_indistinguishability` *", "tv_between_inputs")
    line("`tau` (the bound on the row above)", "tau", "{:.2f}")
    line("`ratio_spread` *", "ratio_spread", "{:.6f}")
    line("exp(`epsilon_tau`) (its bound)", "exp_eps", "{:.6f}")
    line("optimal attacker accuracy *", "adv")
    line("`p_succ_bound` *", "p_succ_bound")
    line("`log_z_p` *", "log_z")
    line("E_Q[u]", "e_u")
    line("Var_Q(u)", "var_u", "{:.3e}")
    line("max u over the support", "u_max")

    print("\n### The same quantities, DIPLOID (the Arm 1 release shape)\n")
    drows = []
    for tau in TAUS:
        eta = C.eta_tau(tau)
        qp, lzp, up = release_and_u(bases, betas, [p_hap, q_hap], eta, 2)
        qq, lzq, uq = release_and_u(bases, betas, [q_hap, p_hap], eta, 2)
        drows.append({"tau": tau, "log_z": lzp, "eta": eta,
                      "tv_between_inputs": exact.tv(qp, qq),
                      "adv": A.exact_advantage(qp, qq),
                      "p_succ_bound": C.p_succ_bound(tau),
                      "e_u": float(qp @ up), "u_max": float(up.max())})
    print("| metric | " + " | ".join(f"{r['tau']}" for r in drows) + " |")
    print("|" + "---|" * (len(drows) + 1))
    for lbl, k, f in (("`log_z_p` (summed over 4 chains) *", "log_z", "{:.6f}"),
                      ("[0, `eta_tau`] — its bracket", "eta", "{:.6f}"),
                      ("`two_target_indistinguishability` *", "tv_between_inputs", "{:.6f}"),
                      ("optimal attacker accuracy *", "adv", "{:.6f}"),
                      ("`p_succ_bound` *", "p_succ_bound", "{:.6f}"),
                      ("E_Q[u]", "e_u", "{:.6f}")):
        print(f"| {lbl} | " + " | ".join(f.format(r[k]) for r in drows) + " |")

    # the sampler against the enumeration -- the end-to-end correctness check
    print(f"\n### The sampler against the enumeration, {a.draws:,} draws per `tau`\n")
    utils = [utility.beta_v(ca.support, share=shares[0], chrom="chrA"),
             utility.beta_v(cb.support, share=shares[1], chrom="chrB")]
    models = [hmm.from_arrays(ca.haplotypes, ca.genetic_pos, ca.positions, cohort_id="chrA"),
              hmm.from_arrays(cb.haplotypes, cb.genetic_pos, cb.positions, cohort_id="chrB")]
    norms = [ffbs.baseline_log_norm(m) for m in models]
    keys = [{row.tobytes(): i for i, row in enumerate(b.paths)} for b in bases]
    print("| `tau` | max abs(empirical - exact) | 4-sigma band | within band? | "
          "genome `log_z_p` exact | sampled |")
    print("|---|---|---|---|---|---|")
    samp = []
    for tau in TAUS:
        eta = C.eta_tau(tau)
        qp, lzp, _ = release_and_u(bases, betas, [p_hap], eta, 1)
        counts = np.zeros(bases[0].M * bases[1].M, dtype=np.int64)
        lz_obs = []
        for s in range(a.draws):
            idx, parts = [], []
            for c, (m, u, nrm) in enumerate(zip(models, utils, norms)):
                # ⚠ A DISTINCT SEED PER CHROMOSOME. The first version of this loop
                # passed the same `seed` to both, and `release` derives its streams
                # from that seed alone -- so chrA and chrB drew from the SAME stream
                # and their outputs were correlated. The marginals stayed correct and
                # the JOINT did not, which is why the check passed at low `tau` and
                # failed at 0.9 and 0.99 where the distributions carry structure for
                # the correlation to distort. The chromosomes of one release are
                # independent chains and must be seeded as such.
                rel = ffbs.release(m, (ca, cb)[c].haplotypes, p_hap[c], u, tau,
                                   seed=90_000 + 2 * s + c, target_id="toygenome",
                                   graph_id="toy", baseline_norm=nrm)
                idx.append(keys[c][np.ascontiguousarray(rel.released_path).tobytes()])
                parts.append(rel)
            counts[idx[0] * bases[1].M + idx[1]] += 1
            lz_obs.append(math.fsum(p.log_z_p for p in parts))
        emp = counts / a.draws
        sd = np.sqrt(np.maximum(qp, 1e-12) * (1 - qp) / a.draws)
        worst = float(np.abs(emp - qp).max())
        band = float((4.0 * sd + 1e-3).max())
        ok = bool((np.abs(emp - qp) <= 4.0 * sd + 1e-3).all())
        samp.append({"tau": tau, "worst": worst, "within": ok,
                     "lz_exact": lzp, "lz_sampled": float(np.mean(lz_obs))})
        print(f"| {tau} | {worst:.5f} | {band:.5f} | {'**yes**' if ok else '**NO**'} | "
              f"{lzp:.9f} | {np.mean(lz_obs):.9f} |")

    if not all(s["within"] for s in samp):
        raise AssertionError(
            "the sampler departs from the exact release distribution. That is the "
            "whole mechanism being wrong, not a tolerance to widen.")

    out = {"shares": shares, "haploid": rows, "diploid": drows, "sampling": samp,
           "support": {"chrA": bases[0].M, "chrB": bases[1].M},
           "mean_beta": 1.0 / (ca.T + cb.T), "draws": a.draws}
    if a.json:
        a.json.write_text(json.dumps(out, indent=1))
        print(f"\nJSON -> {a.json}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
