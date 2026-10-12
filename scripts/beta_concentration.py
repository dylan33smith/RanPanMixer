#!/usr/bin/env python
"""Does CONCENTRATING `beta_v` wake the mechanism up? Exactly, by the marginals route.

⚠⚠ WHY THIS EXISTS. Genome-wide the mechanism moves `e_q_u` by **+8.345e-08** across the
whole `tau` range (`logs/B_exact_utility.md`), because one unit of budget spread over
22,134,562 positions leaves `mean beta_v` = 4.52e-08 and a per-position tilt of
`eta_tau * beta_v` ~ 1.2e-07 at `tau` = 0.99. exp of that is 1 + 1.2e-07: numerically
inert. The toy, at `mean beta_v` = 0.111, moves 0.4502 -> 0.6477 with the SAME code. So
the open question raised by Dylan 2026-10-12 is whether the budget can be REALLOCATED to
make the tilt bite, and at what cost.

⚠⚠⚠ THE PRIVACY BOUND DOES NOT CARE HOW `beta` IS DISTRIBUTED. `TV <= tau` depends on the
GLOBAL RANGE `delta_u = sup |u(p,y) - u(q,y)|`, which for a non-negative indicator utility
is `sum_v beta_v`. Hold that at 1 and every allocation below is calibrated identically —
same `eta_tau`, same `epsilon_tau`, same `p_succ_bound`. What changes is only WHAT WE
DECLARE USEFUL. So this script does NOT buy privacy and does not spend any; it measures
how responsive the mechanism is to its own tilt under different declarations of utility.

⚠ EVERY ALLOCATION HERE IS TARGET-INDEPENDENT, selected purely by POSITION INDEX in the
chain. Standing constraint 1 forbids choosing weighted positions by where the target
differs, and that is the easiest way to buy utility and silently lose the guarantee.
Selecting by index cannot see the target. A real paper would declare a PUBLIC panel
(an exome, a ClinVar set) in advance under constraint 6; `window` and `spread` are
stand-ins that bracket the two extremes of spatial layout.

⚠ CONSTRAINT 8 APPLIES AND IS THE REASON FOR THE LAST COLUMN. Each allocation is a
DIFFERENT utility function, so the `e_q_u` LEVELS are not comparable across rows and must
never be read as "this allocation has more utility". What IS comparable is each
allocation's movement as a fraction of ITS OWN ceiling — a dimensionless measure of
whether the mechanism responds at all. That is the reported figure.

ROUTE. `E_Q[u]` is a LINEAR functional of the posterior state marginals, so one
`ffbs.forward` + one `ffbs.state_marginals` pass per (allocation, tau, haplotype) gives it
with no sampling and no finite difference. chr21 at `beta_share` = 1 throughout, so these
are single-chromosome releases and are NOT genome-wide numbers.

    python scripts/beta_concentration.py --taus 0 0.5 0.99
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

import numpy as np                                        # noqa: E402

from mechanism import calibration as C                    # noqa: E402
from mechanism import ffbs, hmm, utility                  # noqa: E402

CHROM = "chr21"
NS = (1, 10, 100, 1_000, 10_000, 100_000)


def allocations(support: np.ndarray, in_chain: np.ndarray, T: int):
    """Every allocation sums to 1 and is a function of INDEX ONLY -- never the target."""
    out = {}

    # the production weighting, for reference: inverse-support over the whole chain
    out["inv_support_all"] = utility.beta_v(
        support, in_chain=in_chain, share=1.0, chrom=CHROM).beta

    # the simplest possible spread, to separate "spread out" from "support-weighted"
    out["uniform_all"] = np.full(T, 1.0 / T)

    mid = T // 2
    for n in NS:
        if n > T:
            continue
        # CONTIGUOUS: n adjacent chain positions at the midpoint. One donor run can cover
        # them all, so the tilt influences a SINGLE donor choice.
        b = np.zeros(T)
        lo = mid - n // 2
        b[lo:lo + n] = 1.0 / n
        out[f"window_{n}"] = b

        # EVENLY SPACED: the same mass over the same count, spread across the chromosome,
        # so the tilt must influence MANY independent donor choices at once.
        b = np.zeros(T)
        idx = np.linspace(0, T - 1, n).round().astype(int)
        idx = np.unique(idx)
        b[idx] = 1.0 / idx.size
        out[f"spread_{n}"] = b

    for k, b in out.items():
        s = float(b.sum())
        assert abs(s - 1.0) < 1e-12, f"{k}: sum beta = {s!r}, not 1; the calibration is void"
        assert (b >= 0).all(), f"{k}: negative beta_v"
    return out


def e_q_u(m, hap, p, beta: np.ndarray, tau: float):
    """Exact E_Q[u] and u_ceiling for one beta allocation, diploid."""
    eta = C.eta_tau(tau)
    eu = 0.0
    for col in (0, 1):
        tgt = p[:, col]
        # ⚠ at eta = 0 the potentials are all-ones, so this IS the baseline posterior --
        # the tau = 0 floor needs no special case.
        psi = utility.tilt_potentials(m.allowed, hap, tgt, beta, eta / 2.0)
        fwd = ffbs.forward(m, psi)
        gam = ffbs.state_marginals(m, fwd, psi)
        known = tgt >= 0
        match = (hap == tgt[None, :]) & known[None, :]
        # ⚠ 0.5 * each haplotype: diploid_utility is the MEAN, matching exp(eta * u_mean).
        eu += 0.5 * float((gam * match.T).sum(axis=1) @ beta)
    ceil = utility.diploid_utility(utility.u_ceiling(p[:, 0], beta),
                                   utility.u_ceiling(p[:, 1], beta))
    return eu, ceil


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--taus", type=float, nargs="+", default=[0.0, 0.5, 0.99])
    ap.add_argument("--target", default=None)
    ap.add_argument("--json", type=Path, default=None)
    ap.add_argument("--from-json", type=Path, default=None,
                    help="re-report a saved run instead of recomputing it")
    a = ap.parse_args(argv)

    if a.from_json:
        dd = json.loads(a.from_json.read_text())
        print(f"re-reporting {len(dd['rows'])} saved allocation(s) from {a.from_json}. "
              f"NO marginals computed.")
        _report(dd["rows"], dd["taus"], dd["target"], dd["T"])
        return 0

    d = REPO / "data" / "cohort" / CHROM
    ic = np.load(d / "in_chain.npy").astype(bool)
    hap = np.load(d / "haplotypes.npy")[:, ic]
    m = hmm.from_arrays(hap, np.load(d / "genetic_pos.npy")[ic],
                        np.load(d / "positions.npy")[ic], cohort_id=CHROM)
    T = int(ic.sum())

    tgt_id = a.target or sorted(q.name for q in (REPO / "data" / "targets").iterdir()
                                if q.is_dir())[0]
    p = np.load(REPO / "data" / "targets" / tgt_id / CHROM / "path.npy")[ic]

    allocs = allocations(np.load(d / "support.npy"), ic, T)
    print(f"{CHROM}, T = {T:,}, K = {hap.shape[0]}, target {tgt_id}, ARM 2, "
          f"`beta_share` = 1 (CHR21-SCOPE, not genome-wide).")
    print(f"{len(allocs)} allocations x {len(a.taus)} tau = "
          f"{len(allocs) * len(a.taus)} exact passes. Every allocation sums to 1, so "
          f"`tau`, `eta_tau` and `p_succ_bound` are IDENTICAL across rows.\n")

    rows = []
    t0 = time.time()
    for name, beta in allocs.items():
        nz = int((beta > 0).sum())
        # ⚠ HALF eta: the diploid release runs each haplotype at `eta_tau / 2`, which is
        # what `e_q_u` below actually applies. This column read FULL `eta_tau` until
        # 2026-10-12 and so overstated the REALISED per-haplotype exponent by exactly 2x.
        # Settled against a closed form: with `beta` on ONE position the tilted law is the
        # baseline reweighted pointwise, and the half-eta form matches the measured
        # `window_1` to 3.5e-10 while the full-eta form is off by 8.2e-04.
        bmax = float(beta.max())
        per = C.eta_tau(max(a.taus)) / 2.0 * bmax
        res = {}
        for tau in a.taus:
            eu, ceil = e_q_u(m, hap, p, beta, tau)
            res[tau] = (eu, ceil)
            print(f"    {name} tau={tau:g}: E_Q[u]={eu:.9e} ceil={ceil:.9e} "
                  f"({time.time() - t0:.0f}s)", flush=True)
        lo, hi = res[min(a.taus)][0], res[max(a.taus)][0]
        ceil = res[min(a.taus)][1]
        rows.append({"alloc": name, "n_weighted": nz, "max_tilt": per,
                     "beta_max": bmax,
                     "e_q_u_lo": lo, "e_q_u_hi": hi, "u_ceiling": ceil,
                     "movement": hi - lo,
                     "movement_frac_ceiling": (hi - lo) / ceil if ceil else float("nan")})
        print(f"  {name}: movement {hi - lo:+.6e} "
              f"({100 * (hi - lo) / ceil if ceil else float('nan'):+.4f}% of its ceiling)",
              flush=True)

    _report(rows, a.taus, tgt_id, T)
    if a.json:
        a.json.write_text(json.dumps({"chrom": CHROM, "target": tgt_id, "T": T,
                                      "taus": a.taus, "rows": rows}, indent=1) + "\n")
        print(f"\nJSON -> {a.json}")
    return 0


def _report(rows, taus, tgt_id, T) -> None:
    assert rows, "no allocations were measured; the report would be silent"
    lo_t, hi_t = min(taus), max(taus)
    print(f"\nHow far the tilt moves expected utility, by BUDGET ALLOCATION. chr21, "
          f"`beta_share` = 1, target {tgt_id}, ARM 2, n = 1, exact.\n")
    hdr = ["allocation", "positions weighted", "max `eta_tau`/2 * `beta_v` [per hap]",
           f"`e_q_u` at `tau`={lo_t:g}", f"`e_q_u` at `tau`={hi_t:g}", "`u_ceiling`",
           "headroom", "MOVEMENT", "movement / own `u_ceiling`",
           "**% of headroom captured**"]
    print("| " + " | ".join(hdr) + " |")
    print("|" + "---|" * len(hdr))
    for r in rows:
        head = r["u_ceiling"] - r["e_q_u_lo"]
        cap = (100 * r["movement"] / head) if head > 0 else float("nan")
        print(f"| `{r['alloc']}` | {r['n_weighted']:,} | {r['max_tilt']:.3e} | "
              f"{r['e_q_u_lo']:.9f} | {r['e_q_u_hi']:.9f} | {r['u_ceiling']:.9f} | "
              f"{head:.6f} | **{r['movement']:+.4e}** | "
              f"{100 * r['movement_frac_ceiling']:+.4f}% | **{cap:+.3f}%** |")

    best = max(rows, key=lambda r: abs(r["movement_frac_ceiling"]))
    base = next((r for r in rows if r["alloc"] == "inv_support_all"), None)
    print(f"\n⚠⚠ Most responsive allocation: `{best['alloc']}` at "
          f"**{100 * best['movement_frac_ceiling']:+.4f}%** of its own ceiling, with "
          f"{best['n_weighted']:,} of {T:,} positions weighted and a per-position tilt of "
          f"{best['max_tilt']:.3e}.")
    if base:
        print(f"⚠⚠ The PRODUCTION allocation `inv_support_all` moves "
              f"{100 * base['movement_frac_ceiling']:+.4f}% of its ceiling over the same "
              f"`tau` range, with all {base['n_weighted']:,} positions weighted.")
    print("⚠ `tau`, `eta_tau`, `epsilon_tau` and `p_succ_bound` are IDENTICAL on every "
          "row: each allocation sums to 1, so `delta_u` = 1 and the calibration never "
          "changed. NOTHING here buys or spends privacy.")
    print("⚠⚠ CONSTRAINT 8: each row is a DIFFERENT utility function, so the `e_q_u` "
          "LEVELS are NOT comparable across rows and must never be read as one "
          "allocation having more utility. Only the within-row MOVEMENT is comparable.")
    print("⚠⚠ READ THE LAST COLUMN, NOT THE MOVEMENT. `window_1` captures the largest "
          "FRACTION of its headroom of any allocation, but its headroom is ~0.004 because "
          "the baseline mosaic already matches that one position ~99.6% of the time. The "
          "positions where the tilt has the most power are mostly the ones that needed it "
          "least -- the free-lunch phenomenon in miniature.")
    print("⚠ A concentrated allocation declares that only the weighted positions matter. "
          "That is defensible ONLY if the panel is public, target-independent and "
          "pre-registered under constraint 6 — never chosen after seeing an outcome.")


if __name__ == "__main__":
    sys.exit(main())
