#!/usr/bin/env python
"""Exact KL(Q_p || Q_q) with zero Monte-Carlo error, and the bound it implies.

⚠⚠ WHY THIS EXISTS, AND WHY IT SHOULD CARRY THE CLAIM. The sampled attack CANNOT
resolve this quantity at any affordable n: measured 2026-10-11, `sd(Lambda)` per draw is
6.15e-04 against a true value of 1.784e-07, so resolving `E[Lambda]` by sampling needs
about **1.2e7 draws**. The pre-registered 250 trials per cell are 5 orders of magnitude
short. This route computes the same quantity EXACTLY in under a minute per target.

HOW. `Lambda(y) = eta (u_p - u_q) - (log Z_p - log Z_q)` and
`E_{Q_p}[Lambda] = KL(Q_p || Q_q)`. Both `E_{Q_p}[u_p]` and `E_{Q_p}[u_q]` are LINEAR
functionals of the posterior state marginals under `Q_p`, so each is a single pass over
`ffbs.state_marginals` -- no sampling and no differencing. `log Z_p` and `log Z_q` come
from `ffbs.log_z_p`, which is exact.

WHAT IT LICENSES. By Pinsker, `TV <= sqrt(KL/2)`, and the optimal equal-prior attacker's
accuracy is `(1 + TV)/2`. So an exact KL gives an exact UPPER BOUND on what ANY attacker
can achieve -- a strictly stronger statement than "we measured one attacker at chance",
because it covers attackers nobody has written.

✅ DECISION L (Dylan, 2026-10-11): THIS IS ITS OWN NAMED QUANTITY, NOT A SUBSTITUTE FOR
A2. The pre-registration's A2 is "exact enumeration of TV and `ratio_spread` over all
admissible input pairs", and its positive control reads "already passes" -- so A2 was
written as the TOY check, where all 2,304 supported paths can be enumerated, and it
stays satisfied there. On chr21 the support cannot be enumerated at all, so nothing here
replaces A2; `kl_exact` and `acc_ceiling` are new names for new quantities
(`docs/terms.md`).

⚠ WHY INTRODUCING IT AFTER SIGNING IS NOT SELECTIVE REPORTING. It is DETERMINISTIC given
the model and the input pair: there is exactly one value, so there is nothing to choose
among after seeing an outcome. A measured attacker's accuracy would not have that
property, which is why A1 had to be pre-registered and this does not.

⚠⚠ BUT PINSKER IS LOOSE, AND THAT MUST BE SAID EVERY TIME. TV <= sqrt(KL/2) is an
inequality, so the TRUE ceiling is closer to 0.5 than the number printed here. The
looseness runs in the direction that FLATTERS US, which is exactly why it is stated
rather than left implicit.

    python scripts/exact_kl.py --taus 0.25 --targets 2
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

import attacks as A                                       # noqa: E402
from mechanism import calibration as C                    # noqa: E402
from mechanism import ffbs, hmm, utility                  # noqa: E402

CHROM = "chr21"


def load(chrom: str = CHROM) -> dict:
    d = REPO / "data" / "cohort" / chrom
    ic = np.load(d / "in_chain.npy").astype(bool)
    hap = np.load(d / "haplotypes.npy")[:, ic]
    m = hmm.from_arrays(hap, np.load(d / "genetic_pos.npy")[ic],
                        np.load(d / "positions.npy")[ic], cohort_id=chrom)
    util = utility.beta_v(np.load(d / "support.npy"), in_chain=ic)
    paths = {}
    for t in sorted((REPO / "data" / "targets").iterdir()):
        f = t / chrom / "path.npy"
        if f.exists():
            paths[t.name] = np.load(f)[ic]
    return {"hap": hap, "model": m, "util": util, "paths": paths,
            "base": ffbs.baseline_log_norm(m), "ic": ic}


def _hap_terms(st: dict, tgt_hap: np.ndarray, eta_half: float):
    """`log_z` and the posterior marginals for ONE haplotype at `eta/2`."""
    psi = utility.tilt_potentials(st["model"].allowed, st["hap"], tgt_hap,
                                  st["util"].beta, eta_half)
    fwd = ffbs.forward(st["model"], psi)
    lz = ffbs.log_z_p(fwd, st["base"], eta_half,
                      budget=eta_half * float(st["util"].beta.sum()))
    # ⚠ ARGUMENT ORDER IS (hmm, fwd, potentials). I wrote (hmm, psi, fwd) first
    # and the toy validation test is what caught it -- the two arrays have
    # compatible shapes, so a swap does not raise, it returns wrong marginals.
    gamma = ffbs.state_marginals(st["model"], fwd, psi)
    return lz, gamma


def _e_u(st: dict, gamma: np.ndarray, tgt_hap: np.ndarray) -> float:
    """`E[u(tgt, y)]` under the posterior `gamma`, as a linear functional.

    ⚠ `ffbs.expected_mismatches` is UNWEIGHTED and is therefore NOT `E[u]`; the
    utility weights each position by `beta_v`. Verified 2026-10-11: the unweighted
    match mass plus `expected_mismatches` equals the known-cell count to 2.9e-11, which
    is the identity that confirms the marginals are a probability distribution per
    position -- but the weighted sum is what `u` means.
    """
    hap = st["hap"]
    known = tgt_hap >= 0
    match = (hap == tgt_hap[None, :]) & known[None, :]
    return float((gamma * match.T).sum(axis=1) @ st["util"].beta)


def exact_kl(st: dict, p: np.ndarray, q: np.ndarray, tau: float) -> dict:
    """KL(Q_p || Q_q) exactly, plus the Pinsker ceiling on any attacker."""
    eta = C.eta_tau(tau)
    if eta == 0.0:
        # ⚠ tau = 0 is the floor: Q_p == Q_q == R_D, so the KL is exactly 0 and the
        # ceiling is exactly a coin. Returned without computing, because the
        # eta_half = 0 potentials are all-ones and the result is definitional.
        return {"tau": tau, "eta_tau": 0.0, "kl": 0.0, "tv_upper": 0.0,
                "acc_upper": 0.5, "p_succ_bound": 0.5, "log_z_p": 0.0, "log_z_q": 0.0}
    h = eta / 2.0
    lz_p = lz_q = 0.0
    e_up = e_uq = 0.0
    for col in (0, 1):
        lzp, gam_p = _hap_terms(st, p[:, col], h)
        lzq, _ = _hap_terms(st, q[:, col], h)
        lz_p += lzp
        lz_q += lzq
        # ⚠ BOTH expectations are under Q_p's posterior -- that is what makes this
        # E_{Q_p}[.] and therefore the KL. Taking E[u_q] under Q_q's own posterior
        # would compute something else entirely and would not be signed.
        e_up += 0.5 * _e_u(st, gam_p, p[:, col])
        e_uq += 0.5 * _e_u(st, gam_p, q[:, col])
    kl = eta * (e_up - e_uq) - (lz_p - lz_q)
    # ⚠ KL >= 0 is not a convention, it is a theorem. A negative value here means the
    # computation is wrong, not that the distributions are unusual.
    if kl < -1e-12:
        raise AssertionError(
            f"exact KL(Q_p||Q_q) = {kl!r} < 0 at tau={tau}, which is impossible. "
            f"eta*(E u_p - E u_q) = {eta * (e_up - e_uq)!r}, "
            f"log Z_p - log Z_q = {lz_p - lz_q!r}. The two nearly cancel, so suspect "
            f"the normalisers before the expectations.")
    kl = max(kl, 0.0)
    tv_up = min(1.0, math.sqrt(kl / 2.0))          # Pinsker
    return {"tau": tau, "eta_tau": eta, "kl": kl, "tv_upper": tv_up,
            "acc_upper": 0.5 * (1.0 + tv_up), "p_succ_bound": C.p_succ_bound(tau),
            "log_z_p": lz_p, "log_z_q": lz_q, "e_u_p": e_up, "e_u_q": e_uq}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--taus", type=float, nargs="+",
                    default=[0.0, 0.05, 0.10, 0.25, 0.50, 0.90, 0.99])
    ap.add_argument("--targets", type=int, default=10)
    ap.add_argument("--json", type=Path, default=None)
    a = ap.parse_args(argv)

    # ⚠ FIRST, BEFORE ANYTHING IS COMPUTED. Same lock as the sampled attack.
    prereg, rev = A.require_signed_prereg()
    print(f"running under {prereg} (git_rev {rev or 'unrecorded'})")
    print("✅ decision L (2026-10-11): `kl_exact` and `acc_ceiling` are their OWN named "
          "quantities, NOT a substitute for A2. A2 is exact enumeration over all "
          "admissible input pairs — the TOY check, where it already passes — and it is "
          "not computable at chr21 scale.")
    print("⚠ Pinsker is an INEQUALITY: the true ceiling is closer to 0.5 than the "
          "`acc_ceiling` row below. Quote the row WITH this caveat; the looseness runs "
          "in the direction that flatters us.\n")

    st = load()
    tids = sorted(st["paths"])[: a.targets]
    out = []
    t0 = time.time()
    for t in tids:
        p = st["paths"][t]
        others = [k for k in st["paths"] if k != t]
        # the EXTREMAL decoy, built exactly as the sampled attack builds it
        ext = p.copy()
        called = ext >= 0
        ext[called] = np.where(ext[called] == 0, 1, 0)
        for tau in a.taus:
            r = exact_kl(st, p, ext, tau)
            r.update(target=t, decoy="D1_extremal")
            out.append(r)
        print(f"  {t} done ({time.time() - t0:.0f}s)", flush=True)

    print("\nEXACT KL and the ceiling it puts on ANY attacker. D1_extremal, "
          f"{CHROM}, T = {int(st['ic'].sum()):,}, K = {st['hap'].shape[0]}, "
          f"n = {len(tids)} targets, ZERO Monte-Carlo error.\n")
    print("| metric | " + " | ".join(f"{t:g}" for t in a.taus) + " |")
    print("|" + "---|" * (len(a.taus) + 1))
    # ⚠ Row labels are the exact `docs/terms.md` identifiers, gates marked `*`, per
    # the reporting contract. A prose synonym here is what the contract test catches.
    for name, key, fmt in (("`kl_exact` (max over targets)", "kl", "{:.3e}"),
                           ("`tv_pinsker` upper bound", "tv_upper", "{:.3e}"),
                           ("`acc_ceiling` *", "acc_upper", "{:.6f}"),
                           ("`p_succ_bound` *", "p_succ_bound", "{:.4f}")):
        cells = []
        for tau in a.taus:
            vals = [r[key] for r in out if abs(r["tau"] - tau) < 1e-12]
            cells.append(fmt.format(max(vals)))
        print(f"| {name} | " + " | ".join(cells) + " |")
    # ⚠ THE GATE. `acc_ceiling` must not exceed `p_succ_bound`: the first is what this
    # cohort and pair actually permit, the second is what Theorem 1 promises. If the
    # computed ceiling is LOOSER than the theorem, the computation is wrong -- and it
    # is the computation to suspect, not the theorem.
    bad = [(r["target"], r["tau"], r["acc_upper"], r["p_succ_bound"]) for r in out
           if r["acc_upper"] > r["p_succ_bound"] + 1e-12]
    print("\n⚠ The `acc_ceiling` row is what the theorem bounds, computed exactly "
          "rather than sampled, and it is an UPPER bound on EVERY attacker including "
          "unwritten ones — strictly stronger than one measured attacker's accuracy.")
    if bad:
        print(f"\n⚠⚠ {len(bad)} cell(s) where `acc_ceiling` EXCEEDS `p_succ_bound`. "
              f"That is impossible if both are computed correctly, and the "
              f"COMPUTATION is what to suspect:")
        for tgt, tau, a_, b_ in bad[:5]:
            print(f"    {tgt} at tau={tau}: {a_:.6f} > {b_:.4f}")
    print("⚠ Pinsker is not tight; the true TV is smaller, so the real ceiling is "
          "closer to 0.5 than shown.")
    if a.json:
        a.json.write_text(json.dumps(
            {"prereg": prereg, "git_rev": rev, "chrom": CHROM,
             "T": int(st["ic"].sum()), "K": int(st["hap"].shape[0]),
             "method": "exact KL via state marginals; TV bounded by Pinsker",
             "substitutes_for_a2": False,
             "decision": "L (2026-10-11): own named quantity, not an A2 substitute", "results": out}, indent=2) + "\n")
        print(f"\nJSON -> {a.json}")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
