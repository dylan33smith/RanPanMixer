#!/usr/bin/env python
"""Arm 2's attack experiment: the optimal distinguisher, run under a signed plan.

⚠ IT REFUSES TO RUN WITHOUT A SIGNED PRE-REGISTRATION, and that is the first thing it
does. Standing constraint 6 fixes the attacker set, the `tau` grid, the endpoints and
the attack database before any outcome is inspected, because the failure it guards
against is one nobody can self-police: having seen a result, you cannot PROVE you did
not choose the analysis to suit it.

WHAT THE ATTACK IS. A guessing game. The adversary is handed one released path and
told it came from either input `p` or input `q`, holds both genomes, and must say
which. Guessing blind wins half the time. The theorem says they cannot beat
`(1 + tau) / 2`.

WHY THIS ATTACKER AND NOT ONE FROM THE LITERATURE. The proposal specifies no attacker,
but its algebra gives the provably optimal one. With
`Q^tau_p(y) = R_D(y) exp(eta_tau u(p,y)) / Z_p`, the Neyman-Pearson statistic is

    Lambda(y) = eta_tau [u(p,y) - u(q,y)] - [log Z_p - log Z_q]

with `R_D` cancelling, and by Neyman-Pearson its advantage EQUALS `TV(Q_p, Q_q)` — the
quantity the theorem bounds. No attack in any literature can beat it for a given pair.
That matters for reporting a null: in the DP-auditing literature the trade-off function
is unknown and must be approximated, so "was the attacker strong enough" is open. Here
it is closed form, so this attacker provably attains the optimum — but only if the
PRE-REGISTERED attacker is this exact likelihood ratio rather than a heuristic.

⚠ SO STRENGTH IS NOT THE RISK; POWER IS. The repository forecasts a chance-level null
at the operating configuration. A null without a power statement is not evidence
(standing constraint 7). The pre-registered design is 10 targets x 50 draws per `tau`,
giving a minimum detectable advantage of 0.056 — which licenses a claim at `tau` >= 0.25
and NOT at 0.05 or 0.10, where the bound is already tighter than any affordable
experiment. Those rows are still reported, labelled as not resolving below the bound.

    python scripts/arm2_attack.py [--draws 50] [--targets 10] [--json out.json]
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
import runner as R  # noqa: E402
from mechanism import calibration as C  # noqa: E402
from mechanism import ffbs, hmm, utility  # noqa: E402

TAUS = (0.05, 0.10, 0.25, 0.50, 0.90, 0.99)
CHROM = "chr21"
DECOY_SEED = 20261010

_STATE: dict = {}


def load_state(chrom: str = CHROM) -> dict:
    """Cohort arrays, the model, and every prepared target's path, once."""
    if _STATE:
        return _STATE
    d = REPO / "data" / "cohort" / chrom
    hap = np.load(d / "haplotypes.npy")
    ic = np.load(d / "in_chain.npy").astype(bool)
    sup = np.load(d / "support.npy")
    gp = np.load(d / "genetic_pos.npy")
    pos = np.load(d / "positions.npy")
    m = hmm.from_arrays(hap[:, ic], gp[ic], pos[ic], cohort_id=chrom)
    util = utility.beta_v(sup, in_chain=ic)       # scope 'chain', share 1 on one chrom

    paths = {}
    for t in sorted((REPO / "data" / "targets").iterdir()):
        f = t / chrom / "path.npy"
        if f.exists():
            paths[t.name] = np.load(f)[ic]
    if len(paths) < 2:
        raise RuntimeError(
            f"found {len(paths)} prepared target(s) under data/targets/*/{chrom}. The "
            f"experiment needs a target AND a decoy, so at least two.")

    _STATE.update(hap=hap[:, ic], util=util, model=m, paths=paths,
                  base=ffbs.baseline_log_norm(m))
    return _STATE


def concordance(a: np.ndarray, b: np.ndarray) -> float:
    """Share of positions where two diploid paths carry the same allele pair.

    Positions either one does not call are excluded, so a target with more missing
    data is not scored as more different from everyone.
    """
    known = (a >= 0).all(axis=1) & (b >= 0).all(axis=1)
    if not known.any():
        return 0.0
    return float((a[known] == b[known]).all(axis=1).mean())


def pick_decoys(target_id: str, st: dict, rng: np.random.Generator) -> dict:
    """The three pre-registered decoy conditions for one target.

    D1 extremal — constructed, differing from the target at every called position,
      so `u_p - u_q` spans its full range and the ratio bound is saturated.
    D2 nearest — the highest-concordance other target. THE HARDEST REAL DECOY.
    D3 random — one other target drawn uniformly, the easy baseline.

    ⚠ A LIMITATION TO STATE, NOT TO HIDE. D2's pool is the OTHER PREPARED TARGETS —
    nine of them — not the 3,202-sample panel. Nine candidates is a weak version of
    "the hardest real decoy"; the panel version needs a site-aligned projection of
    the panel onto our chain axis, which does not exist yet. So D2 as run here is a
    LOWER BOUND on decoy difficulty, and the real nearest neighbour can only be
    harder. Widening the pool strengthens the attacker, so this limitation flatters
    us and must be reported with the result.
    """
    p = st["paths"][target_id]
    others = [k for k in st["paths"] if k != target_id]

    # D1: flip every called allele to something the target does not carry. The
    # cohort's allele indices start at 0, so 1 - allele is wrong for multi-allelic
    # sites; use "any donor allele that differs", which is what the support admits.
    extremal = p.copy()
    called = extremal >= 0
    extremal[called] = np.where(extremal[called] == 0, 1, 0)

    conc = {k: concordance(p, st["paths"][k]) for k in others}
    nearest = max(conc, key=conc.get)
    random_id = others[int(rng.integers(len(others)))]
    return {
        "D1_extremal": {"path": extremal, "id": "constructed", "concordance": None},
        "D2_nearest": {"path": st["paths"][nearest], "id": nearest,
                       "concordance": conc[nearest]},
        "D3_random": {"path": st["paths"][random_id], "id": random_id,
                      "concordance": conc[random_id]},
    }


def _one_trial(job: tuple) -> dict:
    """One release, scored against all three decoys.

    ⚠ ONE RELEASE SERVES ALL THREE CONDITIONS, which is correct and not a shortcut:
    a release is drawn from `Q_p` and does not depend on `q`. But it makes the three
    conditions CORRELATED through the shared draw, so they must never be pooled into
    a single accuracy.
    """
    target_id, tau, seed = job
    st = load_state()
    eta = C.eta_tau(tau)
    p = st["paths"][target_id]
    rng = np.random.default_rng([DECOY_SEED, hash(target_id) % (1 << 31)])
    decoys = pick_decoys(target_id, st, rng)

    rel = ffbs.release(st["model"], st["hap"], p, st["util"], tau, seed=seed,
                       target_id=target_id, graph_id=CHROM,
                       baseline_norm=st["base"])
    y = rel.released_path
    u_p = utility.diploid_utility(
        *[float(utility.u_of_output(y[:, h], p[:, h], st["util"].beta)[0])
          for h in (0, 1)])

    out = {"target": target_id, "tau": tau, "seed": seed,
           "log_z_p": rel.log_z_p, "u_p": u_p, "conditions": {}}
    for name, d in decoys.items():
        q = d["path"]
        u_q = utility.diploid_utility(
            *[float(utility.u_of_output(y[:, h], q[:, h], st["util"].beta)[0])
              for h in (0, 1)])
        # log Z_q for the SAME release slot: a second forward pass, no sampling
        psi_q = utility.tilt_potentials(st["model"].allowed, st["hap"], q[:, 0],
                                        st["util"].beta, eta / 2.0)
        lz_q0 = ffbs.log_z_p(ffbs.forward(st["model"], psi_q), st["base"], eta / 2.0,
                             budget=(eta / 2.0) * float(st["util"].beta.sum()))
        psi_q1 = utility.tilt_potentials(st["model"].allowed, st["hap"], q[:, 1],
                                         st["util"].beta, eta / 2.0)
        lz_q1 = ffbs.log_z_p(ffbs.forward(st["model"], psi_q1), st["base"], eta / 2.0,
                             budget=(eta / 2.0) * float(st["util"].beta.sum()))
        lam = float(A.lrt_statistic(np.array([u_p]), np.array([u_q]), eta,
                                    rel.log_z_p, lz_q0 + lz_q1)[0])
        out["conditions"][name] = {"decoy": d["id"], "u_q": u_q, "lambda": lam,
                                   "concordance": d["concordance"]}
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--draws", type=int, default=50)
    ap.add_argument("--targets", type=int, default=10)
    ap.add_argument("--json", type=Path, default=None)
    ap.add_argument("--taus", type=float, nargs="+", default=list(TAUS))
    a = ap.parse_args(argv)

    # ⚠ FIRST, BEFORE ANYTHING IS COMPUTED.
    prereg, rev = A.require_signed_prereg()
    print(f"running under {prereg} (git_rev {rev or 'unrecorded'})\n")

    st = load_state()
    tids = sorted(st["paths"])[: a.targets]
    jobs = [(t, tau, 50_000 + i) for tau in a.taus for t in tids
            for i in range(a.draws)]
    n_per = len(tids) * a.draws
    delta = 1.243 / math.sqrt(n_per)

    plan = R.plan_workers(R.array_bytes(st["model"].T, st["model"].K, n_arrays=4),
                          n_jobs=len(jobs))
    print(f"{len(tids)} targets x {a.draws} draws x {len(a.taus)} tau = "
          f"{len(jobs):,} releases")
    print(f"{n_per} trials per (tau, condition) -> minimum detectable advantage "
          f"delta = {delta:.4f}")
    print(f"{plan}\n")

    res = R.run_jobs(_one_trial, jobs, workers=plan.workers, label="trials",
                     progress_every=max(1, len(jobs) // 10))

    print("Accuracy of the optimal distinguisher. Metrics are rows; `tau` is the "
          "column.\n")
    hdr = "| metric | " + " | ".join(str(t) for t in a.taus) + " |"
    print(hdr); print("|" + "---|" * (len(a.taus) + 1))
    table: dict = {}
    for cond in ("D1_extremal", "D2_nearest", "D3_random"):
        cells = []
        for tau in a.taus:
            lams = [r["conditions"][cond]["lambda"] for r in res if r["tau"] == tau]
            # the sign rule, ties split: the release came from p, so Lambda > 0 is right
            acc = (sum(1 for l in lams if l > 0) + 0.5 * sum(1 for l in lams if l == 0)) / len(lams)
            table.setdefault(cond, {})[tau] = acc
            cells.append(f"{acc:.4f}")
        print(f"| {cond} accuracy * | " + " | ".join(cells) + " |")
    print("| `p_succ_bound` * | " +
          " | ".join(f"{C.p_succ_bound(t):.4f}" for t in a.taus) + " |")
    print("| experiment resolves below the bound? | " +
          " | ".join("**yes**" if delta < t / 2 else "no" for t in a.taus) + " |")

    # ⚠ a measured accuracy ABOVE the bound falsifies the implementation
    violations = [(c, t, v) for c, row in table.items() for t, v in row.items()
                  if v > C.p_succ_bound(t) + 1e-9]
    print()
    if violations:
        print("⚠⚠ MEASURED ACCURACY ABOVE p_succ_bound — this falsifies the "
              "IMPLEMENTATION, not the theorem:")
        for c, t, v in violations:
            print(f"    {c} at tau={t}: {v:.4f} > {C.p_succ_bound(t):.4f}")
    else:
        print("no condition exceeded p_succ_bound at any tau")
    print(f"\n⚠ tau = 0.05 and 0.10 do NOT resolve below the bound at delta = "
          f"{delta:.4f}; those rows are reported, not claimed. The PROOF carries them.")
    print("⚠ D2's decoy pool is the other prepared targets, not the 3,202-sample "
          "panel, so it is a LOWER BOUND on decoy difficulty and flatters us.")
    print("⚠ The three conditions share each release and must not be pooled.")

    if a.json:
        a.json.write_text(json.dumps(
            {"prereg": prereg, "git_rev": rev, "draws": a.draws,
             "targets": tids, "taus": a.taus, "trials_per_cell": n_per,
             "delta": delta, "accuracy": {c: {str(k): v for k, v in r.items()}
                                          for c, r in table.items()},
             "trials": res}, indent=1))
        print(f"\nJSON -> {a.json}")
    return 1 if violations else 0


if __name__ == "__main__":
    sys.exit(main())
