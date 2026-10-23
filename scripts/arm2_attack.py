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
import hashlib
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

# ⚠ `tau` = 0 IS MANDATORY and was missing here. The pre-registration calls it
# "the floor, and mandatory" (section 1) and A3's positive control requires it
# explicitly -- the attacker "must sit at chance against a `tau` = 0 release".
# Without it the run has no control, and a chance-level result at every other
# point cannot be distinguished from an attacker that never worked at all.
TAUS = (0.0, 0.05, 0.10, 0.25, 0.50, 0.90, 0.99)
# ⚠ EXPLORATORY, NOT PRE-REGISTERED (requested 2026-10-22, after the pinned grid's
# outcomes were seen). Listed only so each point keeps a FIXED seed index whatever
# subset of `--taus` is run: the index used to be the position in `--taus`, so running
# these three alone would have reused the seeds of 0, 0.05 and 0.10.
EXPLORATORY_TAUS = (0.999, 0.9999, 0.99999)
SEED_GRID = TAUS + EXPLORATORY_TAUS
CHROM = "chr21"
DECOY_SEED = 20261010
# the release-stream base; the per-cell seed is a SeedSequence over
# (SEED_BASE, target index, tau index, draw index) -- see main().
SEED_BASE = 50_000

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
    target_id, tau, seed, one_sided = job
    st = load_state()
    eta = C.eta_tau(tau)
    p = st["paths"][target_id]
    # ⚠⚠ NOT `hash(target_id)`. Python's str hash is randomised per interpreter by
    # PYTHONHASHSEED, so the decoy choice was NOT REPRODUCIBLE -- measured,
    # hash('HG00096') % 2**31 returned 1730677557, 783898805 and 1043054338 on three
    # invocations of the same binary. The recorded run is self-consistent only because
    # `fork` inherits the parent's hash seed; a re-run picks a DIFFERENT D3 decoy, so
    # that cell could never have been reproduced from its recorded seed. A
    # pre-registered experiment whose decoy set is not reproducible is not
    # pre-registered in the sense constraint 6 means.
    # sha256 is stable across interpreters, machines and versions.
    stable = int.from_bytes(hashlib.sha256(target_id.encode()).digest()[:4], "big")
    rng = np.random.default_rng([DECOY_SEED, stable])
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
        lz_q = lz_q0 + lz_q1
        lam = float(A.lrt_statistic(np.array([u_p]), np.array([u_q]), eta,
                                    rel.log_z_p, lz_q)[0])
        cell = {"decoy": d["id"], "u_q": u_q, "lambda": lam,
                "concordance": d["concordance"]}

        # ⚠⚠ THE SECOND HYPOTHESIS, AND WITHOUT IT THERE IS NO BOUND TO COMPARE TO.
        # Theorem 1's ceiling (1 + tau)/2 bounds the EQUAL-PRIOR accuracy
        #     0.5 * P_p(Lambda > 0) + 0.5 * P_q(Lambda < 0),
        # which needs draws from BOTH input paths. The release above comes from `Q_p`
        # only, so on its own it yields the ONE-SIDED rate P_p(Lambda > 0) -- and that
        # quantity is NOT bounded by (1 + tau)/2. From
        #     balanced = 0.5*TPR + 0.5*(1 - FPR) <= (1 + TV)/2  =>  TPR <= TV + FPR,
        # a theorem-respecting attacker with FPR 0.5 reaches TPR = tau + 0.5, which
        # EXCEEDS (1 + tau)/2 at every tau > 0; and a degenerate attacker that always
        # answers "p" scores 1.0 one-sided while being exactly at chance balanced.
        # Comparing the one-sided rate to the bound would therefore report a privacy
        # violation with the theorem fully intact.
        if not one_sided:
            # an INDEPENDENT seed per (condition, hypothesis): reusing the Q_p seed
            # would correlate the two arms of the test and bias the balance.
            q_seed = int(seed) ^ (0x5EED * (1 + sorted(decoys).index(name)))
            rel_q = ffbs.release(st["model"], st["hap"], q, st["util"], tau,
                                 seed=q_seed, target_id=f"{target_id}~{name}",
                                 graph_id=CHROM, baseline_norm=st["base"])
            yq = rel_q.released_path
            u_p_q = utility.diploid_utility(
                *[float(utility.u_of_output(yq[:, h], p[:, h], st["util"].beta)[0])
                  for h in (0, 1)])
            u_q_q = utility.diploid_utility(
                *[float(utility.u_of_output(yq[:, h], q[:, h], st["util"].beta)[0])
                  for h in (0, 1)])
            # ⚠ The SAME statistic, with the same normalisers -- Lambda is a function
            # of y, and log Z_p / log Z_q do not depend on which hypothesis drew y.
            cell["lambda_under_q"] = float(A.lrt_statistic(
                np.array([u_p_q]), np.array([u_q_q]), eta, rel.log_z_p, lz_q)[0])
        out["conditions"][name] = cell
    return out


def _z_for_one_sided(alpha: float) -> float:
    """`z` with `P(Z > z) = alpha`, for a standard normal. No scipy dependency.

    ⚠ Written out rather than imported because `environment.yml` does not promise
    scipy, and a bound check that cannot run for want of a dependency is a bound check
    that is not running. Newton on `erf`, seeded by Winitzki; the test round-trips it
    against `math.erf` so the implementation is checked rather than trusted.
    """
    if not 0.0 < alpha < 0.5:
        raise ValueError(f"alpha must be in (0, 0.5), got {alpha}")
    y = 1.0 - 2.0 * alpha                       # erf(z/sqrt2) = 1 - 2*alpha
    a = 0.147
    ln = math.log(1.0 - y * y)
    t = 2.0 / (math.pi * a) + ln / 2.0
    x = math.sqrt(math.sqrt(t * t - ln / a) - t)
    for _ in range(4):
        x -= (math.erf(x) - y) / (2.0 / math.sqrt(math.pi) * math.exp(-x * x))
    return x * math.sqrt(2.0)


def _report(res: list, a, delta: float, prereg: str, rev: str,
            tids: list, n_per: int) -> int:
    """Build the table from trial records. SHARED by a live run and a re-report.

    ⚠ IT IS A FUNCTION SO THE REPORTING CAN BE FIXED WITHOUT RE-RUNNING. The releases
    cost ~80 minutes and ~110 CPU-hours; if changing the table meant redoing them, the
    table would never get changed -- which is how the defective bound check below
    survived its first run.
    """
    table: dict = {}
    print("Accuracy of the optimal distinguisher. Metrics are rows; `tau` is the "
          "column.\n")
    hdr = "| metric | " + " | ".join(str(t) for t in a.taus) + " |"
    print(hdr); print("|" + "---|" * (len(a.taus) + 1))
    for cond in ("D1_extremal", "D2_nearest", "D3_random"):
        cells = []
        for tau in a.taus:
            rows = [r["conditions"][cond] for r in res if r["tau"] == tau]
            lp = [c["lambda"] for c in rows]
            # ⚠ TIES SPLIT, which is the right convention at equal priors and is not
            # cosmetic: at tau = 0 the tilt vanishes, Lambda is identically 0, and
            # splitting is what makes the control return exactly 0.5 instead of a
            # number decided by a tie-breaking accident.
            tpr = (sum(1 for l in lp if l > 0) + 0.5 * sum(1 for l in lp if l == 0)) / len(lp)
            if a.one_sided:
                table.setdefault(cond, {})[tau] = {"tpr": tpr, "balanced": None}
                cells.append(f"{tpr:.4f}[1-sided]")
                continue
            lq = [c["lambda_under_q"] for c in rows]
            tnr = (sum(1 for l in lq if l < 0) + 0.5 * sum(1 for l in lq if l == 0)) / len(lq)
            # ⚠ THE QUANTITY THEOREM 1 BOUNDS, and it comes from `attacks` rather than
            # being recomputed here. The inline version this replaced was the bug: a
            # correct balanced `empirical_accuracy` already existed in src/attacks.py
            # and this script reimplemented the accuracy one-sided instead of calling
            # it. Duplicated statistics drift; the helper is the single definition.
            bal = A.empirical_accuracy(np.array(lp), np.array(lq))
            assert abs(bal - (0.5 * tpr + 0.5 * tnr)) < 1e-12, (
                "empirical_accuracy disagrees with 0.5*TPR + 0.5*TNR; one of the two "
                "definitions has drifted and no number here is quotable")
            # ⚠ THE SE TRAVELS WITH THE ESTIMATE. A bound comparison that does not
            # carry its uncertainty is the defect fixed below.
            se_bal = 0.5 * math.sqrt(tpr * (1 - tpr) / len(lp)
                                     + tnr * (1 - tnr) / len(lq))
            table.setdefault(cond, {})[tau] = {"tpr": tpr, "tnr": tnr, "balanced": bal,
                                               "se_balanced": se_bal,
                                               "n_per_arm": len(lp)}
            cells.append(f"{bal:.4f}")
        label = "one-sided P_p(L>0)" if a.one_sided else "balanced accuracy"
        print(f"| {cond} {label} * | " + " | ".join(cells) + " |")
    if not a.one_sided:
        for cond in ("D1_extremal", "D2_nearest", "D3_random"):
            print(f"| {cond} one-sided P_p(L>0) [diagnostic, NOT bounded] | " +
                  " | ".join(f"{table[cond][t]['tpr']:.4f}" for t in a.taus) + " |")
    print("| `p_succ_bound` * | " +
          " | ".join(f"{C.p_succ_bound(t):.4f}" for t in a.taus) + " |")
    print("| experiment resolves below the bound? | " +
          " | ".join("**yes**" if delta < t / 2 else "no" for t in a.taus) + " |")

    # ⚠⚠ THE VIOLATION CHECK APPLIES ONLY TO THE BALANCED ACCURACY. (1 + tau)/2 bounds
    # the equal-prior quantity and nothing else; comparing the one-sided rate to it
    # fires with the theorem intact (see the derivation at the Q_q draw above), so in
    # one-sided mode there is NOTHING here to check and saying so beats a false alarm.
    if a.one_sided:
        print("\n⚠ NO BOUND CHECK PERFORMED. The one-sided rate is not bounded by "
              "p_succ_bound, so a comparison would be meaningless in BOTH directions: "
              "exceeding it would not be a violation, and falling below it would not "
              "be evidence. Re-run without --one-sided for a claimable number.")
        violations = []
    else:
        # ⚠⚠ A POINT COMPARISON AGAINST A NOISY ESTIMATE CRIES WOLF, and this one did.
        # At `tau` = 0.05 the bound sits only 0.025 above chance while the balanced
        # accuracy's SE is 0.0223 -- about 1.1 SE -- so if the true accuracy is 0.5 the
        # estimate exceeds the bound about 13% of the time PER CELL, which is ~2.3
        # expected false alarms over 18 cells. MEASURED 2026-10-11 on the re-run:
        # D3_random at `tau` = 0.05 came in at 0.5400 against 0.5250, an excess of
        # +0.015 = +0.67 SE, and the old check called that a falsification of the
        # implementation.
        # ⚠ It was self-contradictory too: the script already prints that `tau` = 0.05
        # and 0.10 do NOT resolve below the bound, and then claimed a falsification
        # there anyway.
        # The test is now one-sided on the excess, Bonferroni-corrected across the
        # cells actually examined. A REAL implementation error still fires, because a
        # true excess grows with n while noise shrinks as 1/sqrt(n); a point excess
        # that is not significant is reported as a NOTE and does not fail the run.
        cells = [(c, t, v) for c, row in table.items() for t, v in row.items()
                 if v.get("balanced") is not None]
        z_crit = _z_for_one_sided(0.05 / max(1, len(cells)))
        violations, notes = [], []
        for c, t, v in cells:
            bal, bound = v["balanced"], C.p_succ_bound(t)
            if bal <= bound + 1e-9:
                continue
            se = v.get("se_balanced") or 0.0
            z = (bal - bound) / se if se > 0 else float("inf")
            (violations if z > z_crit else notes).append((c, t, bal, bound, z))
        if notes:
            print(f"\n⚠ {len(notes)} cell(s) exceed `p_succ_bound` by POINT ESTIMATE "
                  f"but NOT significantly (one-sided 5%, Bonferroni over "
                  f"{len(cells)} cells, z_crit = {z_crit:.2f}). Reported, not claimed:")
            for c, t, bal, bound, z in notes:
                print(f"    {c} at tau={t}: {bal:.4f} vs bound {bound:.4f} "
                      f"(+{bal - bound:.4f}, z = +{z:.2f})")
        violations = [(c, t, bal) for c, t, bal, _, _ in violations]
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


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    # ⚠ 25 IS THE PRE-REGISTERED VALUE (field 6, signed 2026-10-10): 25 draws per
    # `tau` PER HYPOTHESIS, balanced. A different default here would silently run
    # an unregistered design, and the whole point of the lock above is that the
    # design cannot drift after signing.
    ap.add_argument("--draws", type=int, default=25)
    ap.add_argument("--targets", type=int, default=10)
    ap.add_argument("--json", type=Path, default=None)
    ap.add_argument("--from-json", type=Path, default=None,
                    help="re-report from a saved trials JSON; computes "
                         "NO releases. ⚠ The releases cost ~110 CPU-hours, "
                         "so the reporting must be fixable without them.")
    ap.add_argument("--taus", type=float, nargs="+", default=list(TAUS))
    ap.add_argument("--one-sided", action="store_true",
                    help="draw only from Q_p and report the one-sided rate "
                         "P_p(Lambda>0). ⚠ That quantity is NOT bounded by "
                         "(1+tau)/2, so this mode REFUSES to compare it to the "
                         "bound and cannot support a privacy claim. It costs a "
                         "quarter of the compute and is a tilt diagnostic only.")
    a = ap.parse_args(argv)

    # ⚠ FIRST, BEFORE ANYTHING IS COMPUTED.
    prereg, rev = A.require_signed_prereg()
    print(f"running under {prereg} (git_rev {rev or 'unrecorded'})\n")

    if a.from_json:
        saved = json.loads(a.from_json.read_text())
        res = saved["trials"]
        a.taus = sorted({r["tau"] for r in res})
        tids = sorted({r["target"] for r in res})
        n_per = len(res) // max(1, len(a.taus))
        delta = (1.243 if a.one_sided else 0.879) / math.sqrt(max(1, n_per))
        print(f"re-reporting {len(res):,} saved trial(s) from {a.from_json.name}: "
              f"{len(tids)} targets x {len(a.taus)} tau, {n_per} per (tau, condition) "
              f"per arm. NO releases computed.")
        print(f"minimum detectable advantage delta = {delta:.4f}\n")
        return _report(res, a, delta, prereg, rev, tids, n_per)

    st = load_state()
    tids = sorted(st["paths"])[: a.targets]
    # ⚠⚠ THE SEED MUST DEPEND ON THE WHOLE CELL. The previous form was
    #     (t, tau, 50_000 + i, ...)
    # which derives the seed from the DRAW INDEX ALONE, so all 10 targets and all 7
    # `tau` reused seeds 50000..50024. With a per-position tilt of
    # eta * mean(beta) = 0.2554 / 332,521 = 7.7e-07 the tilt cannot move the draw, so
    # the released path was BIT-IDENTICAL across tau at a fixed seed -- MEASURED:
    # u_p = 0.588012460 at tau = 0, 0.05 AND 0.25, to nine digits. The seven tau
    # columns were ONE experiment reported seven times, the ten targets were
    # near-replicates (cross-target correlation +0.305 at fixed seed), and the naive
    # standard error was about 2x too narrow.
    # ⚠⚠ THIS IS THE THIRD OCCURRENCE OF A TWICE-DOCUMENTED FAILURE. `docs/bugs.md`
    # has "A headline finding measured under common random numbers" (2026-09-30), and
    # `docs/memory.md` records seed reuse across chromosomes correlating two
    # chromosomes of one release. Read those before touching this line.
    # The fix: a SeedSequence over the full key, so every (target, tau, draw) cell
    # draws from a provably disjoint stream while staying reproducible from the key.
    def _cell_seed(ti: int, taui: int, i: int) -> int:
        return int(np.random.SeedSequence([SEED_BASE, ti, taui, i]).generate_state(1)[0])

    # the `tau` index is its place in SEED_GRID, NOT in `--taus`, so a subset run
    # draws the same streams as the full grid; the pinned points keep their indices
    off = [t for t in a.taus if t not in SEED_GRID]
    if off:
        raise SystemExit(f"tau {off} is in neither the pinned grid nor "
                         f"EXPLORATORY_TAUS, so it has no fixed seed index")
    if any(t in EXPLORATORY_TAUS for t in a.taus):
        print(f"⚠⚠ EXPLORATORY `tau` {[t for t in a.taus if t in EXPLORATORY_TAUS]}: "
              f"NOT in the signed grid {list(TAUS)}. Report apart from the "
              f"pre-registered result.\n")
    jobs = [(t, tau, _cell_seed(ti, SEED_GRID.index(tau), i), a.one_sided)
            for tau in a.taus
            for ti, t in enumerate(tids)
            for i in range(a.draws)]
    n_per = len(tids) * a.draws
    # ⚠⚠ THE CONSTANT DEPENDS ON THE ESTIMATOR, and using the wrong one misstates
    # what the experiment can conclude. For a ONE-SIDED Bernoulli rate at p ~ 0.5 the
    # minimum detectable advantage is 1.243/sqrt(n); for the BALANCED accuracy this
    # design actually computes it is 0.879/sqrt(n), because averaging two independent
    # arms halves the variance. The pre-registration pins delta = 0.056 from the
    # BALANCED form at n = 250 per arm (0.879/sqrt(250) = 0.0556); the one-sided form
    # would give 0.0786 and would understate the experiment's resolution by 1.41x.
    # ⚠ The one-sided branch is kept because --one-sided still reports that rate, and
    # a delta quoted for one estimator does not apply to the other.
    delta = (1.243 if a.one_sided else 0.879) / math.sqrt(n_per)

    plan = R.plan_workers(R.array_bytes(st["model"].T, st["model"].K, n_arrays=4),
                          n_jobs=len(jobs))
    # ⚠ the balanced design draws from BOTH hypotheses: one Q_p release shared across
    # the three conditions, plus one Q_q release PER condition. So 4 releases a trial,
    # not 1, and the cost the pre-registration quotes must say which design it is for.
    per_trial = 1 if a.one_sided else 1 + 3
    print(f"{len(tids)} targets x {a.draws} draws x {len(a.taus)} tau = "
          f"{len(jobs):,} trials x {per_trial} release(s) = "
          f"{len(jobs) * per_trial:,} releases")
    if a.one_sided:
        print("⚠ ONE-SIDED MODE: reports P_p(Lambda>0) only. This is NOT comparable "
              "to p_succ_bound and no privacy claim may rest on it.")
    print(f"{n_per} trials per (tau, condition) per ARM -> minimum detectable "
          f"advantage delta = {delta:.4f} "
          f"({'one-sided 1.243' if a.one_sided else 'BALANCED 0.879'}/sqrt(n))")
    print(f"{plan}\n")

    res = R.run_jobs(_one_trial, jobs, workers=plan.workers, label="trials",
                     progress_every=max(1, len(jobs) // 10))

    return _report(res, a, delta, prereg, rev, tids, n_per)


if __name__ == "__main__":
    sys.exit(main())
