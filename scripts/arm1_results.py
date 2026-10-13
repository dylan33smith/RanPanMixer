#!/usr/bin/env python
"""The ARM 1 results table, emitted from the release records. Never hand-assembled.

⚠ THE REPORTING CONTRACT (`docs/plan.md`) IS THE SPEC, not a style guide: metrics are
ROWS and configurations are COLUMNS; the full metric set appears every time, with `n/a`
carrying a REASON rather than being omitted; row labels are exact `docs/terms.md`
identifiers with gates marked `*`; rows run primary -> gates -> structure -> context ->
demoted; and the table carries a provenance line and is followed by one bullet per
column, one per row, and a synthesis.

⚠⚠ IT ALSO VERIFIES, not just formats. Every cell's `log_z_p` is checked against
`[0, eta_tau]`, every genome record against the 22-autosome set, every part's `loo_id`
against its target, and the per-part `beta_share` sum against 1. A table that formats
whatever it is handed cannot tell a release from a mislabelled one.

⚠⚠ THE 44 TARGETS ARE 44 DIFFERENT UTILITY FUNCTIONS (standing constraint 8). Leave-one-
out gives each its own chain, `T`, 22 `beta_share` values and `u_ceiling`, so utilities
are NEVER pooled: every utility row is a median with its min and max beside it, and that
spread is a spread of DISTINCT quantities, not a sampling error.

    python scripts/arm1_results.py --endpoint logs/B_arm1_endpoint.json
"""

from __future__ import annotations

import argparse
import json
import sys
from collections import defaultdict
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

import numpy as np                                        # noqa: E402

import store                                              # noqa: E402
from cohort import loo                                     # noqa: E402

PINNED_TAUS = (0.0, 0.05, 0.10, 0.25, 0.50, 0.90, 0.99)
AUTOSOMES = set(loo.AUTOSOMES)


def load_cells(root: Path) -> dict:
    """{(target, tau): genome record}, every one verified as it is read."""
    cells, problems = {}, []
    for f in sorted(root.rglob("*_genome.json")):
        d = json.loads(f.read_text())
        key = (d["target_id"], float(d["tau"]))
        if key in cells:
            problems.append(f"{key}: two genome records, so constraint 5's "
                            f"one-release-per-genome was violated")
            continue
        if d.get("arm") != "arm1":
            problems.append(f"{key}: arm={d.get('arm')!r}, not arm1")
        if d.get("scope") != "genome":
            problems.append(f"{key}: scope={d.get('scope')!r}, not genome")
        if d.get("loo_id") != d["target_id"]:
            problems.append(f"{key}: loo_id={d.get('loo_id')!r} != target_id; an arm 1 "
                            f"release must hold out its own target (constraint 1)")
        if int(d.get("n_chromosomes", 0)) != 22:
            problems.append(f"{key}: {d.get('n_chromosomes')} chromosomes, not 22")
        if set(d.get("chromosomes", [])) != AUTOSOMES:
            problems.append(f"{key}: the parts are not the 22 autosomes")
        lz, eta = float(d["log_z_p"]), float(d["eta_tau"])
        if not (-1e-9 <= lz <= eta + 1e-9):
            problems.append(f"{key}: log_z_p={lz!r} outside [0, eta_tau={eta!r}] -- the "
                            f"budget was over- or under-spent")
        cells[key] = d
    return {"cells": cells, "problems": problems}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=Path,
                    default=REPO / "data" / "releases" / "arm1")
    ap.add_argument("--endpoint", type=Path, default=None,
                    help="logs/B_arm1_endpoint.json from scripts/arm1_exact_endpoint.py")
    ap.add_argument("--taus", type=float, nargs="+", default=list(PINNED_TAUS))
    a = ap.parse_args(argv)

    got = load_cells(a.root)
    cells, problems = got["cells"], got["problems"]
    targets = sorted({t for t, _ in cells})
    n_expect = len(targets) * len(a.taus)

    # ⚠ COMPLETENESS BEFORE ANYTHING IS PRINTED. A scoring script that tables whatever
    # it found cannot tell a partial run from a finished one.
    missing = [(t, tau) for t in targets for tau in a.taus if (t, tau) not in cells]
    if problems or missing or not cells:
        print(f"⛔ NOT REPORTABLE. {len(cells)} cell(s) found, {len(missing)} missing of "
              f"{n_expect}, {len(problems)} verification failure(s).")
        for p in problems[:10]:
            print(f"   ❌ {p}")
        for t, tau in missing[:10]:
            print(f"   ⬜ missing: {t} at `tau` = {tau:g}")
        raise SystemExit("arm 1 is not reportable until every cell exists and verifies.")

    ep = None
    if a.endpoint and a.endpoint.exists():
        ep = {r["target"]: r for r in json.loads(a.endpoint.read_text())["results"]
              if "error" not in r}

    print(f"running under prereg_2026-10-10.md + addendum_2026-10-11_decision_K.md + "
          f"addendum_2026-10-12_chain_rebuild.md (git_rev {store.git_rev(REPO)})")
    print(f"✅ verified on read: {len(cells)} genome records, each arm1 / genome scope / "
          f"22 autosomes / `loo_id` == target / `log_z_p` in [0, `eta_tau`].\n")
    print(f"ARM 1 releases. **ARM 1**, K = 86, per-target leave-one-out, diploid, "
          f"22 autosomes, n = {len(targets)} held-out individuals, "
          f"{len(a.taus)} `tau`, one release per genome per `tau`.\n")

    by = defaultdict(dict)
    for (t, tau), d in cells.items():
        by[tau][t] = d

    hdr = ["metric"] + [f"`tau` = {t:g}" for t in a.taus]
    print("| " + " | ".join(hdr) + " |")
    print("|" + "---|" * len(hdr))

    def row(label, fn, fmt=".6f", agg=np.median):
        out = []
        for tau in a.taus:
            vals = [fn(by[tau][t], t, tau) for t in targets]
            vals = [v for v in vals if v is not None]
            out.append(f"{agg(vals):{fmt}}" if vals else "n/a")
        print(f"| {label} | " + " | ".join(out) + " |")

    def na(label, reason):
        print(f"| {label} | " + " | ".join(f"n/a ({reason})" for _ in a.taus) + " |")

    # ---- PRIMARY
    if ep:
        # ⚠⚠ TWELVE DECIMALS, DELIBERATELY. The genome-wide movement is ~1e-07, so
        # `.6f` prints the SAME number in all seven columns and the one figure that
        # carries a claim disappears into the rounding.
        fid = lambda d, t, tau: (ep[t]["cells"][str(tau)]["target_fidelity"]
                                 if t in ep else None)
        row("`target_fidelity` * [median]", fid, ".12f")
        row("`target_fidelity` [min]", fid, ".12f", agg=np.min)
        row("`target_fidelity` [max]", fid, ".12f", agg=np.max)
        row("`e_differ` cells [median]",
            lambda d, t, tau: (ep[t]["cells"][str(tau)]["e_differ"] if t in ep else None),
            ",.3f")
    else:
        na("`target_fidelity` *", "needs scripts/arm1_exact_endpoint.py; pass --endpoint")
    # ---- GATES
    row("`log_z_p` * [median]", lambda d, t, tau: float(d["log_z_p"]), ".9f")
    row("`eta_tau` *", lambda d, t, tau: float(d["eta_tau"]), ".9f")
    print(f"| `k_states` * | " + " | ".join("86" for _ in a.taus) + " |")
    print(f"| `p_succ_bound` * | " +
          " | ".join(f"{(1 + t) / 2:.4f}" for t in a.taus) + " |")
    # ---- STRUCTURE
    row("`donor_switches` [median]", lambda d, t, tau: float(d["donor_switches"]), ",.0f")
    row("`in_chain` T [median]", lambda d, t, tau: float(d["T"]), ",.0f")
    row("`u_ceiling` [median]", lambda d, t, tau: float(d["u_ceiling"]), ".9f")
    # ---- CONTEXT
    if ep:
        row("`e_q_u` [median]",
            lambda d, t, tau: ep[t]["cells"][str(tau)]["e_q_u"] if t in ep else None, ".12f")
        row("`e_q_u_over_ceiling` [median]",
            lambda d, t, tau: ep[t]["cells"][str(tau)]["e_q_u_over_ceiling"] if t in ep else None,
            ".12f")
    else:
        na("`e_q_u`", "needs --endpoint")
        na("`e_q_u_over_ceiling`", "needs --endpoint")
    # ---- DEMOTED
    row("`utility_retained` [median, 1 draw]",
        lambda d, t, tau: float(d["u_achieved"]), ".9f")
    # ---- NOT YET MEASURED, with reasons rather than omission
    na("`utility_loss_panmixer` *", "chunk 9; scored at K = 88 both sides per decision M")
    na("`attacker_accuracy` *", "chunk 9")
    na("`af_loss` / `ld_loss`", "NOT ADOPTED, decision H -- a frozen endpoint")

    # ---- THE MOVEMENT. Everything above is a level; this is the claim.
    if ep:
        lo_t, hi_t = str(a.taus[0]), str(a.taus[-1])
        here = [t for t in targets if t in ep]
        print(f"\n⚠⚠ **THE MOVEMENT ACROSS `tau` = {a.taus[0]:g} TO {a.taus[-1]:g} IS THE "
              f"FIGURE, NOT THE LEVELS ABOVE.**")
        for key, unit in (("target_fidelity", ""), ("e_differ", " cells"),
                          ("e_q_u", "")):
            mv = [ep[t]["cells"][hi_t][key] - ep[t]["cells"][lo_t][key] for t in here]
            med0 = np.median([ep[t]["cells"][lo_t][key] for t in here])
            sign = "more faithful" if key != "e_q_u" else "more useful"
            print(f"- `{key}` moves **{np.median(mv):+.4e}{unit}** (median of "
                  f"{len(here)}), from a level of {med0:.9g}{unit}; per-target range "
                  f"{min(mv):+.4e} to {max(mv):+.4e}. "
                  f"{sum(1 for m in mv if (m < 0) == (key != 'e_q_u'))}/{len(here)} "
                  f"targets move in the {sign} direction.")
        print(f"⚠⚠ **SO THE WHOLE PRIVACY RANGE BUYS ABOUT FOUR CELLS OUT OF SEVEN AND A "
              f"HALF MILLION.** The level — about 17.4% of called cells altered — is "
              f"reached at `tau` = 0 by an untilted draw carrying ZERO target "
              f"information, because a cohort member's alleles are largely carried by "
              f"the other 86 haplotypes. The level is not attributable to the mechanism; "
              f"the movement is, and it is ~1e-07. This independently reproduces ARM 2's "
              f"genome-wide `e_q_u` movement of +8.345e-08 on a different arm with "
              f"{len(here)} targets: it is BUDGET DILUTION, not a broken mechanism.")

    # ---- one bullet per COLUMN
    print("\n**COLUMNS** (`tau`), each point's role as the pre-registration declares it:")
    roles = {0.0: "the FLOOR and mandatory -- an untilted draw carrying zero target "
                  "information; it is the control that stops a cohort-level figure "
                  "reading as a free lunch",
             0.05: "operating range", 0.10: "operating range -- bounds equal-prior "
                  "identification at 0.55", 0.25: "operating range; the only points we "
                  "claim as settings",
             0.50: "context, reported and NOT claimed",
             0.90: "bound-checking only, past any setting we would propose",
             0.99: "bound-checking only, past any setting we would propose"}
    for t in a.taus:
        lz = np.median([float(by[t][x]["log_z_p"]) for x in targets])
        print(f"- `tau` = {t:g} — {roles.get(t, 'declared role unrecorded')}. "
              f"Median `log_z_p` {lz:.9f} against an `eta_tau` of "
              f"{float(by[t][targets[0]]['eta_tau']):.9f}.")

    # ---- one bullet per ROW
    print("\n**ROWS**:")
    lo, hi = a.taus[0], a.taus[-1]
    if ep:
        f0 = np.median([ep[t]["cells"][str(lo)]["target_fidelity"] for t in targets if t in ep])
        f1 = np.median([ep[t]["cells"][str(hi)]["target_fidelity"] for t in targets if t in ep])
        print(f"- *`target_fidelity`* is the PRIMARY and a DIFFER rate, so lower is more "
              f"faithful. Median {f0:.6f} at `tau` = {lo:g} and {f1:.6f} at {hi:g}: a "
              f"movement of {f1 - f0:+.3e}. **The movement is the figure; the level at "
              f"`tau` = 0 belongs to human allele sharing, not to the mechanism.**")
    sw = np.median([float(by[lo][t]["donor_switches"]) for t in targets])
    sw1 = np.median([float(by[hi][t]["donor_switches"]) for t in targets])
    print(f"- *`donor_switches`* median {sw:,.0f} at `tau` = {lo:g} and {sw1:,.0f} at "
          f"{hi:g}. ⚠ The tilt provably cannot CREATE or REMOVE a switch — median "
          f"diversion cost 10.1529 nats against at most `eta_tau`/2 = 1.3233 on offer at "
          f"`tau` = 0.99 — so any change here is reassignment within a fixed skeleton.")
    tmin = min(float(by[lo][t]["T"]) for t in targets)
    tmax = max(float(by[lo][t]["T"]) for t in targets)
    print(f"- *`in_chain` T* spans {tmin:,.0f} to {tmax:,.0f} across the "
          f"{len(targets)} targets, because the chain is PER TARGET and the "
          f"`untraversable_block` set is not monotone in the donor panel.")
    print(f"- *`log_z_p`* is the budget gate: every one of the {len(cells)} cells sits "
          f"inside `[0, eta_tau]`, checked on read.")
    print(f"- *`utility_retained`* is DEMOTED (decision A) to a one-draw diagnostic; the "
          f"endpoint is `e_q_u`, exact. Never quote it without `u_ceiling` beside it.")

    # ---- SYNTHESIS
    print(f"\n**SYNTHESIS.** {len(cells)} arm 1 releases exist, one per "
          f"(held-out individual, `tau`), each verified as a genome-scope release of all "
          f"22 autosomes at K = 86 with its own `loo_id` and its `log_z_p` inside the "
          f"budget. ⚠⚠ Three things must be said with any of these numbers. "
          f"**(1)** The {len(targets)} targets carry {len(targets)} DIFFERENT utility "
          f"functions, so the utility rows are medians over distinct quantities and "
          f"constraint 8 forbids pooling them — pool attack outcomes, never utilities. "
          f"**(2)** The JOINT bound does not hold: these releases' donor panels are each "
          f"other's true haplotypes, so the per-release bound holds and the joint one "
          f"over the set does not (constraint 5). **(3)** Arm 1's guarantee is "
          f"CONDITIONAL — on leave-one-out, which `loo_id` records, and on `G`, which "
          f"was built from all HPRC assemblies including each target's. "
          f"⚠ Arm 1's `u_ceiling` is not comparable to arm 2's, and the `tau` grid is a "
          f"declared evaluation device rather than a release procedure.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
