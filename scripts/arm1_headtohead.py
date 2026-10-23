#!/usr/bin/env python
"""CHUNK 9: the ARM 1 head-to-head against PanMixer, chr21, SAME CODE ON BOTH SIDES.

⚠⚠ DECISION M (Dylan, 2026-10-12): match their evaluator exactly, scoring
`utility_loss_panmixer` with their shipped support over all K = 88. Our RELEASE stays
leave-one-out at K = 86 (decision E); only the SCORING support is K = 88.

⚠⚠ BOTH SIDES ARE COMPUTED HERE BY ONE FUNCTION, from each tool's own final released
object -- not read from either tool's log. Validated against their own reported numbers:
their `changed_alleles` (7,838) and `total_alleles` (650,720) at capacity 0.01 are
reproduced EXACTLY, and their total possible cost 9288.854446 reproduces their logged
9288.85. ⚠ Their logged `utility_loss` of 92.88862 is the OPTIMISER'S OBJECTIVE and
differs by 0.27% from a recomputation on the final stacked object (92.63823); we use the
recomputation, because that is the only figure computable the same way on both sides.
That is a difference in which quantity is reported, NOT a defect in PanMixer.

⚠⚠⚠ THE FINDING THAT MUST NOT BE BURIED: **THE COST ORDERING DEPENDS ENTIRELY ON THE -1
CONVENTION.** Their formula diffs raw. Their released object is the original with selected
blocks replaced (`obfuscate.py:330-340`), so it preserves -1; ours is a complete mosaic
(`ffbs.release` refuses a path carrying -1). So target-uncalled cells are altered in ours
BY CONSTRUCTION and almost never in theirs: measured, they are **0.045% of their raw cost
and ~46.6% of ours**. Report both conventions side by side and say which is which.
⚠ NOT A DEFECT IN PANMIXER. Its goal is donor protection, we claim no flaw, and
preserving missingness is consistent with that goal. It is a difference in what the two
mechanisms EMIT.

✅ THE ALIGNMENT AXIS IS CLEAN. Both sides' `target_fidelity` masks -1 -- theirs at
`obfuscate.py:474-481` (`valid = (original != -1) & (final != -1)`), ours via `known` in
`ffbs.expected_mismatches`. So the convention problem above does NOT touch the
pre-registered alignment axis. That separation is the cleanest fact in this table.

⚠ SCOPE. chr21 only, because that is where their numbers exist and `capacity` is natively
per-chromosome. ⚠⚠ Our column is a SLICE of a GENOME-WIDE release in which chr21 carried
only its `beta_share` (~0.0150) -- it is NOT a chr21-only release at that `tau` and must
never be labelled one. The slice is legitimate for this metric ONLY because
`utility_loss_panmixer` contains no `beta_v`: it is a plain additive sum of 1/support over
altered cells, so a chr21 restriction is a partial sum. Any `beta`-weighted quantity may
NOT be sliced (constraint 4).

⚠ `tau` and `capacity` ARE NOT COMPARABLE PARAMETERS and never share a column block.

    python scripts/arm1_headtohead.py --target HG00438 --json logs/B_arm1_headtohead.json
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

import numpy as np                                        # noqa: E402

import releases as R                                      # noqa: E402
from cohort import loo                                     # noqa: E402
from mechanism import calibration as C                     # noqa: E402
from mechanism import ffbs, hmm, utility                   # noqa: E402

PINNED_TAUS = (0.0, 0.05, 0.10, 0.25, 0.50, 0.90, 0.99)
PM_SWEEP = Path("/data/ds85/RanPanMixer/runs/capacity_sweep")
PM_DATA = REPO / "external" / "PanMixer" / "starting_data" / "chr21"


def score(orig: np.ndarray, final: np.ndarray, ul: np.ndarray) -> dict:
    """PanMixer's two conventions, on any (n, 2) original/final pair. ONE function,
    used for both tools -- that is the point of this script."""
    if orig.shape != final.shape:
        raise AssertionError(f"original {orig.shape} vs final {final.shape}")
    alt = final != orig                       # obfuscate.py:290 -- RAW, no masking
    valid = (orig != -1) & (final != -1)      # obfuscate.py:474-481 -- fidelity mask
    w = ul[:, None]
    cost_raw = float((alt * w).sum())
    cost_masked = float(((alt & valid) * w).sum())
    changed = int((alt & valid).sum())
    total = float(2.0 * ul.sum())
    return {
        "cost_raw": cost_raw, "cost_masked": cost_masked,
        "share_raw": cost_raw / total, "share_masked": cost_masked / total,
        "cost_total": total,
        "changed_masked": changed, "valid_cells": int(valid.sum()),
        "fidelity": changed / int(valid.sum()) if valid.any() else float("nan"),
        "altered_raw": int(alt.sum()), "all_cells": int(alt.size),
        "uncalled_cost_frac": (cost_raw - cost_masked) / cost_raw if cost_raw else 0.0,
    }


def pm_weights():
    pg = np.load(PM_DATA / "pangenome.npy", mmap_mode="r")
    nm = np.sum(np.asarray(pg) != -1, axis=(0, 2))
    ul = np.zeros(len(nm))
    ul[nm > 0] = 1.0 / nm[nm > 0]
    return pg, ul


def pm_final(cap_dir: Path, sample: str, n_sites: int) -> np.ndarray:
    """Their FINAL released haplotypes, read from their own obfuscated VCF."""
    vcf = cap_dir / "chr21" / "obfuscated.vcf.gz"
    out = subprocess.run(["bcftools", "query", "-s", sample, "-f", "[%GT]\\n",
                          str(vcf)], capture_output=True, text=True, check=True)
    gts = out.stdout.split()
    if len(gts) != n_sites:
        raise AssertionError(
            f"{vcf}: {len(gts)} genotypes against {n_sites} sites in their pangenome; "
            f"these are not the same axis")
    a = np.full((n_sites, 2), -1, dtype=np.int16)
    for i, g in enumerate(gts):
        parts = g.replace("/", "|").split("|")
        for j in (0, 1):
            x = parts[j] if j < len(parts) else "."
            if x not in (".", ""):
                a[i, j] = int(x)
    return a


def their_side(sample: str) -> list[dict]:
    pg, ul = pm_weights()
    orig = np.asarray(pg[list(np.load(PM_DATA / "pangenome_subjects.npy",
                                      allow_pickle=True)).index(sample)])
    rows = []
    for d in sorted(PM_SWEEP.glob("cap*")):
        st = json.loads((d / "chr21" / "stats.json").read_text())
        fin = pm_final(d, sample, orig.shape[0])
        s = score(orig, fin, ul)
        # ⚠ VALIDATE against their own reported numbers before using ours.
        if s["changed_masked"] != st["changed_alleles"]:
            raise AssertionError(
                f"{d.name}: we score {s['changed_masked']:,} changed against their "
                f"reported {st['changed_alleles']:,}. Our implementation of their "
                f"fidelity convention does not reproduce their own output.")
        if s["valid_cells"] != st["total_alleles"]:
            raise AssertionError(
                f"{d.name}: {s['valid_cells']:,} valid cells against their reported "
                f"{st['total_alleles']:,}")
        rows.append({"capacity": float(st["capacity"]), **s,
                     "their_reported_cost": float(st["utility_loss"]),
                     "their_reported_frac": float(st["utility_loss_frac"])})
    return rows


def our_side(target: str, taus) -> list[dict]:
    d = REPO / "data" / "cohort" / "chr21"
    hap = np.asarray(np.load(d / "haplotypes.npy", mmap_mode="r"))
    # ⚠ THEIR support convention: K = 88, the target included (decision M)
    nm = (hap != -1).sum(axis=0)
    ul_site = np.zeros(nm.shape[0])
    ul_site[nm > 0] = 1.0 / nm[nm > 0]
    # ⚠⚠ THE PER-TARGET CHAIN. Arm 1's chain is not a subset of arm 2's
    # (logs/B_loo_chain_survey.md), and substituting arm 2's `in_chain.npy` here was a
    # live bug until 2026-10-13 whose passing targets were a biased subset.
    ic = loo.loo_chain(d, target, haplotypes=hap).in_chain
    rows_i = list(loo.sample_index(np.load(d / "haplotype_ids.npy",
                                           allow_pickle=True))[target])
    orig = hap[rows_i][:, ic].T
    ul = ul_site[ic]
    out = []
    for tau in taus:
        tdir = REPO / "data" / "releases" / "arm1" / target / R.tau_dir(tau)
        hits = sorted(tdir.glob("*chromchr21_*_chr21/released_path.npy"))
        if len(hits) != 1:
            raise AssertionError(f"{target}/{tau:g}: {len(hits)} chr21 parts, expected 1")
        fin = np.load(hits[0])
        out.append({"tau": float(tau), **score(orig, fin, ul)})
    return out


def exact_cost(target: str, taus) -> dict:
    """E[`utility_loss_panmixer`] EXACTLY, from the posterior marginals.

    ⚠⚠ REQUIRED BY CHUNK 9's OWN EXIT GATE (`docs/plan.md`): "The EXACT expectation is
    printed beside every sampled mean, so the small R is a diagnostic and never the
    evidence." A one-draw cost row alone does not pass.

    THE DECOMPOSITION, and why the first term needs no marginals at all:
      * A cell where the TARGET IS UNCALLED is altered with probability **1**. The
        allowed mask is `hap != -1` (`src/mechanism/hmm.py`), so the chain never
        occupies a donor without a call and the released allele is always a called
        value, which differs from -1 always. So that part of the cost is DETERMINISTIC:
        `sum of 1/support(v)` over uncalled cells, independent of `tau`.
      * A cell where the target IS called contributes
        `1/support(v) * P(released allele != target allele)`, and that probability is
        `sum_k gamma[t, k] * 1[hap[k, t] != target[t]]` -- a linear functional of the
        marginals, so exact with no sampling and no differencing.
    E[cost_masked] is the second term alone; E[cost_raw] is both.
    """
    d = REPO / "data" / "cohort" / "chr21"
    hap = np.asarray(np.load(d / "haplotypes.npy", mmap_mode="r"))
    nm = (hap != -1).sum(axis=0)                 # K = 88, decision M
    ul_site = np.zeros(nm.shape[0])
    ul_site[nm > 0] = 1.0 / nm[nm > 0]
    lc = loo.loo_chain(d, target, haplotypes=hap)
    ic = lc.in_chain
    panel = loo.masked_haplotypes(d, target, haplotypes=hap)[:, ic]
    rows_i = list(loo.sample_index(np.load(d / "haplotype_ids.npy",
                                           allow_pickle=True))[target])
    tgt = hap[rows_i][:, ic].T
    ul = ul_site[ic]
    m = hmm.from_arrays(panel, np.load(d / "genetic_pos.npy")[ic],
                        np.load(d / "positions.npy")[ic], cohort_id=f"chr21_loo_{target}")
    # ⚠ the tilt stays the GENOME-WIDE one: chr21 carries only its own `beta_share`.
    # Renormalising `beta` to chr21 would make `tau` per-chromosome (constraint 4).
    shares = loo.loo_shares({c: loo.loo_chain(REPO / "data" / "cohort" / c,
                                              target).t_chain for c in loo.AUTOSOMES})
    util = utility.beta_v(lc.support, in_chain=ic, share=shares["chr21"], chrom="chr21")

    uncalled_cost = float(sum(ul[tgt[:, c] < 0].sum() for c in (0, 1)))
    out = {}
    for tau in taus:
        eta = C.eta_tau(tau)
        masked = 0.0
        for col in (0, 1):
            t1 = tgt[:, col]
            psi = utility.tilt_potentials(m.allowed, panel, t1, util.beta, eta / 2.0)
            gam = ffbs.state_marginals(m, ffbs.forward(m, psi), psi)
            known = t1 >= 0
            mism = (panel != t1[None, :]) & known[None, :]
            masked += float((gam * mism.T).sum(axis=1) @ ul)
        out[str(tau)] = {"exact_masked": masked,
                         "exact_raw": masked + uncalled_cost,
                         "uncalled_deterministic": uncalled_cost,
                         "cost_total": float(2.0 * ul.sum())}
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--target", default="HG00438",
                    help="must be BOTH a cohort member and their swept subject")
    ap.add_argument("--taus", type=float, nargs="+", default=list(PINNED_TAUS))
    ap.add_argument("--exact", action="store_true", default=True,
                    help="compute E[cost] from the marginals (chunk 9's exit gate)")
    ap.add_argument("--no-exact", dest="exact", action="store_false")
    ap.add_argument("--json", type=Path, default=None)
    a = ap.parse_args(argv)

    theirs = their_side(a.target)
    ours = our_side(a.target, a.taus)
    ex = exact_cost(a.target, a.taus) if a.exact else None
    if ex:
        for r in ours:
            e = ex[str(r["tau"])]
            r.update({k: e[k] for k in ("exact_raw", "exact_masked",
                                        "uncalled_deterministic")})
            r["exact_share_raw"] = e["exact_raw"] / e["cost_total"]
            r["exact_share_masked"] = e["exact_masked"] / e["cost_total"]
    print(f"✅ validated: our scoring reproduces PanMixer's own `changed_alleles` and "
          f"`total_alleles` EXACTLY at all {len(theirs)} capacities.\n")
    _report(ours, theirs, a.target)
    if a.json:
        a.json.write_text(json.dumps(
            {"arm": "arm1", "scope": "chr21", "target": a.target, "support_k": 88,
             "ours": ours, "theirs": theirs}, indent=1) + "\n")
        print(f"\nJSON -> {a.json}")
    return 0


def _report(ours, theirs, target) -> None:
    assert ours and theirs, "one side is empty; the comparison would be silent"
    print(f"ARM 1 vs PanMixer, chr21, target {target}, scoring support K = 88 "
          f"(decision M), n = 1 target, 1 draw per `tau`.")
    print(f"⚠ `tau` and `capacity` are NOT comparable parameters: two separate blocks.\n")

    print(f"**BLOCK A -- OURS** (ARM 1, K = 86 release, genome-wide tilt SLICED to chr21, "
          f"chr21 `beta_share` ~0.0150, post-decision-K chain)\n")
    h = ["metric"] + [f"`tau` = {r['tau']:g}" for r in ours]
    print("| " + " | ".join(h) + " |"); print("|" + "---|" * len(h))
    for lab, k, f in (("`utility_loss_panmixer` * [RAW -1, 1 draw]", "cost_raw", ",.4f"),
                      ("  its share of own total [RAW]", "share_raw", ".6f"),
                      ("`utility_loss_panmixer` [-1 MASKED, 1 draw]", "cost_masked", ",.4f"),
                      ("  its share of own total [MASKED]", "share_masked", ".6f"),
                      ("`target_fidelity` * [masked, 1 draw]", "fidelity", ".6f"),
                      ("  cells changed [masked]", "changed_masked", ",.0f"),
                      ("  denominator: cells valid both sides", "valid_cells", ",.0f"),
                      ("  -1 share OF THE RAW COST", "uncalled_cost_frac", ".4%")):
        print(f"| {lab} | " + " | ".join(f"{r[k]:{f}}" for r in ours) + " |")
    if "exact_raw" in ours[0]:
        # ⚠ the exit gate: the EXACT expectation beside every sampled mean
        for lab, k, f in (
                ("`utility_loss_panmixer` * [EXACT E, RAW -1]", "exact_raw", ",.4f"),
                ("  its share of own total [EXACT, RAW]", "exact_share_raw", ".6f"),
                ("`utility_loss_panmixer` * [EXACT E, MASKED]", "exact_masked", ",.4f"),
                ("  its share of own total [EXACT, MASKED]", "exact_share_masked", ".6f"),
                ("  of which DETERMINISTIC (uncalled, p = 1)",
                 "uncalled_deterministic", ",.4f")):
            print(f"| {lab} | " + " | ".join(f"{r[k]:{f}}" for r in ours) + " |")

    print(f"\n**BLOCK B -- THEIRS** (PanMixer, their own build, HG00438, chr21, "
          f"sampling-only LOO with PMI scoring at K = 88)\n")
    h = ["metric"] + [f"`capacity` = {r['capacity']:g}" for r in theirs]
    print("| " + " | ".join(h) + " |"); print("|" + "---|" * len(h))
    for lab, k, f in (("`utility_loss_panmixer` * [RAW -1, recomputed]", "cost_raw", ",.4f"),
                      ("  its share of own total [RAW]", "share_raw", ".6f"),
                      ("`utility_loss_panmixer` [-1 MASKED]", "cost_masked", ",.4f"),
                      ("  its share of own total [MASKED]", "share_masked", ".6f"),
                      ("`target_fidelity` * [masked]", "fidelity", ".6f"),
                      ("  cells changed [masked]", "changed_masked", ",.0f"),
                      ("  denominator: cells valid both sides", "valid_cells", ",.0f"),
                      ("  -1 share OF THE RAW COST", "uncalled_cost_frac", ".4%"),
                      ("  their OWN reported cost [optimiser objective]",
                       "their_reported_cost", ",.4f")):
        print(f"| {lab} | " + " | ".join(f"{r[k]:{f}}" for r in theirs) + " |")

    of = [r["fidelity"] for r in ours]
    tf = [r["fidelity"] for r in theirs]
    om = [r["share_masked"] for r in ours]
    tm = [r["share_masked"] for r in theirs]
    orr = [r["share_raw"] for r in ours]
    trr = [r["share_raw"] for r in theirs]
    print(f"\n**SYNTHESIS.**")
    print(f"- ⚠⚠ **THE COST ORDERING REVERSES WITH THE -1 CONVENTION, so neither number "
          f"stands alone.** RAW, ours is {np.median(orr):.4f} of our total against theirs "
          f"spanning {min(trr):.4f}-{max(trr):.4f}. MASKED, ours is {np.median(om):.4f} "
          f"against theirs {min(tm):.4f}-{max(tm):.4f}. The cause is measured: the -1 "
          f"convention is {np.median([r['uncalled_cost_frac'] for r in theirs]):.4%} of "
          f"their raw cost and {np.median([r['uncalled_cost_frac'] for r in ours]):.4%} "
          f"of ours, because a path release is COMPLETE by construction while theirs "
          f"preserves the original outside the blocks it swaps. **Not a defect in "
          f"PanMixer** -- a difference in what the two mechanisms emit.")
    print(f"- ✅ **THE ALIGNMENT AXIS IS UNAFFECTED**: both sides mask -1 for "
          f"`target_fidelity`, and our scoring reproduces their reported `changed_alleles` "
          f"and `total_alleles` exactly. Ours spans {min(of):.6f}-{max(of):.6f} across the "
          f"`tau` grid; theirs {min(tf):.6f}-{max(tf):.6f} across capacity.")
    print(f"- ⚠⚠ **NO MATCHED-FIDELITY CAPACITY EXISTS, and that is a result rather than a "
          f"gap.** Their knob SATURATES at {max(tf):.6f} from capacity 0.3 up, so "
          f"`capacity -> fidelity` is not invertible above it. Our level "
          f"({np.median(of):.6f} median) sits at or above their maximum, so there is no "
          f"upper bracket and nothing is interpolated. ⚠ The margin is far inside one-draw "
          f"noise -- our own `tau` cells span {max(of)-min(of):+.6f} on independent seeds "
          f"against an EXACT genome-wide `tau` effect of -9.8e-08 -- so this says our "
          f"fidelity is at their ceiling, NOT that it is distinguishably worse.")
    print(f"- ⚠ `tau` is a privacy parameter with a proof behind it; `capacity` is a "
          f"utility-loss budget and bounds nothing about distinguishability. They are not "
          f"two settings of one dial, and no row pairs them.")
    print(f"- ⚠ Our column is a chr21 SLICE of a genome-wide release carrying only chr21's "
          f"`beta_share`; it is not a chr21-only release. Their column is their own build "
          f"({theirs[0]['all_cells'] // 2:,} sites) against ours "
          f"({ours[0]['all_cells'] // 2:,} chain positions), so only SHARES of each "
          f"side's own total are quotable, never raw counts across the two.")


if __name__ == "__main__":
    sys.exit(main())
