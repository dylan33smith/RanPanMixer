#!/usr/bin/env python
"""Emit `new_haplotypes`: our arm 1 release in the format PanMixer's evaluators consume.

⚠⚠ THIS IS ARM 1 PLUMBING AND IT IS NOT OPTIONAL. `docs/terms.md:new_haplotypes` defines
it as "PanMixer's integration artifact: the released allele vector for one target, shape
(n_sites, 2) int16, -1 = missing, row-aligned to that chromosome's pangenome VCF record
order. Every PanMixer attack and utility evaluator consumes exactly this file", and marks
it PLANNED and UNIMPLEMENTED. The head-to-head's exit gate in `docs/plan.md` reads "Both
emit `new_haplotypes`; the shared evaluators run unchanged" -- so without this file their
evaluators cannot run on our release at all, and chunk 9 compared utility only.

⚠⚠⚠ THE HAZARD IS SILENT AND IT IS NAMED IN THE TERM ITSELF: "It is positional -- a
row-order mismatch silently scores the wrong variants rather than erroring. The shape
assertion in PanMixer's gap-score evaluators is VACUOUS (both sides equal
site_mask.sum() by construction), so nothing upstream catches this." So this script
refuses to write unless an alignment control passes, and ships a `sites.tsv` beside every
array as the term requires.

THE AXES DO NOT MATCH AND THAT IS THE WHOLE DIFFICULTY.
  * our site axis:  340,849 chr21 records
  * their site axis: 340,824 chr21 records
  * shared by POSITION: 339,744
  * our RELEASE covers only our CHAIN -- 332,521 positions for a typical arm 1 target,
    and the chain is PER TARGET (`logs/B_loo_chain_survey.md`).
So a typical emitted file carries a released allele at ~331,400 of their 340,824 rows.

⚠⚠ THE REMAINING ROWS ARE -1 (MISSING), NOT THE TARGET'S TRUE ALLELE, AND THE CHOICE IS A
PRIVACY STATEMENT. PanMixer's own stacker leaves the ORIGINAL in place outside the blocks
it swaps (`obfuscate.py:330-340`), so filling our unreleased rows the same way would hand
their linkage attack ~9,400 of the target's TRUE alleles -- alleles our mechanism never
released. That would be a different object from the one our guarantee covers, and it would
make our measured privacy look WORSE than the mechanism is while claiming something we did
not do. Emitting -1 says what is true: at those positions we released nothing.
⚠ Their evaluators mask on both sides being called, so -1 rows drop out of their
statistics rather than scoring as matches. That is the conservative direction for us and
it must be stated with any number computed from this file.

    python scripts/emit_new_haplotypes.py --target HG00438 --tau 0.25
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

import numpy as np                                        # noqa: E402

import releases as R                                      # noqa: E402
from cohort import loo                                     # noqa: E402

PM_DATA = REPO / "external" / "PanMixer" / "starting_data" / "chr21"
OUT_ROOT = REPO / "data" / "panmixer_integration"
PINNED_TAUS = (0.0, 0.05, 0.10, 0.25, 0.50, 0.90, 0.99)
CHROM = "chr21"


def axis_map():
    """Row indices into THEIR axis and OUR site axis for the positions both share."""
    our_pos = np.load(REPO / "data" / "cohort" / CHROM / "positions.npy")
    their_pos = np.load(PM_DATA / "pangenome_positions.npy")
    shared = np.intersect1d(our_pos, their_pos)
    oi = np.searchsorted(our_pos, shared)
    ti = np.searchsorted(their_pos, shared)
    # ⚠ searchsorted returns an insertion point even when the value is absent, so the
    # equality check below is what makes this a mapping rather than a guess.
    assert np.array_equal(our_pos[oi], shared), "our index does not recover the positions"
    assert np.array_equal(their_pos[ti], shared), "their index does not recover the positions"
    return our_pos, their_pos, shared, oi, ti


def alignment_control(oi, ti, their_pos, target: str, hap: np.ndarray) -> dict:
    """THE GATE. Map the target's own original onto their axis and compare it to THEIR
    original. The two builds are known to agree at ~98.66% of jointly-called cells
    (`logs/B_arm1_symmetric.md`); a row-order mismatch collapses that toward the rate two
    unrelated rows agree at, which is what gives this check its power.

    ⚠ This is NOT a check that the two builds agree -- they do not, and that disagreement
    is a stated limitation. It is a check that we are comparing the SAME ROWS.
    """
    pg = np.load(PM_DATA / "pangenome.npy", mmap_mode="r")
    subj = list(np.load(PM_DATA / "pangenome_subjects.npy", allow_pickle=True))
    their_orig = np.asarray(pg[subj.index(target)])[ti]
    rows = list(loo.sample_index(np.load(
        REPO / "data" / "cohort" / CHROM / "haplotype_ids.npy", allow_pickle=True))[target])
    our_orig = hap[rows][:, oi].T
    both = (our_orig != -1) & (their_orig != -1)
    agree = float((our_orig[both] == their_orig[both]).mean())
    # the same comparison with their rows shifted by one -- the misalignment this guards
    shifted = float((our_orig[both][:-1] == their_orig[both][1:]).mean())
    return {"aligned": agree, "shifted_by_one": shifted,
            "n_compared": int(both.sum())}


def emit(target: str, tau: float, out_root: Path, force: bool = False) -> dict:
    our_pos, their_pos, shared, oi, ti = axis_map()
    hap = np.asarray(np.load(REPO / "data" / "cohort" / CHROM / "haplotypes.npy",
                             mmap_mode="r"))

    ctl = alignment_control(oi, ti, their_pos, target, hap)
    # ⚠ THE REFUSAL. 0.90 sits far below the ~0.9866 two correctly aligned builds reach
    # and far above what a one-row shift leaves, so it separates the two cases without
    # pretending the builds agree perfectly.
    if ctl["aligned"] < 0.90 and not force:
        raise SystemExit(
            f"ALIGNMENT CONTROL FAILED for {target}: jointly-called agreement "
            f"{ctl['aligned']:.6%} against the ~98.66% two aligned builds reach. "
            f"Refusing to write -- a row-order mismatch scores the WRONG VARIANTS "
            f"silently in their evaluators, whose shape assertion cannot catch it.")

    # our per-target chain, and the release on it
    ic = loo.loo_chain(REPO / "data" / "cohort" / CHROM, target, haplotypes=hap).in_chain
    tdir = (REPO / "data" / "releases" / "arm1" / target
            / R.tau_dir(tau))
    hits = sorted(tdir.glob(f"*chrom{CHROM}_*_{CHROM}/released_path.npy"))
    if len(hits) != 1:
        raise AssertionError(f"{target}/{tau:g}: {len(hits)} chr21 parts, expected 1")
    rel = np.load(hits[0])                                  # (T_chain, 2) on OUR chain

    chain_of_site = np.full(our_pos.shape[0], -1, dtype=np.int64)
    chain_of_site[ic] = np.arange(int(ic.sum()))
    ci = chain_of_site[oi]                                  # chain row per shared position
    have = ci >= 0                                          # shared AND on our chain

    out = np.full((their_pos.shape[0], 2), -1, dtype=np.int16)
    out[ti[have]] = rel[ci[have]]

    # ⚠⚠ ALLELE CODES ARE INDICES INTO EACH BUILD'S OWN ALT LIST, so a code valid on our
    # axis can be out of range on theirs. Measured 2026-10-13: 5 cells across our 7 taus
    # exceed their `num_alleles` — e.g. row 269,520 carries our code 9 where their
    # encoding declares 4 — at multi-allelic sites where our build called more ALTs than
    # theirs. Their own original has 0 such cells, so this is ours alone.
    # `tools/common/verify_genotypes.py` asserts code < num_alleles, so an unfixed file
    # fails their OWN pre-flight; rendering it to VCF for Beagle would emit an allele
    # index with no corresponding ALT.
    # We MASK rather than remap: our release carries an allele their encoding has no name
    # for, and -1 says exactly that. It is also the conservative direction, since their
    # evaluators drop both-called cells.
    na = np.load(PM_DATA / "num_alleles.npy")
    oor = (out >= na[:, None]) & (out >= 0)
    n_oor = int(oor.sum())
    out[oor] = -1

    d = out_root / "arm1" / target / R.tau_dir(tau) / CHROM
    d.mkdir(parents=True, exist_ok=True)
    np.save(d / "new_haplotypes.npy", out)
    # ⚠ the sidecar docs/terms.md requires: "Ship a `sites.tsv` beside any emitted copy
    # and verify it." Row index, position, and whether we released there.
    with (d / "sites.tsv").open("w") as f:
        f.write("row\tpos\treleased\n")
        rel_mask = np.zeros(their_pos.shape[0], dtype=bool)
        rel_mask[ti[have]] = True
        for i in range(their_pos.shape[0]):
            f.write(f"{i}\t{int(their_pos[i])}\t{int(rel_mask[i])}\n")

    return {"target": target, "tau": tau, "path": str(d / "new_haplotypes.npy"),
            "shape": list(out.shape), "released_rows": int(have.sum()),
            "out_of_range_masked": n_oor,
            "missing_rows": int(their_pos.shape[0] - have.sum()),
            "their_sites": int(their_pos.shape[0]), **ctl}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--target", default="HG00438")
    ap.add_argument("--tau", type=float, nargs="+", default=list(PINNED_TAUS))
    ap.add_argument("--out", type=Path, default=OUT_ROOT)
    ap.add_argument("--force", action="store_true",
                    help="write even if the alignment control fails (do not use)")
    a = ap.parse_args(argv)

    print(f"emitting `new_haplotypes` for {a.target}, {CHROM}, "
          f"{len(a.tau)} `tau` value(s)\n")
    rows = []
    for tau in a.tau:
        r = emit(a.target, tau, a.out, a.force)
        rows.append(r)
        print(f"  tau = {tau:<5g} -> {r['shape']}, {r['released_rows']:,} released rows, "
              f"{r['missing_rows']:,} at -1"
              + (f"  ⚠ {r['out_of_range_masked']} cell(s) masked as OUT OF RANGE for "
                 f"their allele encoding" if r['out_of_range_masked'] else ""))
    c = rows[0]
    print(f"\n✅ ALIGNMENT CONTROL: jointly-called agreement **{c['aligned']:.6%}** over "
          f"{c['n_compared']:,} cells, against **{c['shifted_by_one']:.6%}** with their "
          f"rows shifted by ONE. The gap between those two is what gives the control its "
          f"power -- their own shape assertion is vacuous and cannot catch a shift.")
    print(f"⚠ {c['missing_rows']:,} of {c['their_sites']:,} rows carry -1: positions our "
          f"release does not cover (off-chain, or absent from their axis). They are NOT "
          f"filled with the target's true allele -- doing so would hand a linkage attack "
          f"alleles we never released.")
    print(f"⚠ Their evaluators mask on both sides being called, so those rows drop out "
          f"rather than scoring as matches. State that with any number computed from "
          f"this file.")
    print(f"⚠ `sites.tsv` is written beside each array, as docs/terms.md requires.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
