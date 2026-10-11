#!/usr/bin/env python
"""Check that built cohort artifacts are CORRECT, not merely present.

⚠ WHY THIS IS SEPARATE FROM THE DRIVER. `build_all_cohorts.py` asserts COMPLETENESS —
every chromosome ran, every expected file exists and is non-empty. That is necessary
and says nothing about whether the contents agree with each other. A chromosome whose
`support.npy` is row-misaligned with its `haplotypes.npy` passes every completeness
check and then silently scores the wrong variants, because the arrays are positional
and a shift produces numbers rather than an error.

Each check below corresponds to a way this pipeline can actually go wrong, not to a
property that is true by construction. A check that cannot fail is decoration.

    python scripts/verify_cohorts.py [--chroms 1 2 21]
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]


PER_SITE = ("positions.npy", "support.npy", "an_info.npy", "genetic_pos.npy",
            "in_chain.npy", "cm_source.npy")


def check(c: int, root: Path = None) -> list[str]:
    d = (root or REPO / "data" / "cohort") / f"chr{c}"
    if not (d / "manifest.json").exists():
        return [f"chr{c}: not built"]

    # 0. every artifact this checker reads must EXIST. ⚠ A missing artifact is a
    #    finding, not a crash: a traceback here stops the other 21 chromosomes from
    #    being checked at all, which is exactly the "found nothing vs ran nothing"
    #    confusion the driver was written to avoid. The case this fires on is real --
    #    a chromosome built BEFORE decision C has no `cm_source.npy`, because the
    #    coordinate cut had nothing to record.
    missing = [f for f in ("haplotypes.npy",) + PER_SITE if not (d / f).exists()]
    if missing:
        rule = json.loads((d / "manifest.json").read_text()).get("chain_rule")
        note = ("" if rule == "declared_rate_fill" else
                f" -- manifest says chain_rule={rule!r}, so this is a PRE-DECISION-C "
                f"build and must be rebuilt before any genome-wide number is quoted")
        return [f"chr{c}: missing {', '.join(missing)}{note}"]

    bad: list[str] = []
    hap = np.load(d / "haplotypes.npy")
    pos = np.load(d / "positions.npy")
    sup = np.load(d / "support.npy")
    an = np.load(d / "an_info.npy")
    cm = np.load(d / "genetic_pos.npy")
    ic = np.load(d / "in_chain.npy").astype(bool)
    src = np.load(d / "cm_source.npy")
    K, n = hap.shape

    # 1. every per-site array is the SAME length. A positional array one row short
    #    silently shifts every downstream join.
    for nm, arr in (("positions", pos), ("support", sup), ("an_info", an),
                    ("genetic_pos", cm), ("in_chain", ic), ("cm_source", src)):
        if len(arr) != n:
            bad.append(f"chr{c}: {nm} has {len(arr)} rows, haplotypes has {n}")
    # ⚠ RETURN EARLY on a length disagreement. Every check below is POSITIONAL --
    # it compares arrays element by element -- so once two of them have different
    # lengths the rest are not merely uninformative, they raise on broadcast and
    # take the whole run down with them. Found by the mutation test: truncating
    # `support.npy` by one row crashed check 2 instead of reporting check 1.
    if bad:
        return bad

    # 2. support(v) equals INFO/AN. Two independent routes to one quantity -- one
    #    parsed from genotypes, one read from a VCF field -- so they COULD disagree.
    #    This is the cheapest real correctness check in the pipeline.
    called = (hap != -1).sum(axis=0)
    if not np.array_equal(called, sup):
        bad.append(f"chr{c}: support.npy disagrees with the called-genotype count at "
                   f"{int((called != sup).sum())} sites")
    if not np.array_equal(sup, an):
        bad.append(f"chr{c}: support.npy disagrees with INFO/AN at "
                   f"{int((sup != an).sum())} sites")

    # 3. positions NON-DECREASING -- not strictly increasing.
    #    ⚠ THIS CHECK WAS WRONG WHEN FIRST WRITTEN and chr9 falsified it. It demanded
    #    STRICT increase and flagged chr9's position 37,752,316, where two records
    #    genuinely share a start: a parent bubble (lv=1, 302 bp REF, 4 alts) and a
    #    1 bp deletion NESTED inside it (lv=2, GG->G). Nesting is not an anomaly here
    #    -- 10-18% of records on every autosome carry lv > 0 -- so a shared start
    #    coordinate is legal `vg deconstruct` output, and refusing it would have
    #    condemned a correct chromosome.
    #    What is NOT legal is a DECREASE: it makes a genetic distance negative, which
    #    puts the Li--Stephens switch probability outside [0,1], and `np.interp`
    #    would not complain. That is the error worth catching, so that is what this
    #    checks. Ties are counted and REPORTED in the table rather than failed,
    #    because a count that jumps is worth seeing even when it is not an error.
    dec = int((np.diff(pos) < 0).sum())
    if dec:
        bad.append(f"chr{c}: positions DECREASE at {dec} step(s) -- a negative genetic "
                   f"distance puts the switch probability outside [0,1]")

    # 4. no NaN cM inside the chain, and cM non-decreasing there. A negative step
    #    puts the switch probability outside [0,1].
    if not np.isfinite(cm[ic]).all():
        bad.append(f"chr{c}: {int((~np.isfinite(cm[ic])).sum())} chain positions have "
                   f"no cM -- the declared-rate fill did not run")
    dd = np.diff(cm[ic])
    if (dd < -1e-12).any():
        bad.append(f"chr{c}: cM decreases at {int((dd < -1e-12).sum())} chain steps")

    # 5. THE CHAIN HAS EXACTLY TWO EXCLUSIONS, and both are recomputed here from the
    #    data rather than trusted. Decision C retired the coordinate cut; decision K
    #    adds the untraversable blocks. Any THIRD exclusion is a bug or an
    #    accidentally reintroduced span.
    #
    #    ⚠ The second half is recomputed INDEPENDENTLY of the code that wrote the
    #    mask, which is what gives this check power in both directions: it catches an
    #    arbitrary extra exclusion AND a block that should have been dropped and was
    #    not. Comparing the mask to the manifest's own record of what it dropped would
    #    only confirm the writer agrees with itself.
    usable = sup > 0
    if (ic & ~usable).any():
        bad.append(f"chr{c}: {int((ic & ~usable).sum())} position(s) with support == 0 "
                   f"are IN the chain -- the sampler has nothing to emit there")
    excluded_usable = usable & ~ic
    # ⚠⚠ RECOMPUTED UNCONDITIONALLY, and the first version of this check got that
    # wrong: it only recomputed when something was ALREADY excluded, so a chain that
    # was never rebuilt after decision K -- exactly `in_chain == (support > 0)`, with
    # nothing excluded -- skipped the branch entirely and passed. That is the single
    # case this check exists for, and the forward pass would then raise at RELEASE
    # time instead. Found by the test for the opposite direction, which is why both
    # directions are tested.
    if True:
        import sys as _sys
        _sys.path.insert(0, str(REPO / "src"))
        from cohort import genmap as _G
        allowed = np.asarray(np.load(d / "haplotypes.npy", mmap_mode="r")) != -1
        # the blocks decision K would drop from the support-only chain
        sites = np.flatnonzero(usable)
        blocks = _G.untraversable_blocks(cm[sites], allowed[:, sites])
        should_drop = np.zeros_like(ic)
        for a, b in blocks:
            should_drop[sites[a:b + 1]] = True
        if not np.array_equal(excluded_usable, should_drop):
            over = int((excluded_usable & ~should_drop).sum())
            under = int((should_drop & ~excluded_usable).sum())
            bad.append(
                f"chr{c}: the usable positions excluded from the chain are NOT the "
                f"untraversable blocks -- {over} excluded that decision K would keep, "
                f"{under} kept that it would drop. A third exclusion rule has crept "
                f"in, or the chain was not rebuilt after decision K.")

    # 6. cm_source marks declared-rate positions, and they are only at the ENDS.
    #    numpy.interp cannot leave an interior gap, so an interior declared-rate
    #    position means the fill was applied to something it should not have been.
    if src.max() > 0:
        idx = np.flatnonzero(src == 0)
        if len(idx) and (src[idx[0]:idx[-1] + 1] != 0).any():
            bad.append(f"chr{c}: declared-rate positions appear INSIDE the mapped "
                       f"span, which the fill should never produce")

    # 7. one axis digest across every artifact. Mixing two builds' arrays is the
    #    error `store.check_axis` exists for, and the manifest already records it.
    man = json.loads((d / "manifest.json").read_text())
    if not man.get("axis_consistent", False):
        bad.append(f"chr{c}: manifest reports axis_consistent=False "
                   f"({man.get('axis_digests_seen')})")
    if man.get("chain_rule") != "declared_rate_fill":
        bad.append(f"chr{c}: built under {man.get('chain_rule')!r}, not the "
                   f"post-decision-C declared_rate_fill")
    return bad


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--chroms", type=int, nargs="+", default=list(range(1, 23)))
    ap.add_argument("--root", type=Path, default=REPO / "data" / "cohort",
                    help="cohort tree to check; used by the mutation test")
    a = ap.parse_args(argv)

    built, problems, rows = [], [], []
    for c in a.chroms:
        d = a.root / f"chr{c}"
        if not (d / "manifest.json").exists():
            continue
        built.append(c)
        bad = check(c, a.root)
        problems.extend(bad)
        # ⚠ The row is built from the SAME files `check` reads, so it has the same
        # way of being absent. Guard it here too rather than relying on `check`
        # having returned early -- two readers of one artifact is two crash sites,
        # and fixing only the first is how the traceback came back identical.
        try:
            hap = np.load(d / "haplotypes.npy", mmap_mode="r")
            ic = np.load(d / "in_chain.npy").astype(bool)
            src = np.load(d / "cm_source.npy")
        except FileNotFoundError:
            rows.append({"chrom": c, "records": None, "K": None, "T": None,
                         "declared": None, "shared_pos": None, "ok": False})
            continue
        rows.append({"chrom": c, "records": hap.shape[1], "K": hap.shape[0],
                     "T": int(ic.sum()),
                     "declared": int((src[ic] == 1).sum()),
                     "shared_pos": int((np.diff(np.load(d / "positions.npy")) == 0).sum()),
                     "ok": not bad})

    print(f"| chr | K | records | T (chain) | on declared rate | shared pos | checks |")
    print("|" + "---|" * 7)
    for r in sorted(rows, key=lambda r: r["chrom"]):
        if r["T"] is None:
            print(f"| {r['chrom']} | ? | ? | ? | ? | ? | **FAIL** |")
            continue
        print(f"| {r['chrom']} | {r['K']} | {r['records']:,} | {r['T']:,} | "
              f"{r['declared']:,} ({100*r['declared']/r['T']:.2f}%) | "
              f"{r['shared_pos']} | {'pass' if r['ok'] else '**FAIL**'} |")
    print(f"\n{len(built)} of {len(a.chroms)} chromosome(s) built; "
          f"{sum(1 for r in rows if r['ok'])} passed every check")
    if problems:
        print("\nPROBLEMS — do not quote a number from these:")
        for p in problems:
            print(f"  {p}")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
