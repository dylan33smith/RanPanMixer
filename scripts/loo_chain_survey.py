#!/usr/bin/env python
"""ARM 1's 44 per-target chains, measured: how far does leave-one-out move the chain?

⚠⚠ WHY THIS RUNS BEFORE ANY ARM 1 RELEASE. Decision K created a risk recorded in
`docs/plan.md` before chunk 8 opened: removing both of the target's haplotypes shrinks
every donor-intersection, so a zero-genetic-distance block that some donor spans at
K = 88 may be spanned by NONE at K = 86 — and such a block DEAD-ENDS the forward pass.
That makes arm 1's chain per-target: 44 chains, 44 `T` values, 44 `beta_share`
allocations. This script measures all 968 of them so the release stage is configured from
MEASURED chain lengths rather than from arm 2's, which would reintroduce exactly the dead
end decision K removed.

⚠ IT ALSO SETTLES A CONTAINMENT QUESTION RATHER THAN ASSUMING IT. `support_loo <=
support` is provable and `src/cohort/loo.py` asserts it. Containment of the CHAIN is not
provable: dropping positions SHORTENS the surviving zero-distance runs, and a shorter run
needs a donor to span fewer positions, so a block arm 2 found untraversable can become
traversable under leave-one-out. Reported per target, never assumed.

ROUTE AND WHY IT IS FAST. Chromosome-outer, target-inner, so each `haplotypes.npy` is
read ONCE rather than 44 times (chr1 alone is ~158 MB). Within a chromosome the two
held-out rows are ZEROED in the called-matrix rather than deleted, which is EXACTLY
equivalent for this computation: `untraversable_blocks` asks
`allowed[:, a:b+1].all(axis=1).any()`, and an all-False row can never be all-True, so it
contributes nothing either way. ⚠ That equivalence is spot-checked against the
authoritative `loo.loo_chain` path rather than argued — see `--validate`.

    python scripts/loo_chain_survey.py --json logs/B_loo_chain_survey.json
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

from cohort import genmap, loo                            # noqa: E402

AUTOSOMES = loo.AUTOSOMES


def survey_chrom(chrom: str, members: list[str], validate: int = 0) -> dict:
    d = REPO / "data" / "cohort" / chrom
    positions = np.load(d / "positions.npy")
    cm = np.load(d / "genetic_pos.npy")
    hap = np.asarray(np.load(d / "haplotypes.npy", mmap_mode="r"))
    ids = np.load(d / "haplotype_ids.npy", allow_pickle=True)
    index = loo.sample_index(ids)
    arm2 = np.load(d / "in_chain.npy").astype(bool)

    called = (hap != -1)                                  # (K, n) bool, read once
    tot = called.sum(axis=0).astype(np.int64)

    # ⚠ THE ARM 2 REFERENCE IS RECOMPUTED FROM THIS CHROMOSOME, not quoted from the
    # genome-wide record. Printing a genome-wide 305 / 20 / 172,142 beside a
    # single-chromosome survey is a mislabelled denominator, which is the error class
    # this project keeps hitting. Validated separately: at K = 88 this composition
    # reproduces `in_chain.npy` bit-identically.
    a2_mask = genmap.chain_mask(positions, cm, None, None, tot)
    a2_mask, a2_trav = genmap.drop_untraversable(a2_mask, cm, called, positions)
    assert np.array_equal(a2_mask, arm2), (
        f"{chrom}: the K = 88 reconstruction does not match in_chain.npy, so the arm 2 "
        f"reference column cannot be trusted")
    ref = {"t": int(arm2.sum()), "support_zero": int((tot == 0).sum()),
           "blocks": int(a2_trav["untraversable_blocks_dropped"]),
           "dropped": int(a2_trav["positions_dropped"])}

    out = {"__arm2__": ref}
    for n, m in enumerate(members):
        rows = list(index[m])
        keep_rows = called[rows].copy()                   # save, zero, restore
        called[rows] = False
        support = tot - keep_rows.sum(axis=0).astype(np.int64)
        mask = genmap.chain_mask(positions, cm, None, None, support)
        mask, trav = genmap.drop_untraversable(mask, cm, called, positions)
        left = genmap.untraversable_blocks(cm[mask], called[:, mask])
        assert not left, f"{chrom}/{m}: {len(left)} untraversable block(s) remain"

        if n < validate:
            # the authoritative path, which DELETES the rows instead of zeroing them
            ref = loo.loo_chain(d, m, haplotypes=hap)
            assert np.array_equal(ref.in_chain, mask), \
                f"{chrom}/{m}: the fast path disagrees with loo.loo_chain"
            assert ref.k_states == called.shape[0] - 2

        out[m] = {
            "t_chain": int(mask.sum()),
            "blocks_dropped": int(trav["untraversable_blocks_dropped"]),
            "positions_dropped": int(trav["positions_dropped"]),
            "support_zero": int((support == 0).sum()),
            "outside_arm2": int((mask & ~arm2).sum()),
            "t_arm2": int(arm2.sum()),
        }
        called[rows] = keep_rows
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--chroms", nargs="+", default=list(AUTOSOMES))
    ap.add_argument("--validate", type=int, default=2,
                    help="per chromosome, how many targets to cross-check against "
                         "loo.loo_chain (the row-deleting path)")
    ap.add_argument("--json", type=Path, default=None)
    ap.add_argument("--from-json", type=Path, default=None,
                    help="re-report a saved survey instead of recomputing 968 chains")
    a = ap.parse_args(argv)

    if a.from_json:
        dd = json.loads(a.from_json.read_text())
        print(f"re-reporting a saved survey from {a.from_json}: "
              f"{len(dd['members'])} targets x {len(dd['chroms'])} chromosomes. "
              f"NO chains recomputed.")
        _report(dd["per_chrom"], dd["members"], dd["chroms"])
        return 0

    ids = np.load(REPO / "data" / "cohort" / a.chroms[0] / "haplotype_ids.npy",
                  allow_pickle=True)
    members = sorted(loo.sample_index(ids))
    print(f"ARM 1 leave-one-out chains. {len(members)} cohort members x "
          f"{len(a.chroms)} chromosome(s) = {len(members) * len(a.chroms):,} chains, "
          f"K = {len(ids) - 2} each.")
    print(f"cross-checking {a.validate} target(s) per chromosome against the "
          f"row-deleting path.\n")

    per: dict[str, dict] = {}
    t0 = time.time()
    for c in a.chroms:
        per[c] = survey_chrom(c, members, validate=a.validate)
        d = per[c]
        vals = [v for k, v in d.items() if k != "__arm2__"]
        ta = d["__arm2__"]["t"]
        dt = [v["t_chain"] - ta for v in vals]
        print(f"  {c}: T_arm2 = {ta:,}, "
              f"LOO T spans {min(v['t_chain'] for v in vals):,} to "
              f"{max(v['t_chain'] for v in vals):,} "
              f"(delta {min(dt):+,} to {max(dt):+,}), "
              f"outside-arm2 max {max(v['outside_arm2'] for v in vals):,} "
              f"({time.time() - t0:.0f}s)", flush=True)

    # ⚠ assert completeness: a fan-out that filters cannot tell "found nothing" from
    # "ran nothing".
    assert set(per) == set(a.chroms), f"surveyed {sorted(per)}, asked for {sorted(a.chroms)}"
    for c in a.chroms:
        assert set(per[c]) - {"__arm2__"} == set(members), \
            f"{c}: {len(per[c]) - 1} of {len(members)} members"

    rows = _report(per, members, a.chroms)
    if a.json:
        a.json.write_text(json.dumps(
            {"members": members, "chroms": list(a.chroms), "per_chrom": per,
             "per_target": rows}, indent=1) + "\n")
        print(f"\nJSON -> {a.json}")
    return 0


def _report(per, members, chroms) -> list[dict]:
    full = len(chroms) == 22
    rows = []
    for m in members:
        t_by = {c: per[c][m]["t_chain"] for c in chroms}
        tot = sum(t_by.values())
        rec = {"target": m, "t_total": tot,
               "blocks_dropped": sum(per[c][m]["blocks_dropped"] for c in chroms),
               "positions_dropped": sum(per[c][m]["positions_dropped"] for c in chroms),
               "support_zero": sum(per[c][m]["support_zero"] for c in chroms),
               "outside_arm2": sum(per[c][m]["outside_arm2"] for c in chroms)}
        if full:
            sh = loo.loo_shares(t_by)
            rec["shares"] = sh
            rec["share_chr1"] = sh["chr1"]
        rows.append(rec)

    t_arm2 = sum(per[c]["__arm2__"]["t"] for c in chroms)
    print(f"\nARM 1 per-target chains. {len(members)} targets, "
          f"{'GENOME-WIDE (22 autosomes)' if full else f'{len(chroms)} chromosome(s)'}, "
          f"K = 86, exact counts.\n")
    hdr = ["metric", "min", "median", "max", "ARM 2 (K = 88)"]
    print("| " + " | ".join(hdr) + " |")
    print("|" + "---|" * len(hdr))

    def line(name, vals, ref, fmt=",.0f"):
        v = np.array(vals, dtype=float)
        print(f"| {name} | {v.min():{fmt}} | {np.median(v):{fmt}} | "
              f"{v.max():{fmt}} | {ref} |")

    a2 = {k: sum(per[c]["__arm2__"][k] for c in chroms)
          for k in ("t", "support_zero", "blocks", "dropped")}
    line("`in_chain` T", [r["t_total"] for r in rows], f"{a2['t']:,}")
    line("positions with `support` = 0", [r["support_zero"] for r in rows],
         f"{a2['support_zero']:,}")
    line("`untraversable_block` count", [r["blocks_dropped"] for r in rows],
         f"{a2['blocks']:,}")
    line("positions dropped by decision K", [r["positions_dropped"] for r in rows],
         f"{a2['dropped']:,}")
    line("positions NOT in arm 2's chain", [r["outside_arm2"] for r in rows], "0")
    if full:
        line("`beta_share` of chr1", [r["share_chr1"] for r in rows],
             f"{per['chr1']['__arm2__']['t'] / a2['t']:.9f}",
             fmt=".9f")

    dts = [r["t_total"] - t_arm2 for r in rows]
    print(f"\n⚠⚠ **EVERY TARGET GETS ITS OWN CHAIN.** T spans "
          f"{min(r['t_total'] for r in rows):,} to {max(r['t_total'] for r in rows):,} "
          f"across the {len(members)} targets, i.e. {min(dts):+,} to {max(dts):+,} "
          f"against arm 2's {t_arm2:,}. Arm 1 therefore carries {len(members)} values of "
          f"`T` and {len(members)} `beta_share` allocations; reusing arm 2's would make "
          f"`sum beta_v` != 1 over the scope `tau` covers.")
    out = sum(r["outside_arm2"] for r in rows)
    if out:
        print(f"⚠⚠ CONTAINMENT FAILS: {out:,} (target, position) pairs are in a "
              f"leave-one-out chain but NOT in arm 2's. This is the predicted effect — "
              f"dropping positions shortens the surviving zero-distance runs, so a block "
              f"arm 2 found untraversable can become traversable at K = 86. Arm 1 must "
              f"NOT be built by subsetting arm 2's chain.")
    else:
        print(f"✅ Containment holds on this data: every leave-one-out chain is a subset "
              f"of arm 2's. ⚠ NOT a theorem — it is not provable and must be re-measured "
              f"if the build changes; see the module docstring.")
    print(f"⚠ These are ARM 1 numbers at K = 86. Arm 2 is K = 88 and holds nothing out; "
          f"never read one against the other without saying which (CLAUDE.md).")
    return rows


if __name__ == "__main__":
    sys.exit(main())
