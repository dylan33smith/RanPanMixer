#!/usr/bin/env python
"""ARM 1's releases: the 44 cohort members, diploid, 22 autosomes, K = 86, leave-one-out.

⚠⚠ ARM 1, AND EVERY NUMBER HERE IS CONDITIONAL. The target IS a cohort member, so the
guarantee holds ONLY under leave-one-out — both of its haplotypes leave the donor panel
(`docs/terms.md:loo_cohort`, decision E). And `G` itself was built from all HPRC
assemblies including this target's, so arm 1's guarantee is conditional ON THE GRAPH.
State that; never claim it away. `src/releases.py:_check_arm` refuses an arm 1 release
with no `loo_id`, which is the machine-checked half of this.

⚠⚠ TWO CAVEATS TRAVEL WITH THE OUTPUT AND MUST REACH THE WRITE-UP.
 1. The `tau` grid is a declared EVALUATION DEVICE, never a release procedure
    (`prereg_2026-10-10.md:345`). Seven releases of one genome COMPOSE; we report the
    grid to characterise the mechanism, not as something a deployment would do.
 2. Arm 1's 44 releases have donor panels that are each other's TRUE HAPLOTYPES. The
    PER-RELEASE bound holds; the JOINT bound over the 44 does NOT (standing constraint 5).
    Report it; do not quietly design around it.

⚠ THE CHAIN IS PER-TARGET AND IS NOT A SUBSET OF ARM 2'S. Measured
(`logs/B_loo_chain_survey.md`): `in_chain` T spans 22,129,531 to 22,136,532 against arm
2's 22,134,562, and the `untraversable_block` set is NOT monotone in the donor panel — on
chr16 the block arm 2 excludes becomes traversable once HG02630 is held out. So every
chain, every `T` and all 22 `beta_share` values are recomputed here per target, from the
cohort arrays, rather than read from a file that could go stale.

⚠ ONE STREAM PER CHROMOSOME. `release()` is called with `stream_key=chrom`, because
deriving the stream from `seed` alone gave every chromosome of a genome release the SAME
uniforms, and `_pick` samples by inverse CDF — so every chromosome would have chosen its
donor at the same posterior quantile. The RECORDED seed stays the genome's, so each part
replays from `(seed, chrom)`. See `docs/bugs.md`.

    python scripts/arm1_release.py --workers 24
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
import time
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

import numpy as np                                        # noqa: E402

import releases as R                                      # noqa: E402
import store                                              # noqa: E402
from cohort import loo                                     # noqa: E402
from mechanism import ffbs, hmm, utility                   # noqa: E402

# ✅ FROZEN by prereg_2026-10-10.md section 1 under standing constraint 6, with each
# point's declared role. Not a tunable: all seven run.
PINNED_TAUS = (0.0, 0.05, 0.10, 0.25, 0.50, 0.90, 0.99)
SEED_BASE = 20261012
AUTOSOMES = loo.AUTOSOMES


def genome_seed(target: str, tau_index: int) -> int:
    """A distinct, reproducible seed per (target, `tau`) cell.

    ⚠ hashlib, NEVER `hash()`: Python's is PYTHONHASHSEED-salted and already made a
    decoy choice irreproducible once in this project.
    ⚠ And keyed by BOTH target and `tau` index, because seeding from one of them alone
    is how common random numbers got into this project three times.
    """
    h = int.from_bytes(hashlib.sha256(target.encode()).digest()[:4], "big")
    ss = np.random.SeedSequence([SEED_BASE, h, int(tau_index)])
    return int(ss.generate_state(1, dtype=np.uint32)[0])


def _axis_digest(chains: dict) -> str:
    """One digest for the genome's axis, and it is PER TARGET.

    Arm 1's chain depends on who is held out, so the axis a release is indexed by does
    too. Folding each chromosome's site-axis digest together with that target's own
    `in_chain` means a release can never be joined to a chain built for someone else.
    """
    h = hashlib.sha256()
    for c in AUTOSOMES:
        man = json.loads((REPO / "data" / "cohort" / c / "manifest.json").read_text())
        seen = man.get("axis_digests_seen") or []
        if len(seen) != 1:
            raise AssertionError(
                f"{c}: {len(seen)} site-axis digests in the manifest, expected exactly "
                f"1. A release cannot be stamped with an ambiguous axis.")
        h.update(c.encode())
        h.update(seen[0].encode())
        h.update(hashlib.sha256(
            np.ascontiguousarray(chains[c].in_chain).tobytes()).digest())
    return h.hexdigest()


def _inputs() -> dict[str, str]:
    d = {"git_rev": store.git_rev(REPO),
         "prereg": store.sha256_file(REPO / "paper" / "prereg_2026-10-10.md")}
    for f in sorted((REPO / "paper").glob("addendum_*.md")):
        d[f"addendum:{f.stem}"] = store.sha256_file(f)
    for c in AUTOSOMES:
        d[f"manifest:{c}"] = store.sha256_file(
            REPO / "data" / "cohort" / c / "manifest.json")
    return d


def run_target(target: str, taus, root: str | None) -> dict:
    """All seven `tau` for one held-out individual. Chromosome-outer, so each
    chromosome's arrays and its `baseline_log_norm` are computed ONCE and reused
    across the grid -- `baseline_log_norm` is target-independent but the HMM is not,
    since K = 86 and the chain are both per-target."""
    t0 = time.time()
    todo = [(i, t) for i, t in enumerate(taus)
            if not R.existing_seeds(target, t, "arm1", root)]
    if not todo:
        return {"target": target, "skipped": True, "cells": len(taus),
                "seconds": 0.0}

    chains = {c: loo.loo_chain(REPO / "data" / "cohort" / c, target)
              for c in AUTOSOMES}
    shares = loo.loo_shares({c: chains[c].t_chain for c in AUTOSOMES})
    axis = _axis_digest(chains)
    inputs = _inputs()

    parts: dict[int, list] = {i: [] for i, _ in todo}
    # ⚠⚠ THE EXACT ENDPOINT RIDES ALONG. `E[differ]` and `E_Q[u]` need the state
    # marginals, which are a backward pass over the SAME forward messages the draw
    # already computes. Measured 2026-10-12: computing them in a separate run doubles
    # `psi` and `forward` and cost ~6.5 h extra across arm 1's 308 cells. Pinned inert
    # with respect to the draw by tests/test_exact_out.py.
    ex: dict[int, dict] = {i: {"e_differ": 0.0, "e_q_u": 0.0, "called": 0,
                               "u_ceiling": 0.0} for i, _ in todo}
    for n_c, c in enumerate(AUTOSOMES, 1):
        d = REPO / "data" / "cohort" / c
        lc = chains[c]
        ic = lc.in_chain
        hap_all = np.asarray(np.load(d / "haplotypes.npy", mmap_mode="r"))
        panel = loo.masked_haplotypes(d, target, haplotypes=hap_all)[:, ic]
        rows = loo.sample_index(
            np.load(d / "haplotype_ids.npy", allow_pickle=True))[target]
        # ARM 1's target path IS its own two haplotype rows -- it is a cohort member,
        # so no target VCF exists or is needed.
        tgt = hap_all[list(rows)][:, ic].T
        m = hmm.from_arrays(panel, np.load(d / "genetic_pos.npy")[ic],
                            np.load(d / "positions.npy")[ic],
                            cohort_id=f"{c}_loo_{target}")
        util = utility.beta_v(lc.support, in_chain=ic, share=shares[c], chrom=c)
        base = ffbs.baseline_log_norm(m)
        # ⚠ `u_ceiling` is `tau`-INDEPENDENT, so it is accumulated once per chromosome.
        # Adding it inside the `tau` loop would multiply it by len(taus) and still print
        # a plausible table -- the bug caught in arm1_exact_endpoint.py before it ran.
        ceil_c = utility.diploid_utility(
            utility.u_ceiling(tgt[:, 0], util.beta),
            utility.u_ceiling(tgt[:, 1], util.beta))
        for i, tau in todo:
            box: dict = {}
            parts[i].append(ffbs.release(
                m, panel, tgt, util, tau, seed=genome_seed(target, i),
                target_id=target, loo_id=target, baseline_norm=base,
                stream_key=c, exact_out=box))
            assert len(box["u_haps"]) == 2, (
                f"{target}/{c}/tau={tau:g}: {len(box['u_haps'])} haplotype "
                f"contribution(s), expected 2 -- a diploid release must give two")
            ex[i]["e_differ"] += box["e_differ"]
            # the diploid utility is the MEAN of the two haplotypes, matching the
            # density exp(eta * u_mean) the sampler draws from
            ex[i]["e_q_u"] += 0.5 * (box["u_haps"][0] + box["u_haps"][1])
            ex[i]["called"] += box["called"]
            ex[i]["u_ceiling"] += ceil_c
        print(f"      {target} [{n_c}/22] {c} done ({time.time() - t0:.0f}s)",
              flush=True)
        del panel, hap_all, m, util

    written, cells = [], {}
    for i, tau in todo:
        gr = ffbs.GenomeRelease(parts=tuple(parts[i]))
        p = R.save_genome_release(gr, arm="arm1", axis_digest=axis, inputs=inputs,
                                  root=root, scope="genome")
        written.append(str(p))
        e = ex[i]
        cells[str(tau)] = {
            "e_differ": e["e_differ"], "called_cells": e["called"],
            "target_fidelity": e["e_differ"] / e["called"] if e["called"] else float("nan"),
            "e_q_u": e["e_q_u"], "u_ceiling": e["u_ceiling"],
            "e_q_u_over_ceiling": e["e_q_u"] / e["u_ceiling"] if e["u_ceiling"] else float("nan"),
        }
    return {"target": target, "skipped": False, "cells": len(todo),
            "written": written, "seconds": time.time() - t0,
            "k_states": chains[AUTOSOMES[0]].k_states,
            "t_total": sum(chains[c].t_chain for c in AUTOSOMES),
            "exact": cells}


def _one(args):
    target, taus, root = args
    try:
        return run_target(target, taus, root)
    except Exception as e:                      # a worker must report, never vanish
        return {"target": target, "error": f"{type(e).__name__}: {e}",
                "traceback": traceback.format_exc()[-1500:]}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--targets", nargs="+", default=None)
    ap.add_argument("--taus", type=float, nargs="+", default=list(PINNED_TAUS))
    ap.add_argument("--workers", type=int, default=24)
    ap.add_argument("--root", default=None)
    ap.add_argument("--json", type=Path, default=None)
    ap.add_argument("--endpoint-json", type=Path, default=None,
                    help="write the exact endpoint in the shape "
                         "scripts/arm1_results.py --endpoint expects")
    a = ap.parse_args(argv)

    if set(a.taus) != set(PINNED_TAUS):
        print(f"⚠⚠ RUNNING A NON-PINNED `tau` SET {a.taus}. The pre-registration freezes "
              f"{list(PINNED_TAUS)} under standing constraint 6. This is legitimate only "
              f"as a PILOT and its output must not be reported as the arm 1 result.\n")

    ids = np.load(REPO / "data" / "cohort" / "chr1" / "haplotype_ids.npy",
                  allow_pickle=True)
    members = sorted(loo.sample_index(ids))
    targets = a.targets or members
    unknown = [t for t in targets if t not in members]
    if unknown:
        raise SystemExit(f"not cohort members: {unknown}")

    n_cells = len(targets) * len(a.taus)
    print(f"ARM 1 releases. {len(targets)} held-out individuals x {len(a.taus)} `tau` "
          f"= {n_cells} genome-wide releases, K = 86, diploid, 22 autosomes.")
    print(f"{a.workers} workers, one held-out individual each; chromosome-outer so each "
          f"`baseline_log_norm` is computed once per (target, chromosome).")
    print(f"⚠ The `tau` grid is an EVALUATION DEVICE, not a release procedure, and the "
          f"JOINT bound over the 44 does not hold — only the per-release one.\n")

    t0 = time.time()
    out = []
    with ProcessPoolExecutor(max_workers=a.workers) as ex:
        futs = {ex.submit(_one, (t, a.taus, a.root)): t for t in targets}
        for n, f in enumerate(as_completed(futs), 1):
            r = f.result()
            out.append(r)
            if "error" in r:
                print(f"  [{n}/{len(targets)}] ❌ {r['target']}: {r['error']}",
                      flush=True)
            elif r["skipped"]:
                print(f"  [{n}/{len(targets)}] ⏭ {r['target']}: already written",
                      flush=True)
            else:
                print(f"  [{n}/{len(targets)}] ✅ {r['target']}: {r['cells']} release(s), "
                      f"T = {r['t_total']:,} ({r['seconds'] / 60:.1f} min) "
                      f"[{(time.time() - t0) / 60:.0f} min elapsed]", flush=True)

    errs = [r for r in out if "error" in r]
    # ⚠ ASSERT COMPLETENESS FROM THE DISK, not from this loop's bookkeeping: a fan-out
    # that reports its own success cannot tell "wrote nothing" from "ran nothing".
    missing = [(t, tau) for t in targets for tau in a.taus
               if not R.existing_seeds(t, tau, "arm1", a.root)]
    print(f"\n{len(out) - len(errs)}/{len(targets)} targets completed, "
          f"{len(errs)} failed, {len(missing)} of {n_cells} cells missing on disk.")
    if errs:
        for r in errs[:5]:
            print(f"\n❌ {r['target']}\n{r['traceback']}")
    if a.json:
        a.json.write_text(json.dumps(
            {"targets": targets, "taus": a.taus, "results": out,
             "missing": [[t, tau] for t, tau in missing]}, indent=1) + "\n")
        print(f"JSON -> {a.json}")
    if a.endpoint_json:
        # ⚠ Only targets that actually produced exact values go in. A skipped target
        # (its cells already written by an earlier run) has none, and inventing an
        # entry for it would make a resumed run look complete when its endpoint is not.
        res = [{"target": r["target"], "cells": r["exact"],
                "t_total": r.get("t_total"), "k_states": r.get("k_states")}
               for r in out if "exact" in r and r.get("exact")]
        a.endpoint_json.write_text(json.dumps(
            {"arm": "arm1", "taus": a.taus, "k_states": 86, "results": res},
            indent=1) + "\n")
        skipped = [r["target"] for r in out if r.get("skipped")]
        print(f"ENDPOINT JSON -> {a.endpoint_json} ({len(res)} target(s))")
        if skipped:
            print(f"⚠ {len(skipped)} target(s) were SKIPPED as already-released and so "
                  f"carry NO exact endpoint here: {skipped[:5]}. Their releases exist; "
                  f"their endpoint must come from scripts/arm1_exact_endpoint.py or "
                  f"from a re-run with the releases removed.")
    if errs or missing:
        raise SystemExit(
            f"INCOMPLETE: {len(errs)} worker error(s), {len(missing)} cell(s) missing. "
            f"Arm 1 is not releasable until every cell exists; re-run to resume.")
    print(f"\n✅ all {n_cells} arm 1 releases written in "
          f"{(time.time() - t0) / 60:.1f} min.")
    return 0


if __name__ == "__main__":
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    sys.exit(main())
