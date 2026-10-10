"""Run the mechanism on chr21 and emit the results table. Operational, not mechanism.

Metrics are ROWS, `tau` values are COLUMNS -- never the transpose, because a wide
table of arms invites cherry-picking which metrics to show. The table is EMITTED
by this script rather than hand-assembled, and every row that could not be
computed appears as `n/a` with a reason rather than being omitted.

⚠ WHAT THIS IS AND IS NOT. It is a diagnostic run at real scale: does the sampler
work at T = 305,886 and K = 88, what does the utility curve look like, and do the
runtime assertions hold on real data. It is NOT a pre-registered result. Standing
constraint 6 freezes `tau`, the utility and the endpoints before any attack
outcome is inspected, and `A-LCK-preregister` has not run -- nor have the four
deferred decisions R1-R4 been settled. EVERY NUMBER FROM HERE IS A v1 NUMBER AND
MUST SAY SO.

    python scripts/run_mechanism_chr21.py --target HG00096 --seeds 3
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

from mechanism import calibration as C   # noqa: E402
from mechanism import ffbs as F          # noqa: E402
from mechanism import hmm as H           # noqa: E402
from mechanism import utility as U       # noqa: E402

TAU_GRID = (0.0, 0.05, 0.10, 0.25, 0.50, 0.90)


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--target", default="HG00096")
    p.add_argument("--chrom", default="chr21")
    p.add_argument("--root", type=Path, default=REPO / "data")
    p.add_argument("--taus", type=float, nargs="+", default=list(TAU_GRID))
    p.add_argument("--seeds", type=int, default=3,
                   help="independent draws per tau. ⚠ These are an EVALUATION "
                        "device: a real release is ONE draw per genome, because "
                        "independent re-releases compose (constraint 5).")
    p.add_argument("--beta-scope", default="chain",
                   choices=list(U.BETA_SCOPES),
                   help="deferred decision R1; stamped on every number")
    p.add_argument("--out", type=Path, default=None)
    a = p.parse_args(argv)

    cohort = a.root / "cohort" / a.chrom
    tdir = a.root / "targets" / a.target / a.chrom

    hmm = H.from_cohort_dir(cohort)
    in_chain = np.load(cohort / "in_chain.npy")
    support = np.load(cohort / "support.npy")
    hap = np.load(cohort / "haplotypes.npy")[:, in_chain]
    readable = (np.load(tdir / "readable.npy")
                if (tdir / "readable.npy").exists() else None)
    util = U.beta_v(support, in_chain=in_chain, readable=readable,
                    scope=a.beta_scope)
    target = np.load(tdir / "path.npy")[in_chain]
    if target.ndim != 2 or target.shape[1] != 2:
        raise ValueError(f"path.npy is {target.shape}; expected (n_sites, 2)")

    # Target-independent, so computed ONCE and reused across every tau and seed.
    t0 = time.time()
    baseline = F.baseline_log_norm(hmm)
    t_base = time.time() - t0

    graph_id = f"{a.chrom}/hprc-v1.0-pggb/K{hmm.K}/T{hmm.T}"
    rows: dict[str, dict] = {}
    releases: list[dict] = []

    for tau in a.taus:
        eta = C.eta_tau(tau)
        draws = []
        for s in range(a.seeds):
            t0 = time.time()
            r = F.release(hmm, hap, target, util, tau, seed=1000 + s,
                          target_id=a.target, graph_id=graph_id,
                          baseline_norm=baseline)
            rec = r.record()
            rec["seconds"] = round(time.time() - t0, 2)
            releases.append(rec)
            draws.append(r)

        key = f"tau={tau:g}"
        u = np.array([d.u_achieved for d in draws])
        sw = np.array([int((np.diff(d.state_path, axis=0) != 0).sum()) for d in draws])
        lz = np.array([d.log_z_p for d in draws])
        rows[key] = {
            "utility_retained": f"{u.mean():.4f}",
            "u_ceiling": f"{draws[0].u_ceiling:.4f}",
            "utility_vs_ceiling": f"{u.mean()/draws[0].u_ceiling:.1%}",
            "p_succ_bound": f"{C.p_succ_bound(tau):.4f}",
            "log_z_p": f"{lz.mean():.4f}",
            "log_z_p_bracket": f"[0, {eta:.4f}]",
            "log_z_p_frac_of_eta": "0" if eta == 0 else f"{lz.mean()/eta:.1%}",
            "donor_switches": f"{sw.mean():.0f}",
            "eta_tau": f"{eta:.6f}",
            "epsilon_tau": f"{C.epsilon_tau(tau):.6f}",
            "seconds per release, diploid": f"{np.mean([r['seconds'] for r in releases[-a.seeds:]]):.1f}",
        }

    # Rows that CANNOT be computed here. Named, with the reason, never omitted.
    unavailable = {
        "attacker_accuracy":
            "n/a -- no attacker implemented. Pre-registration (A-LCK-preregister) "
            "must fix the attacker set and the attack database first; and a weak "
            "attacker would prove nothing anyway (constraint 7).",
        "two_target_indistinguishability":
            "n/a at this scale -- plug-in TV over a support of this size is biased "
            "UPWARD and cannot falsify anything. Exact at K=4, T=7 in "
            "tests/test_toy_enumeration.py; the worst pair there reached 43.1% of "
            "tau at tau=0.9.",
        "tv_empirical":
            "n/a -- same reason. DIAGNOSTIC only, never a headline number.",
        "target_fidelity":
            "n/a -- needs its ceiling computed on this graph, and the ceiling is "
            "target-side: about 9.9% of a target's chr21 non-reference calls have "
            "no record in G at all, carrying 20.5% of its -log f information. "
            "u_path is fidelity-on-what-the-graph-can-see and is NOT a substitute.",
        "af_loss / ld_loss":
            "n/a BY DEFINITION -- cohort-distortion metrics, undefined for an "
            "external-target release that never edits the cohort.",
    }

    order = ["utility_retained", "u_ceiling", "utility_vs_ceiling",
             "attacker_accuracy", "p_succ_bound", "log_z_p", "log_z_p_bracket",
             "log_z_p_frac_of_eta", "two_target_indistinguishability",
             "target_fidelity", "donor_switches", "eta_tau", "epsilon_tau",
             "tv_empirical", "af_loss / ld_loss", "seconds per release, diploid"]
    gates = {"utility_retained", "attacker_accuracy", "p_succ_bound", "log_z_p",
             "two_target_indistinguishability", "target_fidelity"}

    cols = list(rows)
    out = [
        f"### v1 diagnostic run — {a.target}, {a.chrom} — NOT a pre-registered result",
        "",
        "| metric | " + " | ".join(cols) + " |",
        "|---" * (len(cols) + 1) + "|",
    ]
    for name in order:
        # Backtick only real `docs/terms.md` identifiers. A prose row label in
        # backticks would claim to be a defined term and fail the docs contract,
        # which checks that every backticked label in a table resolves to an entry.
        label = (f"`{name}`" if name.replace("_", "").isalnum() else name)
        label += " *" if name in gates else ""
        if name in unavailable:
            out.append(f"| {label} | " + unavailable[name] +
                       " |" + "".join(" — |" for _ in cols[1:]))
        else:
            out.append(f"| {label} | " +
                       " | ".join(rows[c].get(name, "—") for c in cols) + " |")

    prov = (f"**Provenance.** {graph_id} · cohort 44 individuals / {hmm.K} haplotypes · "
            f"utility `{F.utility_id(util)}` (beta over {util.n_scope}/{util.n_chain} "
            f"chain positions) · {a.seeds} draws per tau · seeds 1000..{999+a.seeds} · "
            f"baseline log_norm {baseline:.4f} in {t_base:.1f}s · "
            f"chain_span PROVISIONAL · deferred decisions R1-R4 UNSETTLED")
    out += ["", prov, "",
            "⚠ The floor is the `tau=0` column: an untilted draw from `baseline_model`, "
            "which is a real cohort haplotype mosaic carrying ZERO target information. "
            "The ceiling is `u_ceiling`, which is BELOW 1 under "
            f"`beta={a.beta_scope}` because the target has no call at some weighted "
            "positions (deferred decision R1)."]

    text = "\n".join(out)
    print(text)
    if a.out:
        a.out.parent.mkdir(parents=True, exist_ok=True)
        a.out.write_text(text + "\n")
        (a.out.parent / (a.out.stem + ".releases.json")).write_text(
            json.dumps(releases, indent=2) + "\n")
        print(f"\nwrote {a.out} and {a.out.stem}.releases.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
