"""Where does PanMixer's capacity frontier cross OUR target_fidelity?

⚠ Standing constraint 5: privacy is comparable only at MATCHED UTILITY. If a crossing
exists, that capacity is where the privacy comparison is legitimate. docs/plan.md has
carried BOTH "a matched-fidelity capacity exists at ~0.17" and "no matched-fidelity
capacity exists" -- both measured against different releases of ours. This settles it
against the live chr21 release with 5 draws on their side.
"""
import json, glob, math, statistics as st
from collections import defaultdict

CAPS = [0.01, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 0.9]
draws = defaultdict(dict)
for pat in ("/home/ds85/pm_ourbuild/obf_out/*/cap0p*/chr21/stats.json",
            "/home/ds85/pm_ourbuild/obf_draws/*/cap0p*/s*/chr21/stats.json"):
    for f in glob.glob(pat):
        d = json.load(open(f))
        draws[(d["subject_name"], round(d["capacity"], 2))][d["seed"]] = d["changed_frac"]
th = {k: st.mean(list(v.values())) for k, v in draws.items()}

ep = json.load(open("logs/B_arm1_chr21/endpoint.json"))
ours = {r["target"]: {l: c["target_fidelity"] for l, c in r["cells"].items()}
        for r in ep["results"]}
subs = sorted(ours)

def cross(subj, our_val):
    """The capacity at which their 5-draw-mean frontier reaches our fidelity,
    by linear interpolation between the bracketing capacities."""
    xs = [(c, th[(subj, c)]) for c in CAPS]
    for (c0, v0), (c1, v1) in zip(xs, xs[1:]):
        if (v0 - our_val) * (v1 - our_val) <= 0 and v1 != v0:
            return c0 + (our_val - v0) * (c1 - c0) / (v1 - v0)
    return None

for LAB in ("0.50", "0.65", "0.80"):
    xs = [cross(s, ours[s][LAB]) for s in subs]
    ok = [x for x in xs if x is not None]
    print(f"\n  OUR acc_ceiling {LAB}  (our median fidelity "
          f"{st.median([ours[s][LAB] for s in subs]):.6f})")
    print(f"    subjects with a crossing: {len(ok)}/{len(subs)}")
    if ok:
        print(f"    matched-fidelity capacity: median {st.median(ok):.4f}  "
              f"min {min(ok):.4f}  max {max(ok):.4f}")
    nope = [s for s, x in zip(subs, xs) if x is None]
    if nope:
        print(f"    NO crossing for {len(nope)}: e.g. {nope[:4]}")
        s = nope[0]
        print(f"      {s}: ours {ours[s][LAB]:.6f}, their range "
              f"{th[(s,0.01)]:.6f}..{th[(s,0.9)]:.6f}")

print(f"\n  ⚠ A crossing is where their UTILITY equals ours. It is the capacity at which "
      f"a PRIVACY comparison becomes legitimate under constraint 5 -- it is NOT itself a "
      f"privacy result, and nothing about who is more private follows from it.")
print(f"  ⚠ Their frontier SATURATES from capacity ~0.3, so above the plateau the map "
      f"capacity -> fidelity is not invertible and no crossing can be interpolated there.")
