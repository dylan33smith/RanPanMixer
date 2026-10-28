import json, glob, math, statistics as st
from collections import defaultdict

# seed 123 lives in obf_out/, seeds 124-127 in obf_draws/
draws = defaultdict(dict)          # (subject, cap) -> {seed: changed_frac}
for f in glob.glob("/home/ds85/pm_ourbuild/obf_out/*/cap0p*/chr21/stats.json"):
    d = json.load(open(f)); draws[(d["subject_name"], round(d["capacity"],2))][d["seed"]] = d["changed_frac"]
for f in glob.glob("/home/ds85/pm_ourbuild/obf_draws/*/cap0p*/s*/chr21/stats.json"):
    d = json.load(open(f)); draws[(d["subject_name"], round(d["capacity"],2))][d["seed"]] = d["changed_frac"]

cells = sorted(draws)
ndraw = {len(v) for v in draws.values()}
print(f"  cells: {len(cells)} (expect 352)   draws per cell: {sorted(ndraw)}")
assert len(cells) == 352 and ndraw == {5}, "incomplete draw matrix"

CAPS = [0.01,0.05,0.1,0.2,0.3,0.5,0.7,0.9]
print(f"\n**THEIR DRAW NOISE, 5 seeds per cell, n = 44 subjects per capacity.**\n")
print("| metric | " + " | ".join(f"{c:g}" for c in CAPS) + " |")
print("|" + "---|"*(len(CAPS)+1))
for lab, fn in (("within-cell sd, median", st.median), ("within-cell sd, max", max)):
    print(f"| {lab} | " + " | ".join(
        f"{fn([st.stdev(list(draws[(s,c)].values())) for s in {k[0] for k in cells}]):.2e}"
        for c in CAPS) + " |")
print(f"| between-subject sd | " + " | ".join(
    f"{st.stdev([st.mean(list(draws[(s,c)].values())) for s in {k[0] for k in cells}]):.2e}"
    for c in CAPS) + " |")

sd_w = st.median([st.stdev(list(draws[(s,0.9)].values())) for s in {k[0] for k in cells}])
sd_b = st.stdev([st.mean(list(draws[(s,0.9)].values())) for s in {k[0] for k in cells}])
print(f"\n  AT SATURATION: within-cell (draw) sd {sd_w:.3e} vs between-subject sd "
      f"{sd_b:.3e} -- draw noise is {100*sd_w/sd_b:.1f}% of the biological spread")
print(f"  SEM of a 5-draw mean: {sd_w/math.sqrt(5):.3e} (one draw: {sd_w:.3e})")

# ---- the paired test, now against their 5-draw MEAN ----
ep = json.load(open("logs/B_arm1_chr21/endpoint.json"))
ours = {r["target"]: r["cells"]["0.80"]["target_fidelity"] for r in ep["results"]}
subs = sorted(ours)
for label, pick in (("ONE DRAW (seed 123)", lambda s: draws[(s,0.9)][123]),
                    ("MEAN of 5 DRAWS",     lambda s: st.mean(list(draws[(s,0.9)].values())))):
    diffs = [ours[s] - pick(s) for s in subs]
    neg = sum(1 for x in diffs if x < 0); n = len(diffs)
    k = min(neg, n-neg)
    p = 2*sum(math.comb(n,i) for i in range(k+1))/2**n
    z = st.mean(diffs)/(st.stdev(diffs)/math.sqrt(n))
    print(f"\n  PAIRED vs {label}:")
    print(f"    more faithful in {neg}/{n}   median {st.median(diffs):+.6f}  "
          f"mean {st.mean(diffs):+.6f}  sd {st.stdev(diffs):.6f}")
    print(f"    sign test p = {p:.3e}      paired z = {z:+.2f}")
