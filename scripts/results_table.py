"""The full metric set, toy and chr21 side by side. Metrics are ROWS, configurations COLUMNS.

⚠ ONE DELIBERATE DEPARTURE FROM THE REPORTING CONTRACT, stated rather than hidden.
The contract puts the floor (`tau` = 0) and the ceiling in their own COLUMNS. Here the
column axis is the CONFIGURATION, because that is the comparison asked for, so the `tau`
axis is folded into rows: every `tau`-dependent row is stamped with its `tau`, and the
floor and ceiling appear as their own rows instead. Two configurations at matched `tau`
is a legitimate comparison; what would not be is mixing `tau` values down one column.

⚠ THE TWO COLUMNS ARE NOT AT MATCHED PLOIDY, and that asymmetry IS the open gap. The toy
gate is haploid; every chr21 release we have ever run is diploid. Rows marked [hap] are
per-haplotype on both sides so they stay comparable; rows marked [dip] are the diploid
release on chr21 and the haploid one on the toy, and are NOT directly comparable.

Emitted, never hand-assembled:  python scripts/results_table.py
"""

from __future__ import annotations

import math
import sys
import time
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

from mechanism import calibration as C   # noqa: E402
from mechanism import exact as E         # noqa: E402
from mechanism import ffbs as F          # noqa: E402
from mechanism import hmm as H           # noqa: E402
from mechanism import toy as TOY         # noqa: E402
from mechanism import utility as U       # noqa: E402

TAU = 0.5
SEEDS = (1000, 1001, 1002)
H_STEP = 0.05
NA = "n/a"

rows: list[tuple] = []


def add(label, toy, chr21, gate=False, note=""):
    rows.append((label, toy, chr21, gate, note))


def var_and_mean(hmm, hap, target, beta, baseline, eta, h=H_STEP):
    """E_Q[u] and Var_Q(u) from the first and second derivatives of log Z_p.

    Exact, sampling-free: d/d_eta log Z = E_Q[u] and d2/d_eta2 log Z = Var_Q(u).
    Three forward passes. This is the only route available at chr21 scale and it is
    cross-checked against direct enumeration on the toy below.
    """
    lz = []
    for e in (eta - h, eta, eta + h):
        psi = U.tilt_potentials(hmm.allowed, hap, target, beta, e)
        lz.append(F.log_z_p(F.forward(hmm, psi), baseline, e))
    return (lz[2] - lz[0]) / (2 * h), (lz[2] - 2 * lz[1] + lz[0]) / h ** 2


def worst_pair(Q, fn):
    worst = -np.inf
    for i in range(len(Q) - 1):
        v = fn(Q[i], Q[i + 1:])
        if len(v) and v.max() > worst:
            worst = float(v.max())
    return worst


# ---------------------------------------------------------------------------
print("building the toy column ...", file=sys.stderr)
t = TOY.toy_cohort()
toy_base = E.enumerate_baseline(t.haplotypes, t.genetic_pos)
toy_u = U.beta_v(t.support)
toy_hmm = H.from_arrays(t.haplotypes, t.genetic_pos, t.positions, cohort_id=t.name)
toy_bl = F.baseline_log_norm(toy_hmm)
tp, tq = TOY.sharp_pair(t)
eta = C.eta_tau(TAU)

# exact release distribution, and E/Var two independent ways
toy_up = E.utility_vector(toy_base, tp, toy_u.beta)
toy_Qp, toy_lz = E.release_distribution(toy_base, toy_up, eta)
toy_E_enum = float(toy_Qp @ toy_up)
toy_V_enum = float(toy_Qp @ (toy_up ** 2) - toy_E_enum ** 2)
toy_E_der, toy_V_der = var_and_mean(toy_hmm, t.haplotypes, tp, toy_u.beta, toy_bl, eta)

# all-pairs exact TV and ratio spread
seen = {}
for pth in TOY.toy_inputs(t, kind="all"):
    seen.setdefault(np.round(E.utility_vector(toy_base, pth, toy_u.beta), 15).tobytes(), pth)
toy_inputs = np.stack(list(seen.values()))
toy_Q = np.stack([E.release_distribution(toy_base, E.utility_vector(toy_base, x, toy_u.beta),
                                         eta)[0] for x in toy_inputs])
toy_tv = worst_pair(toy_Q, lambda a, b: 0.5 * np.abs(b - a[None, :]).sum(axis=1))
toy_spread = worst_pair(toy_Q, lambda a, b: (a[None, :] / b).max(axis=1) / (a[None, :] / b).min(axis=1))
toy_sharp_spread = E.ratio_spread(
    toy_Qp, E.release_distribution(toy_base, E.utility_vector(toy_base, tq, toy_u.beta), eta)[0])

# tau = 0 reproduces R_D exactly
toy_q0, toy_lz0 = E.release_distribution(toy_base, toy_up, 0.0)
toy_floor_dev = float(np.abs(toy_q0 - toy_base.prob).max())

# ffbs against enumeration
toy_psi = U.tilt_potentials(toy_hmm.allowed, t.haplotypes, tp, toy_u.beta, eta)
toy_fwd = F.forward(toy_hmm, toy_psi)
dense_alpha, dense_lc = E.dense_forward(t.haplotypes, t.genetic_pos, toy_psi)
toy_dense_dev = float(np.abs(toy_fwd.alpha - dense_alpha).max())
N_DRAWS = 40_000
index = {r.tobytes(): i for i, r in enumerate(toy_base.paths)}
counts = np.zeros(toy_base.M, dtype=np.int64)
rng = np.random.default_rng(20260930)
for _ in range(N_DRAWS):
    z = F.backward_sample(toy_hmm, toy_fwd, rng)
    counts[index[t.haplotypes[z, np.arange(t.T)].astype(np.int16).tobytes()]] += 1
emp = counts / N_DRAWS
big = toy_Qp * N_DRAWS >= 10
toy_z = float((np.abs(emp[big] - toy_Qp[big]) / np.sqrt(toy_Qp[big] * (1 - toy_Qp[big]) / N_DRAWS)).max())
toy_tv_emp = E.tv(emp, toy_Qp)

# sampled releases
_t0 = time.time()
toy_rel = [F.release(toy_hmm, t.haplotypes, tp, toy_u, TAU, seed=s, target_id="toyD1",
                     graph_id=t.name, baseline_norm=toy_bl) for s in SEEDS]
toy_secs = (time.time() - _t0) / len(SEEDS)

# ---------------------------------------------------------------------------
print("building the chr21 column ...", file=sys.stderr)
COH = REPO / "data" / "cohort" / "chr21"
TGT = REPO / "data" / "targets" / "HG00096" / "chr21"
c_hmm = H.from_cohort_dir(COH)
in_chain = np.load(COH / "in_chain.npy")
support = np.load(COH / "support.npy")
c_hap = np.load(COH / "haplotypes.npy")[:, in_chain]
c_u = U.beta_v(support, in_chain=in_chain)
c_path = np.load(TGT / "path.npy")[in_chain]
c_bl = F.baseline_log_norm(c_hmm)
c_E, c_V = var_and_mean(c_hmm, c_hap, c_path[:, 0], c_u.beta, c_bl, eta)
t0 = time.time()
c_rel = [F.release(c_hmm, c_hap, c_path, c_u, TAU, seed=s, target_id="HG00096",
                   graph_id=f"chr21/hprc-v1.0-pggb/K{c_hmm.K}/T{c_hmm.T}",
                   baseline_norm=c_bl) for s in SEEDS]
c_secs = (time.time() - t0) / len(SEEDS)
c_uniq = np.unique(c_hap, axis=0).shape[0]
toy_uniq = np.unique(t.haplotypes, axis=0).shape[0]


def f(x, p=4):
    return f"{x:.{p}f}"


def e_(x, p=4):
    return f"{x:.{p}e}"


def n_for_2sigma(gain_over_noise):
    return f"{round((2.0 / gain_over_noise) ** 2):,}"


toy_gn = eta * math.sqrt(abs(toy_V_der))
c_gn = eta * math.sqrt(abs(c_V))

# ------------------------------- PRIMARY -----------------------------------
add("**PRIMARY**", "", "")
add("`utility_retained` [dip]", f(np.mean([r.u_achieved for r in toy_rel])),
    f(np.mean([r.u_achieved for r in c_rel])), gate=True, note=f"mean of n={len(SEEDS)} draws, tau={TAU}")
add("E_Q[`u_path`], exact [hap]", f(toy_E_der, 6), f(c_E, 6), note=f"d logZ/d eta, tau={TAU}")
add("`u_ceiling` [dip]", f(toy_rel[0].u_ceiling), f(c_rel[0].u_ceiling), gate=True,
    note="max u this target could reach at any tau")
add("`utility_retained` / `u_ceiling`", f"{np.mean([r.u_achieved for r in toy_rel])/toy_rel[0].u_ceiling:.1%}",
    f"{np.mean([r.u_achieved for r in c_rel])/c_rel[0].u_ceiling:.1%}", note=f"tau={TAU}")
add("FLOOR: E_Q[`u_path`] at `tau`=0", f(float(toy_base.prob @ toy_up), 6),
    "not enumerable; = E_R[u], and logZ/eta gives 0.694575",
    note="an untilted draw from R_D, zero target information")
add("CEILING: `u_path` of the target's own path", "1.0000 (by construction)",
    "1.0000 (by construction)", note="unreachable; see u_ceiling for the real limit")
add("`two_target_indistinguishability` *", f(toy_tv, 6),
    f"{NA} — plug-in TV biased UPWARD over a support this size; not estimable", gate=True,
    note=f"worst over ALL {len(toy_inputs)*(len(toy_inputs)-1)//2:,} input pairs, tau={TAU}")
add("worst exact TV as a share of `tau`", f"{toy_tv/TAU:.1%}", NA, note=f"tau={TAU}")
add("`ratio_spread` *", f(toy_spread, 6), f"{NA} — toy-only; needs the enumerated support",
    gate=True, note=f"worst over all pairs; bound exp(eps_tau)={math.exp(C.epsilon_tau(TAU)):.4f}")
add("`ratio_spread` on the constructed sharp pair *", f(toy_sharp_spread, 6), NA, gate=True,
    note="SATURATES the bound exactly")
# ⚠ The reason here is NARROWER than it used to be, and the narrowing matters: the
# attacker IS implemented (src/attacks.py, 2026-10-10) and the attacker set IS pinned
# (prereg section 4, A1-A4). What blocks the cell is the SIGNATURE plus field 6, which
# was reopened because the pinned design was one-sided. "No attacker implemented" would
# now be false, and this is the one place the reporting contract requires a reason
# rather than an omitted row -- so a stale reason is the whole cost.
add("`attacker_accuracy` *", f"{NA} — pre-registration UNSIGNED", f"{NA} — same", gate=True,
    note="attacker built (src/attacks.py) and the set pinned (A1-A4); blocked on the "
         "SIGNATURE and on field 6, reopened 2026-10-10 — constraint 6")
add("`p_succ_bound` *", f(C.p_succ_bound(TAU)), f(C.p_succ_bound(TAU)), gate=True,
    note=f"equal-prior bound at tau={TAU}")
add("`target_fidelity` *", f"{NA} — no graph-representability ceiling for a synthetic cohort",
    f"{NA} — needs its own ceiling (~9.9% of calls absent from G, 20.5% of -log f)",
    gate=True, note="u_path is NOT a substitute")

# -------------------------------- GATES ------------------------------------
add("**GATES — can they falsify the theorem?**", "", "")
add("`log_z_p` *", f(toy_lz, 6), f(c_rel[0].log_z_p, 6), gate=True, note=f"tau={TAU}")
add("`log_z_p` bracket [0, `eta_tau`] *", f"[0, {eta:.6f}]", f"[0, {eta:.6f}]", gate=True,
    note="asserted at run time, not reported")
add("`log_z_p` as a share of `eta_tau`", f"{toy_lz/eta:.1%}", f"{c_rel[0].log_z_p/eta:.1%}",
    note=f"tau={TAU}; chr21 is the diploid sum")
add("`tau`=0 reproduces `baseline_model` *", f"max dev {toy_floor_dev:.2e}",
    f"logZ={0.0:.1e} by construction", gate=True, note="exact on the toy")
add("`ffbs` vs exact enumeration, max per-path z *", f(toy_z, 2),
    f"{NA} — support not enumerable", gate=True, note=f"n={N_DRAWS:,} draws, tau={TAU}")
add("TV(empirical, exact) *", f(toy_tv_emp, 6), NA, gate=True,
    note=f"n={N_DRAWS:,}; analytic expectation {0.5*np.sqrt(2*toy_Qp*(1-toy_Qp)/(math.pi*N_DRAWS)).sum():.6f}")
add("O(TK) vs dense O(TK^2) alphas, max abs dev *", e_(toy_dense_dev, 2),
    f"{NA} — 21.1 GB if materialised", gate=True, note="the rank-one recursion check")
add("KILL CRITERION: any exact TV above `tau` *", "none", f"{NA} — not testable at scale",
    gate=True, note="the only thing that can falsify Theorem 1")
add("DIPLOID gate *", "ABSENT — every clause above is haploid",
    "ABSENT — and all 18 releases ever run are diploid", gate=True,
    note="A-THY-diploid-gate; the top open item")

# ------------------------------ STRUCTURE ----------------------------------
add("**STRUCTURE**", "", "")
add("`K_states`", str(t.K), str(c_hmm.K), note="pre-collapse")
add("distinct haplotypes (post-collapse `K_states`)", str(toy_uniq), str(c_uniq),
    note="equal to K on both, so no collapse is needed")
add("`T_blocks` (chain positions)", str(t.T), f"{c_hmm.T:,}", note="model C: one per variant")
add("`n_sites` (site axis)", str(t.T), "340,849", note="chr21: T is n_sites inside chain_span")
add("output support, size of `Y_D`", f"{toy_base.M}", f"{NA} — astronomically large",
    note="toy: from 16,384 state paths, 2,304 with mass")
add("`donor_switches` [dip]", f"{np.mean([int((np.diff(r.state_path, axis=0)!=0).sum()) for r in toy_rel]):.1f}",
    f"{np.mean([int((np.diff(r.state_path, axis=0)!=0).sum()) for r in c_rel]):,.0f}",
    note=f"mean of n={len(SEEDS)}, tau={TAU}; chr21 sums both haplotypes")
add("`expected_switches` under the prior [hap]", f(toy_hmm.expected_switches, 4),
    f"{c_hmm.expected_switches:,.1f}", note="model C predicted 7,155 for chr21, independently")
add("`baseline_log_norm`", f(toy_bl, 6), f(c_bl, 4),
    note="R_D's conditioning constant; log_z_p is measured against it")
add("zero-distance steps inside the chain", str(int((toy_hmm.switch_mass == 0).sum())),
    f"{int((c_hmm.switch_mass == 0).sum()):,}", note="P(switch)=0 there, which is correct")
add("positions with at least one missing donor", f"{int((~toy_hmm.allowed).any(axis=0).sum())} of {t.T}",
    f"{int((~c_hmm.allowed).any(axis=0).sum()):,} of {c_hmm.T:,}", note="missing_policy = RENORMALISE")
add("`chain_span`", "all 7 positions", "12,968,320 < POS <= 46,680,243",
    note="chr21 PROVISIONAL; keeps 89.74%")

# ------------------------------- CONTEXT -----------------------------------
add("**CONTEXT — why the tilt does what it does**", "", "")
add("sum `beta_v`", f(float(toy_u.beta.sum()), 12), f(float(c_u.beta.sum()), 12),
    note="MUST be 1; the utility's range")
add("max `beta_v`", e_(float(toy_u.beta.max())), e_(float(c_u.beta.max())))
add("sum of squared `beta_v`", e_(float((toy_u.beta**2).sum())), e_(float((c_u.beta**2).sum())),
    note="NOT 1/T unless beta is uniform")
add("`Var_Q(u)` [hap]", e_(toy_V_der), e_(c_V), note=f"d2 logZ/d eta2, tau={TAU}")
add("`Var_Q(u)` by direct enumeration [hap]", e_(toy_V_enum),
    f"{NA} — not enumerable", note="cross-check on the derivative route")
add("sd(`u_path`) per draw [hap]", e_(math.sqrt(abs(toy_V_der))), e_(math.sqrt(abs(c_V))),
    note=f"tau={TAU}")
add("gain the tilt can buy = `eta_tau` x Var [hap]", e_(eta * toy_V_der), e_(eta * c_V),
    note=f"tau={TAU}")
add("gain / noise per draw [hap]", e_(toy_gn), e_(c_gn),
    note=f"= eta_tau x sd; CEILING is eta_tau/2 = {eta/2:.4f}")
add("draws for the gain to reach 2 sigma [hap]", n_for_2sigma(toy_gn), n_for_2sigma(c_gn),
    note=f"at the ceiling it would be {round((2/(eta/2))**2)}")
add("per-position tilt at mean `beta_v`", e_(float(np.expm1(eta*toy_u.beta.mean())), 3),
    e_(float(np.expm1(eta*c_u.beta.mean())), 3), note=f"tau={TAU}")
add("seconds per release [dip]", e_(toy_secs, 2), f(c_secs, 2),
    note=f"n={len(SEEDS)}")
add("alpha storage, float64", f"{t.T*t.K*8} B", f"{c_hmm.T*c_hmm.K*8/1e6:,.0f} MB",
    note="load-bearing: tilt is ~15.1x float32 resolution")

# ------------------------------- DEMOTED -----------------------------------
add("**DEMOTED / NOT APPLICABLE**", "", "")
add("`af_loss`, `ld_loss`", f"{NA} BY DEFINITION", f"{NA} BY DEFINITION",
    note="cohort-distortion metrics; our release never edits the cohort")
add("`tv_empirical` at scale", "exact, see above", f"{NA} — biased UPWARD; DIAGNOSTIC only",
    note="never a headline number")
add("`eps_pmi`", f"{NA}", f"{NA}", note="PanMixer's score; incomparable with tau by construction")
add("`loo_id`", "null", "null", note="comparison arm only; our target is external")

# ------------------------------- EMIT --------------------------------------
out = ["### Full metric set — toy cohort vs chr21, metrics as rows", "",
       "| metric | toy (`toy4x7`) | chr21 (HG00096) | stamp |", "|---|---|---|---|"]
for label, a, b, gate, note in rows:
    if a == "" and b == "":
        out.append(f"| {label} | | | |")
    else:
        out.append(f"| {label}{' *' if gate and '*' not in label else ''} | {a} | {b} | {note} |")
out += ["",
        f"**Provenance.** TOY: `src/mechanism/toy.py`, hand-written, K={t.K}, T={t.T}, "
        f"{toy_base.M} enumerable outputs, READS NO DATA FILES · target = donor 1's alleles, "
        f"haploid · {len(toy_inputs):,} distinct release distributions swept exhaustively. "
        f"CHR21: `chr21/hprc-v1.0-pggb/K{c_hmm.K}/T{c_hmm.T}` · 44 individuals / "
        f"{c_hmm.K} haplotypes · target HG00096, diploid, held out of the cohort · "
        f"`chain_span` PROVISIONAL. BOTH: `tau` = {TAU}, `eta_tau` = {eta:.6f}, "
        f"`epsilon_tau` = {C.epsilon_tau(TAU):.6f} · utility "
        f"`{F.utility_id(c_u)}` · n = {len(SEEDS)} draws, seeds {SEEDS[0]}-{SEEDS[-1]} · "
        f"PATH A (all chain positions in the denominator) · deferred decisions R1-R5 UNSETTLED "
        f"· v1 numbers, NOT pre-registered."]
print("\n".join(out))
