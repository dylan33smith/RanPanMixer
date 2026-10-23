#!/usr/bin/env python
"""Every PanMixer attack that can run here, on our ARM 1 releases, beside their published values.

⚠⚠ WHAT "DIRECTLY COMPARABLE" CAN AND CANNOT MEAN HERE, because the paper decides that and not
us. Established 2026-10-14 by two independent extractions of the published PDF plus the preprint,
reconciled against both (`/data/ds85/RanPanMixer/runs/beagle_attack/published/PROVENANCE.md`):

  * **The published paper prints NO number for any attack except the linkage threshold.**
    `eps_private` = 0.002 and average utility loss 0.28 are printed; Fig 3B prints only its
    group sizes (All 55/45 Can Link / Cannot Link over 100 observations). Supplementary Figures
    S1-S7, which hold the MIA, Beagle, ablation and unique-variant results, are all CITED and
    all ABSENT from the file.
  * **The Beagle numbers exist only in the PREPRINT's Fig S2, and only as a curve**, which we
    digitized from its vector path (49 points, axes and caption independently verified).
    Comparing to it is comparing to OUR DIGITIZATION of their figure, never to a table.
  * **The MIA attack has no published number at all**, in either document, so our MIA row stands
    alone against their words ("the red line falls within the bulk") and nothing else.

So this table reports OUR numbers in full and places a published value beside them only where one
exists. A blank in the published column is the paper's silence, not a gap in this run.

⚠ THE PAPER'S BEAGLE PROTOCOL AND THEIR OWN CODE DISAGREE, so both are run and both reported:
  published Methods: "...with the 1000 Genomes dataset (excluding the target individuals) as the
    reference panel"  -> the `exclude-self` column, PAPER-FAITHFUL
  their `beagle_refinement.py`: passes no `excludesamples=` at all -> the `protocol` column
⚠ 39 of our 44 targets ARE in that panel, so without exclusion Beagle may copy the target's own
true row out of the reference. The 5 out-of-panel targets are a natural control on exactly that.

⚠ THE PRIMARY SUBSET IS `snp_only`, not `all`. The paper: "we first extracted SNPs overlapping
with those in the 1000 Genomes dataset. Using only the genotypes at these overlapping SNPs..."
That is their `snp_only` subset, 258,610 of 340,824 sites.

⚠ `xsol` IS ALL ONES, a declared substitution. Their metrics restrict to the cells their OPTIMIZER
selected; a path release has no selection step, so every released cell is in scope. Our
denominator is therefore the larger one.

    python scripts/arm1_attacks_results.py --json logs/B_arm1_attacks.json
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import statistics as st
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
ROOT = Path("/data/ds85/RanPanMixer/runs/beagle_attack/root")
PUB = Path("/data/ds85/RanPanMixer/runs/beagle_attack/published")
PM = REPO / "external" / "PanMixer"
CHROM = "21"
TAUS = ("0p00", "0p05", "0p10", "0p25", "0p50", "0p90", "0p99")
TAU_LABEL = {"0p00": "0", "0p05": "0.05", "0p10": "0.1", "0p25": "0.25",
             "0p50": "0.5", "0p90": "0.9", "0p99": "0.99"}
# exp_1 = beagle, exp_2 = gap_score, exp_3 = diploid_gap_score, exp_4 = MIA_privacy
EXP = {"beagle": "1", "gap_score": "2", "diploid_gap_score": "3", "MIA_privacy": "4"}


def rowdir(exp: str, row: int) -> Path:
    return ROOT / "experiments" / f"exp_{exp}" / "data" / f"chr{CHROM}" / str(row)


def load_changed(man: dict) -> dict:
    """Cells our release altered, per row. Their Fig S2 x-axis is a count of obfuscated variants.

    ⚠ A cell counts as ALTERED only where BOTH the truth and the release call an allele. An
    uncalled release row is not an obfuscation -- it is a position we never covered -- and
    counting it would inflate our x-coordinate against their curve.
    """
    orig_all = np.load(PM / "starting_data" / f"chr{CHROM}" / "pangenome.npy")
    subs = np.load(PM / "starting_data" / f"chr{CHROM}" / "pangenome_subjects.npy").tolist()
    # ⚠⚠ `pangenome_mask.npy` IS AN INDEX ARRAY, NOT A BOOLEAN MASK: 258,610 site
    # INDICES into the 340,824-site axis. This read it as `.astype(bool)`, which is
    # wrong twice -- it is True at every nonzero index (so 258,609 of 258,610 entries)
    # and it is the wrong LENGTH, which is what raised ValueError here. Their own
    # `accuracy_stats.py:43` builds the boolean the way it is built below, and
    # `gap_score.py:101` uses it as indices. Nothing published came from this
    # function: it crashed before returning, and `--subset snp_only` is its default.
    snp_idx = np.load(PM / "starting_data" / f"chr{CHROM}" / "pangenome_mask.npy")
    snp = np.zeros(orig_all.shape[1], dtype=bool)
    snp[snp_idx] = True
    out = {}
    for r in man["rows"]:
        i = r["row_id"]
        rel = np.load(rowdir(EXP["beagle"], i) / "new_haplotypes.npy")
        o = orig_all[subs.index(r["subject"])]
        both = (o >= 0) & (rel >= 0)
        diff = both & (o != rel)
        s2 = np.column_stack((snp, snp))
        out[i] = {"cells_called_both": int(both.sum()),
                  "cells_changed": int(diff.sum()),
                  "sites_changed": int(diff.any(axis=1).sum()),
                  "cells_changed_snp_only": int((diff & s2).sum()),
                  "sites_changed_snp_only": int((diff & s2).any(axis=1).sum())}
    return out


def collect(man: dict, changed: dict) -> dict:
    """Read every result file. ⚠ A missing file RAISES: a partial table is worse than none."""
    per = {}
    for r in man["rows"]:
        i, tau = r["row_id"], r["tau"]
        rec = {"subject": r["subject"], "tau": tau, **changed[i]}

        f = rowdir(EXP["beagle"], i) / "accuracy_results.json"
        if not f.exists():
            raise SystemExit(f"missing {f} -- the Beagle matrix is incomplete; refusing to "
                             f"report a partial comparison")
        rec["beagle"] = json.loads(f.read_text())
        f2 = rowdir(EXP["beagle"], i) / "accuracy_results_selfout.json"
        rec["beagle_selfout"] = json.loads(f2.read_text()) if f2.exists() else None

        for name in ("gap_score", "diploid_gap_score"):
            f = rowdir(EXP[name], i) / "gap_score_all.json"
            if not f.exists():
                raise SystemExit(f"missing {f}")
            rec[name] = json.loads(f.read_text())
        f = rowdir(EXP["MIA_privacy"], i) / "MIA_privacy.json"
        if not f.exists():
            raise SystemExit(f"missing {f}")
        rec["MIA_privacy"] = json.loads(f.read_text())
        per[i] = rec
    return per


def _med(v):
    v = [x for x in v if x is not None and not (isinstance(x, float) and math.isnan(x))]
    return st.median(v) if v else float("nan")


def published_curve():
    f = PUB / "beagle_S2_digitized.csv"
    if not f.exists():
        return []
    with f.open() as fh:
        return [{k: float(v) for k, v in row.items()} for row in csv.DictReader(fh)]


def interp_published(x: float, curve: list) -> float:
    if not curve:
        return float("nan")
    xs = [c["n_obfuscated_variants"] for c in curve]
    ys = [c["accuracy_difference_beagle_minus_majorallele"] for c in curve]
    return float(np.interp(x, xs, ys))


def main(argv=None) -> int:
    global ROOT, TAUS
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=Path, default=ROOT)
    ap.add_argument("--json", type=Path, default=None)
    ap.add_argument("--subset", default="snp_only",
                    help="PRIMARY is snp_only -- the paper's 'SNPs overlapping with 1000G'")
    a = ap.parse_args(argv)

    # ⚠ `rowdir` reads the module-level ROOT, so `--root` used to change ONLY the
    # manifest: a second mirror's labels would have been printed over the FIRST
    # mirror's result files, with no error. ROOT and the grid now both come from it.
    ROOT = a.root
    man = json.loads((a.root / "manifest.json").read_text())
    TAUS = tuple(man["taus"])
    TAU_LABEL.update({t: f"{float(t.replace('p', '.')):g}" for t in TAUS})
    changed = load_changed(man)
    per = collect(man, changed)
    curve = published_curve()
    S = a.subset

    by_tau = {t: [v for v in per.values() if v["tau"] == t] for t in TAUS}
    n_t = len(man["targets"])
    print(f"PanMixer's OWN attacks, run on OUR ARM 1 releases. {n_t} targets x "
          f"{len(TAUS)} `tau` = {len(per)} releases, chr21, K = 86 release on their K = 88 axis.")
    print(f"Metrics are ROWS, `tau` is the COLUMN. Median over {n_t} targets unless stated.")
    print(f"⚠ PRIMARY subset `{S}` -- the paper's 'SNPs overlapping with 1000 Genomes'.")
    print(f"⚠ Provenance `logs/B_arm1_attacks.md`; published values "
          f"`/data/ds85/RanPanMixer/runs/beagle_attack/published/PROVENANCE.md`.\n")

    hdr = ["metric"] + [TAU_LABEL[t] for t in TAUS] + ["PanMixer PUBLISHED"]
    print("| " + " | ".join(hdr) + " |")
    print("|" + "---|" * len(hdr))

    def row(label, fn, fmt="{:.6f}", pub=""):
        cells = []
        for t in TAUS:
            cells.append(fmt.format(fn(by_tau[t])) if by_tau[t] else "n/a")
        print(f"| {label} | " + " | ".join(cells) + f" | {pub} |")

    # ---------------- the Beagle reconstruction attack
    print(f"| **BEAGLE RECONSTRUCTION — paper-faithful (target EXCLUDED from panel)** | | | | | | | | |")
    has_self = any(v["beagle_selfout"] for v in per.values())
    if has_self:
        row("`beagle_accuracy` − `major_allele_accuracy` *  [THE PUBLISHED STATISTIC]",
            lambda g: _med([v["beagle_selfout"][S]["beagle_accuracy"]
                            - v["beagle_selfout"][S]["major_allele_accuracy"]
                            for v in g if v["beagle_selfout"]]),
            pub="**no number published**; preprint Fig S2 digitized, see below")
        row("`reconstruction_accuracy` *",
            lambda g: _med([v["beagle_selfout"][S]["reconstruction_accuracy"]
                            for v in g if v["beagle_selfout"]]),
            pub="no number published")
        row("`beagle_accuracy`",
            lambda g: _med([v["beagle_selfout"][S]["beagle_accuracy"]
                            for v in g if v["beagle_selfout"]]), pub="no number published")
        row("`major_allele_accuracy` [their baseline]",
            lambda g: _med([v["beagle_selfout"][S]["major_allele_accuracy"]
                            for v in g if v["beagle_selfout"]]), pub="no number published")
    else:
        print("| (exclude-self runs not yet present) | " + " | ".join(["n/a"] * len(TAUS))
              + " | |")

    print(f"| **BEAGLE RECONSTRUCTION — their CODE's protocol (no exclusion)** | | | | | | | | |")
    row("`beagle_accuracy` − `major_allele_accuracy`",
        lambda g: _med([v["beagle"][S]["beagle_accuracy"]
                        - v["beagle"][S]["major_allele_accuracy"] for v in g]),
        pub="their code passes no `excludesamples`")
    row("`reconstruction_accuracy`",
        lambda g: _med([v["beagle"][S]["reconstruction_accuracy"] for v in g]),
        pub="no number published")
    row("`obfuscated_accuracy` [the attacker's floor]",
        lambda g: _med([v["beagle"][S]["obfuscated_accuracy"] for v in g]),
        pub="no number published")

    # ---------------- the x-axis their figure is drawn against
    print(f"| **WHERE WE SIT ON THEIR FIG S2 x-AXIS** | | | | | | | | |")
    row("cells altered, `snp_only` [their '# of Obfuscated Variants']",
        lambda g: _med([v["cells_changed_snp_only"] for v in g]), fmt="{:,.0f}",
        pub="their curve spans 1,657 to 232,613")
    row("sites altered, `snp_only`",
        lambda g: _med([v["sites_changed_snp_only"] for v in g]), fmt="{:,.0f}", pub="")
    row("their DIGITIZED curve at our cell count [interpolated]",
        lambda g: interp_published(_med([v["cells_changed_snp_only"] for v in g]), curve),
        pub="**OUR digitization of their Fig S2**")

    # ---------------- linkage
    print(f"| **GAP SCORE — rarity-weighted GENOTYPE linkage** | | | | | | | | |")
    row("`genotype_g_to_gstar`",
        lambda g: _med([v["gap_score"]["genotype_g_to_gstar"] for v in g]), fmt="{:,.1f}",
        pub="no box statistic printed")
    row("targets LINKED (Gap Score >= 0) of 44 *",
        lambda g: sum(1 for v in g if v["gap_score"]["genotype_g_to_gstar"] >= 0),
        fmt="{:.0f}", pub="**55 of 100 Can Link** across their eps grid")
    print(f"| **DIPLOID GAP SCORE — the HAPLOTYPE version** | | | | | | | | |")
    row("`haplotype_g_to_gstar`",
        lambda g: _med([v["diploid_gap_score"]["haplotype_g_to_gstar"] for v in g]),
        fmt="{:,.1f}", pub="ablation, no numbers published")
    row("targets LINKED (>= 0) of 44 *",
        lambda g: sum(1 for v in g if v["diploid_gap_score"]["haplotype_g_to_gstar"] >= 0),
        fmt="{:.0f}", pub="")

    # ---------------- MIA
    print(f"| **MEMBERSHIP INFERENCE** | | | | | | | | |")
    row("`hap_original_rank` of 44 *",
        lambda g: _med([v["MIA_privacy"]["hap_original_rank"] for v in g]), fmt="{:.1f}",
        pub="**no number published**, either document")
    row("targets at rank 1 of 44 *",
        lambda g: sum(1 for v in g if v["MIA_privacy"]["hap_original_rank"] == 1),
        fmt="{:.0f}", pub="'the red line falls within the bulk'")
    row("`hap_bin_match_frac_vs_original`",
        lambda g: _med([v["MIA_privacy"]["hap_bin_match_frac_vs_original"] for v in g]),
        pub="no number published")

    # ---------------- the in-panel / out-of-panel control
    panel = set((PUB / "panel_samples.txt").read_text().split()) \
        if (PUB / "panel_samples.txt").exists() else set()
    if panel:
        ins = [t for t in man["targets"] if t in panel]
        outs = [t for t in man["targets"] if t not in panel]
        print(f"\n**THE PANEL-MEMBERSHIP CONTROL.** {len(ins)} of {len(man['targets'])} targets "
              f"are IN the 1000G reference panel and {len(outs)} are not "
              f"({', '.join(outs)}). Their code's protocol lets Beagle copy an in-panel "
              f"target's own true row out of the reference; the out-of-panel targets cannot be "
              f"copied at all, so the gap between the two groups is how much of the "
              f"'reconstruction' comes from the PANEL rather than from our release.\n")
        print("| group | " + " | ".join(TAU_LABEL[t] for t in TAUS) + " |")
        print("|" + "---|" * (len(TAUS) + 1))
        for nm, grp in (("IN panel (n=%d)" % len(ins), set(ins)),
                        ("OUT of panel (n=%d)" % len(outs), set(outs))):
            cells = []
            for t in TAUS:
                g = [v for v in by_tau[t] if v["subject"] in grp]
                cells.append(f"{_med([v['beagle'][S]['reconstruction_accuracy'] for v in g]):.6f}"
                             if g else "n/a")
            print(f"| `reconstruction_accuracy`, {nm} | " + " | ".join(cells) + " |")

    print(f"\n⚠ `sampled_gap_score` is NOT runnable here: it reads a PRE-stacker object from "
          f"`starting_data/chr{CHROM}/subjects/`, which is absent from the shipped tree. "
          f"THREE of their four privacy evaluators ran, never four.")
    print(f"⚠ `wged` and `pmi` are NOT Beagle statistics in their paper. 'WGED' appears once, as "
          f"*Weighted Graph Edit Distance*, a conceptual foil; PMI is the mechanism's own privacy "
          f"metric. Their `tools/beagle/beagle_stats.py`, which computes both, is DEAD CODE that "
          f"cannot run as shipped. The published statistic is an accuracy DIFFERENCE.")
    print(f"⚠ NOT a same-build comparison (decision F): their published numbers come from their "
          f"own vg-1.36 pipeline and ours from our vg-1.68 build. Dylan's instruction for this "
          f"run was to compare anyway and name the gap, which is what this line does.")
    print(f"⚠ Gap Score uses `>=` per Eq 6 and Fig 3A. The published Results text says '>0'. "
          f"⚠ The PDF's TEXT LAYER returns '>' for the Fig 3A glyph, which renders as '>=' -- "
          f"anyone grepping that file for the rule gets the wrong answer.")
    print(f"⚠ `eps_private` is 0.002 published and **0.001** in the preprint, with the definition "
          f"flipping from 'largest' to 'smallest'. No comparison here depends on it.")

    if a.json:
        a.json.write_text(json.dumps(
            {"manifest": {k: man[k] for k in ("n_sites", "n_blocks", "taus", "targets")},
             "subset_primary": S, "published_curve": curve, "per_row": per}, indent=1) + "\n")
        print(f"\nJSON -> {a.json}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
