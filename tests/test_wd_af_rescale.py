"""GATES on `wd_af`, the number PanMixer's figures actually report.

`wd_af` is `af_loss` times two. That is a one-line transform, which is exactly why it
needs a gate: the failure mode is not a crash but a column that is silently off their
published scale by 2x, or a stratum mapped to the wrong panel. Either makes a
head-to-head table wrong while every number in it looks plausible.

The mapping is NOT ours to choose — it is read off
`external/PanMixer/tools/common/gather_results.py`, and these tests check OUR table
against THEIR source rather than against our own belief about it.
"""

from __future__ import annotations

import re
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
TABLER = REPO / "scripts" / "arm1_af_ld_results.py"
THEIRS = REPO / "external" / "PanMixer" / "tools" / "common" / "gather_results.py"
AF_LOSS = (REPO / "external" / "PanMixer" / "tools" / "downstream" / "utility_in"
           / "af_loss.py")

sys.path.insert(0, str(REPO / "scripts"))


def _their_pairings() -> dict[str, str]:
    """Parse `df["wd_af_x"] = 2 * np.array(df["af_loss_y"])` out of THEIR source.

    Returns {wd column -> the af_loss DATAFRAME column it doubles}.
    """
    src = THEIRS.read_text()
    pairs = dict(re.findall(
        r'df\["(wd_af[a-z0-9_]*)"\]\s*=\s*2\s*\*\s*np\.array\(df\["([a-z0-9_]+)"\]\)',
        src))
    assert pairs, "found no wd_af derivations in their gather_results.py"
    return pairs


def _their_json_keys() -> dict[str, str]:
    """{their af_loss dataframe column -> the af_loss.json key it was filled from}."""
    src = THEIRS.read_text()
    # af_loss_all.append(results["all"]) ... then df["af_loss"] = af_loss_all
    listvar = dict(re.findall(r'(\w+)\.append\(results\["([^"]+)"\]\)', src))
    colvar = dict(re.findall(r'df\["([a-z0-9_]+)"\]\s*=\s*(\w+)\s*$', src, re.M))
    return {col: listvar[var] for col, var in colvar.items() if var in listvar}


def test_the_parsers_actually_found_their_derivations():
    """⚠ A gate whose parser returns nothing passes vacuously. Pin the floor."""
    pairs, keys = _their_pairings(), _their_json_keys()
    assert len(pairs) >= 7, f"expected >=7 wd_af columns in their source, got {pairs}"
    assert len(keys) >= 7, f"expected >=7 af_loss columns in their source, got {keys}"
    assert "wd_af" in pairs and pairs["wd_af"] == "af_loss"


def test_our_map_matches_their_source_stratum_for_stratum():
    """THE GATE. Every row we label `wd_af_x` must double the stratum THEY double.

    A mismapping here would put, say, the complex-variant distortion under the
    SNPs-only panel and compare it to the wrong published number.
    """
    from arm1_af_ld_results import WD_MAP

    pairs, keys = _their_pairings(), _their_json_keys()
    wrong = []
    for their_col, our_stratum, _label in WD_MAP:
        if their_col not in pairs:
            wrong.append(f"{their_col}: we emit it, their source does not derive it")
            continue
        their_json_key = keys.get(pairs[their_col])
        if their_json_key != our_stratum:
            wrong.append(f"{their_col}: we double af_loss[{our_stratum!r}], they "
                         f"double af_loss[{their_json_key!r}]")
    assert not wrong, "WD_MAP disagrees with their gather_results.py:\n  " + \
        "\n  ".join(wrong)


def test_we_cover_every_wd_column_they_derive():
    """If they plot a panel we do not emit, our table is missing a column."""
    from arm1_af_ld_results import WD_MAP

    ours = {t for t, _, _ in WD_MAP}
    missing = set(_their_pairings()) - ours
    assert not missing, f"their source derives {sorted(missing)} and we do not emit it"


def test_the_published_reference_is_only_set_for_panels_they_plot():
    """Their Fig 4A plots three AF panels plus ld_loss. A published value attached to
    a panel they never plot would be invented."""
    from arm1_af_ld_results import THEIR_FIG4A

    plotted = {"wd_af", "wd_af_snp_only", "wd_af_snps_0_05"}
    assert set(THEIR_FIG4A) == plotted, (
        f"published references are set for {sorted(THEIR_FIG4A)}; their figure plots "
        f"{sorted(plotted)} (plot_util/utility_in.py subplot_columns)")


def test_the_maf_strata_are_disjoint_bins_not_thresholds():
    """⚠ Their own figure labels `wd_af_snps_0_05` as "MAF < 0.05". It is NOT.

    `af_loss.py` defines half-open BINS, so that panel is [0.01, 0.05) — the middle
    bin — not every SNP below 0.05, which would be the union of two bins. Our labels
    must say the bin, or we inherit their mislabel.
    """
    src = AF_LOSS.read_text()
    assert "snps_maf_0_01 = global_maf < 0.01" in src
    assert "snps_maf_0_05 = (global_maf >= 0.01) & (global_maf < 0.05)" in src, (
        "their MAF binning changed; re-read af_loss.py before trusting any stratum "
        "label in our table")

    from arm1_af_ld_results import WD_MAP

    lbl = {t: l for t, _, l in WD_MAP}
    assert "[0.01, 0.05)" in lbl["wd_af_snps_0_05"], (
        f"our label for the middle bin is {lbl['wd_af_snps_0_05']!r}; it must state "
        f"the BIN, not reproduce their 'MAF < 0.05'")


def test_the_emitted_table_really_doubles(tmp_path):
    """End to end: run the tabler and check each wd row is 2x its af_loss row.

    ⚠ This is the one test that would catch a transform that is right in the map and
    wrong in the arithmetic.
    """
    mirror = Path("/data/ds85/RanPanMixer/runs/beagle_attack_eta_uniform/root")
    if not (mirror / "manifest.json").exists():
        pytest.skip("eta mirror not present")
    r = subprocess.run([sys.executable, str(TABLER), "--root", str(mirror)],
                       capture_output=True, text=True, timeout=900)
    assert r.returncode == 0, f"tabler failed: {r.stderr[-1500:]}"

    af_rows, wd_rows = {}, {}
    for line in r.stdout.splitlines():
        m = re.match(r"\| `af_loss` \[([^\]]+)\] \| (.+) \|$", line)
        if m and m.group(1) not in af_rows:          # OURS block comes first
            af_rows[m.group(1)] = [float(x) for x in m.group(2).split(" | ")]
        m = re.match(r"\| `(wd_af[a-z0-9_]*)` = 2x \[([^\]]+)\][^|]*\| (.+) \|$", line)
        if m and m.group(1) not in wd_rows:
            wd_rows[m.group(1)] = (m.group(2),
                                   [float(x) for x in m.group(3).split(" | ")])
    assert len(wd_rows) >= 7, f"parsed only {len(wd_rows)} wd rows from the table"

    for wd, (stratum, vals) in wd_rows.items():
        base = af_rows.get(stratum)
        assert base, f"{wd} cites stratum {stratum!r} with no af_loss row"
        assert len(base) == len(vals), f"{wd}: {len(vals)} cells vs {len(base)}"
        for i, (b, v) in enumerate(zip(base, vals)):
            # ⚠ TOLERANCE, and why it is not a weakened gate. Both rows are printed at
            # 6 dp, so each carries up to 5e-7 of rounding. Doubling the ROUNDED
            # af_loss can therefore differ from the rounded DOUBLE by up to
            # 2*5e-7 + 5e-7 = 1.5e-6 with the arithmetic perfectly correct, which is
            # what a 1e-6 bound tripped on. The defects this test exists to catch are
            # orders of magnitude bigger: a missing factor of 2 is off by ~3.7e-3, and
            # a mis-mapped stratum by ~1e-3 or more. 2e-6 separates the two cleanly.
            assert abs(v - 2.0 * b) < 2e-6, (
                f"{wd} column {i}: printed {v}, but 2 x af_loss[{stratum}] = "
                f"{2.0 * b} -- off by {abs(v - 2.0 * b):.2e}, far beyond 6-dp rounding")
