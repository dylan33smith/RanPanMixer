"""The site axis: what row `i` means, and the digest that proves two artifacts agree.

Every preprocessing artifact is indexed by VCF record order. That is the join key
for the entire pipeline and nothing in the upstream tooling checks it, so a
re-sorted or re-called VCF renumbers everything silently. `sites.tsv` is the one
artifact that says what each row actually IS, and `axis_digest` is the fingerprint
every other artifact carries so a mis-join becomes an exception instead of a
plausible wrong number.

`sites.tsv` is also the only place the graph's nesting structure survives:
`VCFtoNP`-style converters keep position and genotype and drop the whole INFO
column, which is where LV and PS live.
"""

from __future__ import annotations

import csv
import hashlib
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np

COLUMNS = ["chrom", "pos", "ref", "alt", "n_alt", "lv", "ps"]

# The hypervariable chr21 record at 14,569,980 declares 90 alleles and its ALT
# string is ~28 MB -- one record holding a quarter of the chromosome's total ALT
# bytes. Python's csv module refuses fields over 128 KB by default, so reading
# sites.tsv fails on real data without this. Found by the dev slice, which was
# chosen to contain exactly that record.
csv.field_size_limit(min(sys.maxsize, 2**31 - 1))


@dataclass(frozen=True)
class Sites:
    """The site table, plus the two index spaces kept deliberately distinct.

    `n_sites` is the AXIS -- every array's first dimension, always the full
    record count. `T` is the CHAIN -- how many of those rows the sampler
    actually steps over, i.e. `in_chain.sum()`. Conflating them is how a
    provisional coordinate cut turns into a silently wrong denominator, so they
    are separate attributes and `T` is never stored, only derived.
    """

    chrom: np.ndarray       # (n_sites,) <U16
    pos: np.ndarray         # (n_sites,) int64
    ref: np.ndarray         # (n_sites,) object
    alt: np.ndarray         # (n_sites,) object  -- full comma-joined ALT string
    n_alt: np.ndarray       # (n_sites,) int32
    lv: np.ndarray          # (n_sites,) int16   -- snarl-tree level, -1 if absent
    ps: np.ndarray          # (n_sites,) object  -- parent snarl id, "" if absent
    in_chain: np.ndarray    # (n_sites,) bool    -- the chain_span MASK
    digest: str

    @property
    def n_sites(self) -> int:
        return len(self.pos)

    @property
    def T(self) -> int:
        return int(self.in_chain.sum())


def axis_digest(chrom, pos, ref, alt) -> str:
    """A canonical fingerprint of the site axis.

    Over (chrom, pos, ref, alt), NOT over file bytes: recompressing a VCF must
    not change it, while reordering two records must. That is the sensitivity we
    actually want, and it is why this is not just `sha256_file`.
    """
    h = hashlib.sha256()
    for c, p, r, a in zip(chrom, pos, ref, alt):
        h.update(f"{c}\t{int(p)}\t{r}\t{a}\n".encode())
    return h.hexdigest()


def write_sites_tsv(path: str | Path, rows) -> str:
    """Write sites.tsv and return its axis digest.

    `rows` is an iterable of (chrom, pos, ref, alt, n_alt, lv, ps).
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    h = hashlib.sha256()
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(COLUMNS)
        for chrom, pos, ref, alt, n_alt, lv, ps in rows:
            w.writerow([chrom, pos, ref, alt, n_alt, lv, ps])
            h.update(f"{chrom}\t{int(pos)}\t{ref}\t{alt}\n".encode())
    return h.hexdigest()


def load_sites(path: str | Path, in_chain_path: str | Path | None = None) -> Sites:
    """Read sites.tsv, and the chain mask beside it unless told otherwise.

    The mask is loaded BY DEFAULT. Forgetting to apply the coordinate cut should
    be the unusual path, not the easy one -- the cut is provisional and will
    move, and a caller who silently used the full axis would produce a number
    that looks fine and is against the wrong denominator.
    """
    path = Path(path)
    chrom, pos, ref, alt, n_alt, lv, ps = [], [], [], [], [], [], []
    with open(path, newline="") as fh:
        r = csv.DictReader(fh, delimiter="\t")
        if r.fieldnames != COLUMNS:
            raise ValueError(f"{path}: expected columns {COLUMNS}, got {r.fieldnames}")
        for row in r:
            chrom.append(row["chrom"]); pos.append(int(row["pos"]))
            ref.append(row["ref"]);     alt.append(row["alt"])
            n_alt.append(int(row["n_alt"]))
            lv.append(int(row["lv"]));  ps.append(row["ps"])

    chrom = np.asarray(chrom, dtype="<U16")
    pos = np.asarray(pos, dtype=np.int64)
    ref = np.asarray(ref, dtype=object)
    alt = np.asarray(alt, dtype=object)

    if in_chain_path is None:
        cand = path.parent / "in_chain.npy"
        in_chain = np.load(cand) if cand.exists() else np.ones(len(pos), dtype=bool)
    else:
        in_chain = np.load(in_chain_path)

    if len(in_chain) != len(pos):
        raise ValueError(
            f"in_chain has {len(in_chain)} entries but sites.tsv has {len(pos)}. "
            f"These were built against different VCFs."
        )

    return Sites(
        chrom=chrom, pos=pos, ref=ref, alt=alt,
        n_alt=np.asarray(n_alt, dtype=np.int32),
        lv=np.asarray(lv, dtype=np.int16),
        ps=np.asarray(ps, dtype=object),
        in_chain=in_chain,
        digest=axis_digest(chrom, pos, ref, alt),
    )
