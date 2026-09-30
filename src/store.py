"""Writing artifacts, with provenance attached whether you want it or not.

This module exists for one reason. `CLAUDE.md` says every number must be
traceable to an artifact, and `docs/data.md` says a file not registered there is
not quotable. That is a rule someone has to remember. Routing every write
through `save()` turns it into something the code does on your behalf: you
cannot produce an artifact without also producing the record of where it came
from, because they are the same function call.

It is deliberately small. If this file grows past a couple of hundred lines it
has stopped paying for itself -- the point is a chokepoint, not a framework.
"""

from __future__ import annotations

import hashlib
import json
import subprocess
from pathlib import Path

import numpy as np

PROV_SUFFIX = ".prov.json"


def sha256_file(path: str | Path, chunk: int = 1 << 20) -> str:
    """Content digest of a file on disk."""
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        while block := fh.read(chunk):
            h.update(block)
    return h.hexdigest()


def _repo_root(start: Path) -> Path | None:
    """Walk up until a .git is found. Counting `parents[n]` is brittle: the first
    version used parents[2] from src/store.py, which lands OUTSIDE the repo, so
    every artifact recorded git_rev "unknown" while looking fine."""
    for d in [start, *start.parents]:
        if (d / ".git").exists():
            return d
    return None


def git_rev(repo: str | Path | None = None) -> str:
    """The state of the code that ran -- NOT just HEAD.

    ⚠ The first version returned `git rev-parse HEAD` with no dirty check, so an
    artifact built from a modified tree was attributed to the preceding commit.
    That is worse than recording nothing: 64 of 84 artifacts named a commit whose
    code provably could not have produced them (they carried parameters the
    committed code cannot emit). Provenance that confidently names the wrong cause
    is a trap.

    Returns `<sha>` for a clean tree and `<sha>-dirty+<12 hex of the diff>` for a
    modified one, so an artifact built mid-edit is self-evidently so and two such
    artifacts are distinguishable.
    """
    root = Path(repo) if repo else _repo_root(Path(__file__).resolve())
    if root is None:
        return "unknown"
    try:
        sha = subprocess.run(["git", "rev-parse", "HEAD"], cwd=str(root),
                             capture_output=True, text=True, check=True).stdout.strip()
        dirty = subprocess.run(["git", "status", "--porcelain"], cwd=str(root),
                               capture_output=True, text=True, check=True).stdout.strip()
        if not dirty:
            return sha
        diff = subprocess.run(["git", "diff", "HEAD"], cwd=str(root),
                              capture_output=True, text=True, check=True).stdout
        return f"{sha}-dirty+{hashlib.sha256(diff.encode()).hexdigest()[:12]}"
    except Exception:
        return "unknown"


def save(
    path: str | Path,
    obj,
    *,
    inputs: dict[str, str],
    params: dict,
    axis_digest: str | None = None,
    note: str = "",
) -> Path:
    """Write one artifact and its provenance sidecar.

    inputs      name -> sha256 (or a bare identifier for a non-file input).
                This is what makes a stale artifact detectable.
    params      every value that shaped the output. A parameter missing here is
                a parameter nobody can later prove was used -- and several of
                ours (`chain_span` above all) are explicitly provisional, so
                this is the field that matters most.
    axis_digest the site-axis digest this artifact is indexed by, when it has
                one. Row order is the join key for everything and nothing
                upstream checks it; carrying the digest is what lets a reader
                refuse to join two artifacts built against different VCFs.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    if isinstance(obj, np.ndarray):
        np.save(path, obj, allow_pickle=False)
        written = path.with_suffix(".npy") if path.suffix != ".npy" else path
        shape, dtype = list(obj.shape), str(obj.dtype)
    elif isinstance(obj, (dict, list)):
        written = path
        written.write_text(json.dumps(obj, indent=2, sort_keys=True) + "\n")
        shape, dtype = None, "json"
    elif isinstance(obj, str):
        written = path
        written.write_text(obj)
        shape, dtype = None, "text"
    else:
        raise TypeError(f"save() does not handle {type(obj)!r}; add it deliberately")

    prov = {
        "artifact": written.name,
        "sha256": sha256_file(written),
        "shape": shape,
        "dtype": dtype,
        "inputs": inputs,
        "params": params,
        "axis_digest": axis_digest,
        "git_rev": git_rev(),
        "note": note,
    }
    Path(str(written) + PROV_SUFFIX).write_text(
        json.dumps(prov, indent=2, sort_keys=True) + "\n"
    )
    return written


def attach_prov(path: str | Path, *, inputs: dict, params: dict,
                axis_digest: str | None = None, note: str = "") -> Path:
    """Attach a provenance sidecar to a file written elsewhere.

    `sites.tsv` is streamed by sites.write_sites_tsv so its digest can be
    computed in the same pass, which means it does not go through save(). It is
    also the most important artifact in the tree -- the one that defines what
    every row means -- so it must not be the one thing without a record. This is
    the narrow, deliberate exception, and it still produces the same sidecar.
    """
    path = Path(path)
    prov = {
        "artifact": path.name,
        "sha256": sha256_file(path),
        "shape": None,
        "dtype": "tsv",
        "inputs": inputs,
        "params": params,
        "axis_digest": axis_digest,
        "git_rev": git_rev(),
        "note": note,
    }
    Path(str(path) + PROV_SUFFIX).write_text(
        json.dumps(prov, indent=2, sort_keys=True) + "\n")
    return path


def load_prov(path: str | Path) -> dict:
    """The provenance record for an artifact, or a clear failure."""
    p = Path(str(path) + PROV_SUFFIX)
    if not p.exists():
        raise FileNotFoundError(
            f"{path} has no provenance sidecar. It was not written by store.save(), "
            f"so nothing records what produced it -- treat it as not quotable."
        )
    return json.loads(p.read_text())


def verify(path: str | Path) -> None:
    """Refuse an artifact whose bytes no longer match its recorded digest.

    ⚠ The sidecar's sha256 was write-only until 2026-09-29: nothing compared it
    against the file, so a silently mutated artifact passed every guard. A digest
    nobody checks is decoration.
    """
    prov = load_prov(path)
    actual = sha256_file(path)
    if actual != prov["sha256"]:
        raise ValueError(
            f"{path} does not match its provenance record:\n"
            f"  recorded sha256 {prov['sha256']}\n"
            f"  actual sha256   {actual}\n"
            f"The file changed after it was written. Rebuild it; do not override."
        )


def check_axis(path: str | Path, expected_digest: str) -> None:
    """Refuse to use an artifact built against a different site axis.

    The failure this prevents is the quiet one: a re-sorted or re-filtered VCF
    silently re-labels which variant is which, every downstream number is then
    about the wrong variants, and nothing errors. PanMixer has exactly this
    hazard and its own guard is vacuous -- both sides of its shape assertion
    equal `site_mask.sum()` by construction, so it can never fire.
    """
    got = load_prov(path).get("axis_digest")
    if got != expected_digest:
        raise ValueError(
            f"site-axis mismatch for {path}:\n"
            f"  artifact was built against {got}\n"
            f"  you are joining it against  {expected_digest}\n"
            f"These are different variant orderings. Rebuild, do not override."
        )
