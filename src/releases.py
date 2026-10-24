"""Writing a release to disk, with its provenance, and refusing to write a wrong one.

⚠ WHY THIS MODULE EXISTS. `src/mechanism/ffbs.py` does not import `src/store.py` at
all, so until now a release lived only in memory: `Release.record()`'s docstring says
"the arrays go to disk through `store.save()`", which describes what a CALLER was
expected to do and not what any code did. The one script that persisted releases
(`scripts/run_mechanism_chr21.py`) used a bare `json.dumps`, so the first release on
record has no provenance sidecar, no `git_rev` and no axis digest. An artifact with no
provenance is not quotable by this project's own rule, so no genome-wide release was
quotable. That is what this closes.

⚠ IT IS A CHOKEPOINT, NOT A HELPER. `src/store.py` exists so that producing an artifact
and producing its record are the same function call; this does the same for the one
object the whole paper is about. The invariants below are the ones that are TRUE OF THE
WRITE and therefore invisible to `GenomeRelease`, which already checks everything
visible in memory -- field agreement, the utility family, duplicate chromosomes, the
share sum, and both bounds.

WHAT THIS ADDS THAT THE IN-MEMORY TYPES CANNOT SEE:

  * **A COLLISION MUST RAISE.** `Release`'s own docstring says the key is
    (`target_id`, `tau`, `utility_id`, `seed`) and "a collision must raise rather than
    overwrite". Nothing enforced it. Overwriting is the worst available failure here,
    because the second write looks exactly like the first and the guarantee is about
    the path that was RELEASED, not the one left on disk.
  * **THE ARM IS EXPLICIT, never inferred.** `CLAUDE.md`: never state a property of the
    target without naming its arm. The arm also decides the directory, so a mislabelled
    release lands in the wrong tree and is then pooled with incomparable numbers.
    `loo_id` must agree: arm 1 requires leave-one-out, arm 2 forbids it.
  * **THE AXIS DIGEST IS REQUIRED.** A release is a vector indexed by the site axis.
    Without the digest it can be silently joined to a different cohort build, which is
    the row-misalignment failure `store.check_axis` exists for.
  * **GENOME SCOPE IS ASSERTED AGAINST THE AUTOSOME SET.** ⚠ This is a real gap in the
    in-memory type: a `GenomeRelease` holding ONE chromosome with `share` 1.0 passes
    every check it makes -- the shares sum to 1, there are no duplicates, and the
    bounds hold -- and is still called a genome release. Only the write knows whether
    all 22 are there.
  * **CONSTRAINT 5 IS MADE VISIBLE ON DISK.** One release per genome. A `tau` grid and
    multiple seeds are an EVALUATION DEVICE and never a release procedure, so they are
    permitted but counted: every sidecar records how many sibling seeds already existed,
    and `assert_one_release_per_genome` is what anything quoting a release as DEPLOYED
    must call.
"""

from __future__ import annotations

import json
import re
from pathlib import Path

import store

ARMS = ("arm1", "arm2")
AUTOSOMES = tuple(f"chr{c}" for c in range(1, 23))
SCOPES = ("genome", "chromosome")


def release_root(root: str | Path | None = None) -> Path:
    if root is not None:
        return Path(root)
    here = Path(__file__).resolve().parent
    return here.parent / "data" / "releases"


def tau_dir(tau: float) -> str:
    """`0.1` -> `tau_0p10`, the form `CLAUDE.md`'s naming rule asks for.

    Two decimals wherever that is exact, so every pinned-grid directory keeps its
    name; otherwise the shortest exact decimal, so `0.999` -> `tau_0p999` and
    `0.99999` -> `tau_0p99999` (the exploratory high-`tau` points of 2026-10-22).

    ⚠ THE NAME MUST READ BACK AS EXACTLY `tau`. At two decimals alone, 0.995 and
    0.99 would both become `tau_0p99` and the second release would land in the
    first one's directory -- a silent collision between two different privacy
    levels, which is the one collision no reader could spot afterwards. This used
    to REFUSE any `tau` past two decimals; it now names it losslessly instead, and
    still raises on any name that does not round-trip.
    """
    two = f"{float(tau):.2f}"
    num = two if float(two) == float(tau) else repr(float(tau))
    name = "tau_" + num.replace(".", "p")
    back = float(name.removeprefix("tau_").replace("p", "."))
    if back != float(tau):
        raise ValueError(
            f"tau={tau!r} does not survive its directory name {name!r} (it reads "
            f"back as {back!r}). Two tau values would share a directory.")
    return name


def eta_dir(eta: float) -> str:
    """`500.0` -> `eta_500`, the directory for an ETA-DRIVEN release.

    ⚠⚠ WHY THIS CANNOT REUSE `tau_dir`. An eta-driven release has
    `tau = tanh(eta)`, which is exactly 1.0 in float64 past `eta` ~ 19 -- so EVERY
    high-eta release would land in one `tau_1p00` directory and overwrite the last.
    That is the same silent-collision failure `tau_dir`'s round-trip guard exists
    to stop, one level up, and it is why the eta grid gets its own namespace.

    The name is the shortest exact decimal of `eta` (so `eta_0`, `eta_500`,
    `eta_2821`, `eta_1p5`), and it must read back as exactly `eta` or this raises.
    `tau_*` and `eta_*` directories can never collide: the prefixes differ.
    """
    f = float(eta)
    if f < 0.0:
        raise ValueError(f"eta must be non-negative, got {eta!r}")
    num = repr(int(f)) if f == int(f) else repr(f)
    name = "eta_" + num.replace(".", "p")
    back = float(name.removeprefix("eta_").replace("p", "."))
    if back != f:
        raise ValueError(
            f"eta={eta!r} does not survive its directory name {name!r} (it reads "
            f"back as {back!r}). Two eta values would share a directory.")
    return name


def _name_for(rel) -> str:
    """The directory for a Release record, inferred when no name was given.

    ⚠ THIS IS THE FALLBACK. A caller that ran an `eta` grid should pass
    `point_name=eta_dir(eta)` explicitly, because the record does not store WHICH
    knob drove it and inference splits the namespace by accident: at `eta` = 1 the
    equivalent `tau` is representable, so inference would file it under `tau_*`
    while `eta` = 500 lands under `eta_*`. The inference below only stops the worst
    case -- every saturated release silently sharing `tau_1p00`.

    ⚠ `tau` == 1.0 is the tell: `ffbs.release(eta=...)` records `tanh(eta)`, which
    saturates to exactly 1.0 past `eta` ~ 19. `tau_dir(1.0)` is a legal name, so
    without this every eta-driven release would silently share `tau_1p00`.
    """
    eta = getattr(rel, "eta_tau", None)
    if eta is not None and float(rel.tau) >= 1.0:
        return eta_dir(float(eta))
    return tau_dir(float(rel.tau))


def point_dir(tau: float | None = None, *, eta: float | None = None) -> str:
    """The release directory for ONE operating point, named by whichever knob drove
    it. Exactly one of `tau` / `eta`, mirroring `ffbs.release`."""
    if (tau is None) == (eta is None):
        raise ValueError("pass exactly one of tau= or eta=")
    return tau_dir(tau) if eta is None else eta_dir(eta)


def _check_arm(arm: str, loo_id: str | None) -> None:
    if arm not in ARMS:
        raise ValueError(f"arm must be one of {ARMS}, not {arm!r}")
    # ⚠ Arm 1's guarantee holds ONLY under leave-one-out (standing constraint 1), and
    # arm 2's target is external by construction so there is nothing to hold out. A
    # release whose arm and `loo_id` disagree is one of those two claims being made
    # falsely, and which one depends on the direction.
    if arm == "arm1" and not loo_id:
        raise ValueError(
            "an arm 1 release has no loo_id. Arm 1's target IS a cohort member, so "
            "the guarantee is conditional on BOTH of its haplotypes leaving the donor "
            "panel (standing constraint 1). Without a loo_id this release was drawn "
            "against a panel containing its own target and the bound does not apply.")
    if arm == "arm2" and loo_id:
        raise ValueError(
            f"an arm 2 release carries loo_id={loo_id!r}. Arm 2's target is external "
            f"to G, D, every panel and the attack database, so there is nothing to "
            f"hold out and a loo_id means this is an arm 1 release mislabelled.")


def _dir_for(rel, arm: str, root: str | Path | None = None,
             point_name: str | None = None) -> Path:
    # ⚠ A record whose `tau` is 1.0 was driven by `eta` (tanh saturates), so it must
    # be named in the eta namespace or every high-eta release lands in one directory.
    return (release_root(root) / arm / rel.target_id / (point_name or _name_for(rel))
            / rel.release_id)


def existing_seeds(target_id: str, tau: float | None, arm: str,
                   root: str | Path | None = None, *,
                   eta: float | None = None) -> list[int]:
    """Seeds already written for this (target, `tau`, arm). The constraint-5 counter."""
    d = release_root(root) / arm / target_id / (
        eta_dir(eta) if eta is not None else tau_dir(tau))
    if not d.is_dir():
        return []
    seeds = set()
    for child in d.iterdir():
        m = re.search(r"_seed(\d+)", child.name)
        if m:
            seeds.add(int(m.group(1)))
    return sorted(seeds)


def save_release(rel, *, arm: str, axis_digest: str, inputs: dict[str, str],
                 root: str | Path | None = None, note: str = "",
                 point_name: str | None = None) -> Path:
    """Write one `Release` -- both arrays and its record -- each with provenance.

    Returns the directory. ⚠ RAISES if that directory already exists.
    """
    _check_arm(arm, rel.loo_id)
    if not axis_digest or not isinstance(axis_digest, str):
        raise ValueError(
            "axis_digest is required. A release is a vector indexed by the site axis, "
            "and without the digest nothing can refuse to join it to a different "
            "cohort build -- the row-misalignment failure store.check_axis exists for.")
    if not inputs:
        raise ValueError(
            "inputs is empty. An artifact with no recorded inputs cannot be shown to "
            "be stale, which is the whole purpose of the sidecar.")

    d = _dir_for(rel, arm, root, point_name)
    # ⚠ THE COLLISION CHECK, and it is the reason this function exists rather than a
    # bare np.save. The key is (target, tau, utility, seed, chrom); two releases with
    # that key are the SAME release, so a second write is either a redundant recompute
    # or -- worse -- a different draw silently replacing the one the guarantee is
    # about. Refuse and let the caller decide.
    if d.exists():
        raise FileExistsError(
            f"{d} already exists. The release key (target_id, tau, utility_id, seed"
            f"{', chrom' if rel.chrom else ''}) is already on disk, so this write "
            f"would replace a release rather than add one. If the intent is to REUSE "
            f"it, read it -- standing constraint 5 says one release per genome, cached "
            f"and reused. If the intent is a genuinely new draw, it needs a new seed, "
            f"and that makes this an EVALUATION sweep rather than a release.")

    siblings = existing_seeds(rel.target_id, rel.tau, arm, root)
    d.mkdir(parents=True)
    params = {
        "arm": arm, "scope": "chromosome" if rel.chrom else "whole",
        "tau": rel.tau, "eta_tau": rel.eta_tau, "seed": rel.seed,
        "utility_id": rel.utility_id, "graph_id": rel.graph_id,
        "chrom": rel.chrom, "share": rel.share, "ploidy": rel.ploidy,
        "loo_id": rel.loo_id,
        # ⚠ ON DISK, DELIBERATELY. Standing constraint 5 permits ONE release per
        # genome; more seeds at one tau is an evaluation device and must be legible
        # as such rather than inferred by counting directories later.
        "sibling_seeds_before_this_write": siblings,
        "is_the_only_seed_at_this_tau": not siblings,
    }
    for name, arr in (("state_path", rel.state_path),
                      ("released_path", rel.released_path)):
        store.save(d / f"{name}.npy", arr, inputs=inputs, params=params,
                   axis_digest=axis_digest,
                   note=note or f"{name} of release {rel.release_id}")
    store.save(d / "release.json", rel.record(), inputs=inputs, params=params,
               axis_digest=axis_digest,
               note="the quotable release record, docs/data.md section 4 schema")
    return d


def save_genome_release(gr, *, arm: str, axis_digest: str, inputs: dict[str, str],
                        root: str | Path | None = None,
                        scope: str = "genome",
                        point_name: str | None = None) -> Path:
    """Write every part of a `GenomeRelease`, then the genome-level record.

    ⚠ `scope` IS REQUIRED TO BE EXPLICIT because `GenomeRelease` cannot tell the
    difference. One chromosome with `share` 1.0 satisfies every check it makes --
    the shares sum to 1, there are no duplicates, both bounds hold -- and would be
    recorded as a release of a genome. Only the write knows the autosome set.
    """
    if scope not in SCOPES:
        raise ValueError(f"scope must be one of {SCOPES}, not {scope!r}")
    first = gr.parts[0]
    _check_arm(arm, first.loo_id)

    chroms = [p.chrom for p in gr.parts]
    if scope == "genome":
        missing = sorted(set(AUTOSOMES) - set(chroms),
                         key=lambda c: int(c[3:]))
        extra = sorted(set(chroms) - set(AUTOSOMES))
        if missing or extra:
            raise AssertionError(
                f"scope='genome' but the parts are not the 22 autosomes: "
                f"{len(gr.parts)} part(s), missing {missing or 'none'}, "
                f"unexpected {extra or 'none'}.\n"
                f"  ⚠ The shares summing to 1 does NOT catch this: a single "
                f"chromosome carrying the whole budget sums to 1 too, and `tau` "
                f"bounds the distinguishability of the WHOLE released path. A "
                f"partial set is a chromosome-scope release -- pass "
                f"scope='chromosome' and say so in every table.")

    # ⚠⚠ THE PARTS MUST INHERIT `point_name`. Without it each part falls back to
    # `_name_for`, and at `eta` = 0 the equivalent `tau` is 0.0 -- representable --
    # so the 22 parts would land in `tau_0p00/` (mixing with the pinned-grid
    # releases already there) while the genome record went to `eta_0/`. A release
    # split across two namespaces, with the genome record pointing at parts that
    # are not beside it. Measured 2026-10-24, before the 264-cell run.
    paths = [save_release(p, arm=arm, axis_digest=axis_digest, inputs=inputs,
                          root=root, point_name=point_name) for p in gr.parts]
    rec = dict(gr.record())
    rec["scope"] = scope
    rec["arm"] = arm
    out = (release_root(root) / arm / first.target_id
           / (point_name or _name_for(first))
           / f"{gr.release_id}.json")
    if out.exists():
        raise FileExistsError(
            f"{out} already exists; refusing to overwrite a genome release record.")
    store.save(out, rec, inputs=inputs,
               params={"arm": arm, "scope": scope, "tau": first.tau,
                       "seed": first.seed, "n_parts": len(gr.parts),
                       "parts": [str(p.name) for p in paths]},
               axis_digest=axis_digest,
               note=f"genome-level record over {len(gr.parts)} part(s); `tau` bounds "
                    f"the whole released path, so THIS is the quotable object and a "
                    f"per-chromosome log_z_p or u_achieved is not 'the release's'")
    return out


def assert_one_release_per_genome(target_id: str, tau: float, arm: str,
                                  root: str | Path | None = None) -> None:
    """⚠ Call this before quoting a release as DEPLOYED, never before an evaluation.

    Standing constraint 5 permits one release per genome because independent
    re-releases of the same path compose and degrade the guarantee. A `tau` grid and
    multiple seeds are a declared evaluation device over TEST targets; this is the
    check that stops an evaluation sweep being presented as a deployment.
    """
    seeds = existing_seeds(target_id, tau, arm, root)
    if len(seeds) > 1:
        raise AssertionError(
            f"{target_id} has {len(seeds)} releases at tau={tau} in {arm}: seeds "
            f"{seeds}. Independent re-releases of one path COMPOSE, so their joint "
            f"distinguishability is not bounded by tau (standing constraint 5). This "
            f"is a legitimate EVALUATION sweep over a test target, but it may not be "
            f"quoted as a release of a person -- report the joint_release_gap instead.")


def load_release(path: str | Path) -> dict:
    """Read a release record and VERIFY its sidecar before returning it."""
    p = Path(path)
    if p.is_dir():
        p = p / "release.json"
    store.verify(p)
    return json.loads(p.read_text())
