# data.md — the map

Read before touching data, runs or paths.

**Rule.** Every number is traceable to a (`<ARTIFACT>`, `<GRAPH+COHORT>`,
`<TAU>`, `<UTILITY>`, n) tuple. **If it is not in this file, it is not quotable.**

**State values:** `OK` (present, current) · `PLANNED` (not built) · `MISSING`
(referenced but absent) · `STALE` (present but superseded) · `DEPRECATED`
(renamed, retained, never quote).

Verified against disk 2026-08-29. `tests/test_docs_contract.py` re-verifies.

> **STANDING NOTE (2026-08-29).** This project has **no data and no code yet**.
> Every row below except the documentation set is `PLANNED`. The layout is
> recorded now so that artifacts land in registered places from the first run
> instead of being retro-fitted. A `PLANNED` row is a commitment to a path, not
> a claim that it exists.

---

## 1. Storage layout

| root | contents | state |
|---|---|---|
| `docs` | the six-file documentation set | `OK` |
| `paper` | the proposal PDF and manuscript drafts | `OK` |
| `archive_docs` | read-only source material; never edited | `OK` |
| `tests` | code tests and the docs contract | `OK` |
| `src/ranpanmixer` | the mechanism: HMM, utility, sampler, attacks, eval | `PLANNED` |
| `data/graph` | pangenome graph + cohort haplotype paths | `PLANNED` |
| `data/target` | external target genomes and their mapped paths — **private input** | `PLANNED` |
| `artifacts/releases` | sampled sanitized paths, one directory per (graph, tau) | `PLANNED` |
| `artifacts/eval` | utility and attack measurements | `PLANNED` |
| `logs` | tmux run logs and status sentinels | `PLANNED` |
| `external/PanMixer` | upstream checkout, pinned by commit | `PLANNED` |

## 2. Upstream dependencies

| source | what | state |
|---|---|---|
| `https://github.com/G2Lab/PanMixer` | Reference implementation of the cohort HMM, LD-block segmentation and stay–switch transitions. We reuse its cohort model as the target-independent prior and **replace its final release rule**. Not yet cloned; no commit pinned. | `PLANNED` |
| Blindenbach, Soni & Gürsoy, *Privacy-preserving pangenome graphs*, bioRxiv 2026 | The PanMixer paper. Protects cohort contributors by editing the graph — a different threat model from ours. | reference |
| Eizenga et al., *Pangenome graphs*, Annu. Rev. Genomics Hum. Genet. 21:139–162, 2020 | Background on graph construction and path representation. | reference |

⚠ **No pangenome graph, cohort, or target genome has been selected.** The choice
of `G` and `D` fixes `T_blocks`, `K_states` and the entire runtime story, and it
is the first open decision in `docs/plan.md`.

## 3. Record schema — released path

**PLANNED.** Fields fixed now so that releases are self-describing from the first
run. A release file that cannot answer "at what tau, under what utility" is not
quotable, and unstamped output is the failure this schema exists to prevent.

| field | type | notes |
|---|---|---|
| `release_id` | str | ⚠ **not a unique key on its own.** One target may be released under several tau. The unique key is (`target_id`, `tau`, `utility_id`, `seed`). Collisions on `release_id` alone must raise, never overwrite. |
| `target_id` | str | pseudonymous target identifier |
| `tau` | float | the privacy parameter. Absent = the record is unusable. |
| `eta_tau` | float | stored, not recomputed at read time, so a calibration change is detectable |
| `utility_id` | str | identifies phi_t, its k, w, and the beta_t scheme |
| `graph_id` | str | graph + cohort + block segmentation |
| `seed` | int | the private random draw |
| `state_path` | int[T] | sampled Z_1:T |
| `log_z_p` | float | must satisfy 0 <= log_z_p <= eta_tau |
| `released_at` | date | one release per genome; re-releases compose |

## 4. Datasets — LIVE

| name | n | construction | leakage control |
|---|---|---|---|
| — | — | none yet | — |

## 5. Datasets — DEPRECATED (DO NOT USE)

| path | why | state |
|---|---|---|
| — | — | — |

## 6. Artifact registry

**Status key:** LIVE (quotable) · SUPERSEDED · INVALID (never quote)

| artifact | date | contents | status |
|---|---|---|---|
| `paper/Private_Genome_Path_Release.pdf` | 2026-08-29 | The proposal: mechanism, Theorem 1, Corollary 1, Algorithm 1, runtime analysis. The specification this project implements. | LIVE |
| `archive_docs/The_Six_File_Lab_Record.pdf` | 2026-08-29 | The documentation framework this repository follows. | LIVE |

## 7. Reproducibility contract

A result is reproducible here only if all five are recorded: the graph+cohort id,
`tau`, the utility id (including k and w), the number of draws, and the seed.
Two results are comparable only when the first three match — see the reporting
contract in `docs/plan.md`.

**The private input is the target path.** The `PLANNED` target directory must
never be committed, copied into an artifact directory, or quoted in any document.
Nothing in it is covered by the release guarantee; only the sampled output path
is. Written as a rule now, before the directory exists, because the first run is
where it would be broken.
