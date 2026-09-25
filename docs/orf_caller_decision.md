# DECISION: the ORF/CDS caller is Pyrodigal anon. One caller, one config.

Decided 2026-09-02. Implementation: `scripts/lib/orf_caller.py`.
Evidence: `docs/orf_caller_benchmark/` · `docs/slot_hypothesis_GO.md`.

```
meta / anon = True      min_gene = 90 nt      closed = True      table = 11
```

Changing any of those invalidates the benchmark and must be justified there.
Alternatives are not a parameter of this module; adding one is a decision to
re-open with new evidence.

## Why anon rather than trained

Trained mode has a **hard floor, not graceful degradation**:
`sequence must be at least 20000 characters`. It raises on `U09868.1`
(12,416 bp) and `J01566.1` (6,646 bp). Fragment records and small plasmids are
routine input, so trained mode would fail on part of any batch — silently if the
exception were caught. `assert_anon_only()` makes that unreachable. Trained mode
also had the worst recall in the benchmark (7/10 vs 9/10).

## Why not six-frame as the primary layer

The go/no-go test settled this. In all four MPF_T controls the exclusion gene is
an **ordinary skeleton call** — Pyrodigal anon found RP4 `trbK`, R751 `trbK`,
pKM101 `eex` and the R388 candidate at exact coordinates. The canonical slot does
not need enumeration, so the enumerative layer earns nothing on the main path.

## Limits accepted deliberately, not overlooked

| limit | case | consequence |
|---|---|---|
| one ORF per stop codon | R64 `excA`/`excB` share a stop at `AP005147.1:76983`, starts 219 nt apart | `excB` is not expressible at any threshold |
| nested different-frame gene | ColE1 `mob4` inside `mob3` | missed entirely |
| boundary-case recall ~50% | 2/4 on genes a production pipeline missed | vs 90% on genes one had already found |

Accepted because the primary use — skeleton definition and slot occupancy — does
not depend on them, and because determinism is worth more here than marginal
recall: on two deposits of one **identical** sequence, every caller tested was
byte-identical while the **GenBank annotations differed by 4 of 47 CDS**. The
annotation is the unreliable layer, not the caller.

## Circularity caveat

Pyrodigal has no notion of a circular replicon. With `closed=True` a gene
spanning the origin is not called; with `closed=False` it appears as two partials
at the ends. For a gene of interest near coordinate 1 or near the end of a
plasmid, re-run on a rotated sequence rather than trusting the default.

## Regression sentinel

`python3 scripts/lib/orf_caller.py` re-calls the four MPF_T slot genes and
requires exact coordinates and lengths on all four. Currently **PASS**.
This is the sentinel that must be green before any batch run.

## Cleanup performed

| | |
|---|---|
| kept as evidence | `docs/orf_caller_benchmark/` — RESULTS.md + 6 result TSVs |
| archived, recoverable | `archive/orfbench_2026-09-02/` — 8 harness modules (orfipy, StORF, six-frame adapters, test drivers) |
| removed | `scripts/orfbench/`, all `__pycache__` |
| kept as infrastructure | the HMM models in `/global/scratch/users/kh36969/exclusion_gene/hmm/` — slot anchors (`skeleton_anchors.hmm`) and exclusion families (`exclusion_anchors.hmm`) are pipeline components, not benchmark leftovers |

`orfipy` and `StORF-Reporter` remain installed in `claude-env` but the project no
longer imports them. Nothing in the main pipeline ever referenced the benchmark,
so its removal touched no other module.
