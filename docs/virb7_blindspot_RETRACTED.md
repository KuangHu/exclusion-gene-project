# RETRACTION: VirB7 is not a Pfam blind spot. Pyrodigal misses the gene.

## What was recorded

`data/anchors/anchor_set.tsv` and `docs/pipeline_v2.md` §1.1a Finding 3 stated:

> "VirB7 — no Pfam family hits VirB7 in any control (pKM101 `traN`, 48 aa — too
> small and divergent). VirB7 must not be used for slot definition."

That is wrong, and the reasoning built on it — that small T4SS components are
generally invisible to Pfam — was wrong with it.

## What is actually true

pKM101 TraN, 48 aa, scanned directly against Pfam-A at GA:

```
pKM101_traN   48 aa   HIT   P_T4SS_TraN (PF20898)   score 33.1
```

The family exists: **`PF20898` "P-type Type IV secretion system, TraN"**, and
there is a sister family **`PF06986` "F-type Type IV secretion system, TraN"**.

## Why the anchor scan missed it

`scripts/80_determine_anchors.py` labels **Pyrodigal calls** by exact coordinate
match against the deposited annotation, then reads off which Pfam family hits
each label. pKM101 `traN` is annotated at `U09868.1:7899..8045 (+)`.

**Pyrodigal does not call that ORF.** Its nearest call is `8035..8733`, a 232 aa
gene overlapping traN's 3' end. With no Pyrodigal call carrying the label
`traN`, no protein was scanned as VirB7, and the read-off reported "no family".

So the finding was an artefact of the **labelling step**, and underneath it is a
genuine **caller miss** — `docs/orf_caller_benchmark/` already documents that
Pyrodigal's boundary-case recall is ~50% and that it cannot express certain
short/overlapping genes. This is that limitation hiding a real T4SS component,
which is a sharper example than the benchmark's own cases.

## Consequences

**1. VirB7 has anchors.** `PF20898` (MPF_T) and `PF06986` (MPF_F). It may be used
for annotation. It still should not bound the exclusion slot, but for a different
reason than recorded: not "no family", rather that Pyrodigal frequently will not
supply a protein to scan.

**2. §7's pessimism was misplaced.** TraN — the anchor a `sfx`/traN-slot
extension would need — is well covered:

| protein | aa | Pfam at GA |
|---|---|---|
| F TraN | 602 | `PF06986` score **268.0** |
| IncC TraN | 936 | `PF06986` score **93.3** |
| pKM101 TraN | 48 | `PF20898` score **33.1** |
| R100 TraN | 300 | **no hit** — see below |

The traN slot is therefore **technically reachable**, contrary to what was
recorded. That does not by itself put `sfx` back in scope (it is surface
exclusion, and the entry-only cut still applies), but the stated technical
obstacle was not real.

**3. The two failure mechanisms must not be conflated.** An earlier note treated
VirB5's 21.4% failure as evidence that "small-protein anchors are the weak
point". They are different problems:

| | VirB5 / TrbJ | VirB7 / TraN |
|---|---|---|
| size | 226–258 aa — **not small** | 48 aa |
| failure | PF07996 built on IncN/IncW (87, 123); fails on IncP (17.7, 22.2) | Pyrodigal does not call the ORF |
| class | **family-coverage bias** | **caller recall** |
| fix | lower threshold + length conjunct (done) | six-frame within the slot window, or annotation-guided calling |

**4. R100 TraN is a loose end.** 300 aa with no Pfam hit, against F TraN at 602 aa
scoring 268. Two orthologues of the same gene should not differ two-fold in
length. Either the R100 annotation is truncated or my extraction took the wrong
feature. Flagged, not resolved.

## Bookkeeping

The `anchor_set.tsv` VirB7 row and `pipeline_v2.md` §1.1a Finding 3 both need
correcting. The blind-spot claim also propagated into
`docs/positional_assignment_NOGO.md` (VirB7 excluded from positional assignment)
and into §6.3's contamination argument — the §6.3 reasoning still holds, because
an unrecognised component in the slot is still possible, but the specific claim
"VirB7 has no family" must not be repeated.
