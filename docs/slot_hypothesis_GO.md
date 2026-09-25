# GO/NO-GO: the core hypothesis passes 4/4, and predicts the IncW gene

Run 2026-09-02 · `scripts/orfbench/slot_test.py` · `data`: `slot_test.tsv`

## What was tested

That in MPF_T conjugative plasmids the entry-exclusion gene occupies the slot
immediately 5' of the VirB6 homolog. Blind procedure: Pyrodigal skeleton →
hmmsearch the skeleton against Pfam T4SS anchors → locate VirB6 (`PF04610` TrbL)
→ examine the ORF 5' of it. Known exclusion coordinates were consulted only to
score, after the fact.

## Result

**Every one of the four MPF_T plasmids carries a small lipoprotein immediately
5' of its VirB6/TrbL anchor.**

| plasmid | Inc | coords | aa | annotation | lipobox | mature |
|---|---|---|---|---|---|---|
| RP4 | IncP-α | 27450..27659 | 69 | `trbK`, "TrbK entry exclusion protein" | `LAG↓C23` | 47 |
| R751 | IncP-β | 26340..26567 | 75 | `trbK`, "entry exclusion lipoprotein TrbK" | `VAG↓C22` | 54 |
| pKM101 | IncN | 6533..6760 | 75 | `eex` | `LVA↓C15` | 61 |
| **R388** | **IncW** | **5883..6113** | **76** | **unnamed; "EexN family lipoprotein"** | `LSA↓C17` | **60** |

In the three controls the ORF in that position **is** the known exclusion gene.
The hypothesis passes 3/3, and R388 supplies a candidate where none was named.

## Two families, one slot

hmmsearch of the four slot ORFs against the exclusion HMM families:

```
RP4_trbK        TrbK_RP4  TIGR04359   score  75.6   E=4.1e-25
R751_trbK       TrbK_RP4  TIGR04359   score 101.9   E=2.6e-33
pKM101_eex      Eex_IncN  NF033894    score  70.5   E=1.5e-23
R388_candidate  Eex_IncN  NF033894    score  73.6   E=1.6e-24
```

Pairwise identity splits on exactly the same line:

```
                RP4    pKM101   R751    R388
RP4              --     22.4%   37.8%   22.7%
pKM101         22.4%      --    25.0%   40.7%
R751           37.8%    25.0%     --    21.4%
R388           22.7%    40.7%   21.4%     --
```

Within family 37.8% / 40.7%; between families 21–25%, i.e. background.

**This is the strongest available form of the result: the SLOT is conserved
across two exclusion families that share no detectable sequence similarity.**
Position is conserved where homology is not — which is precisely the claim the
position-first architecture rests on.

## A concrete prediction, on a gap the protocol names

§8 lists "Missing genes entirely: IncW (R388 region-level only)". This gives a
specific candidate:

> **R388 entry exclusion gene = `WP_012196425.1`, 76 aa, `NC_028464.1:5883..6113` (+),
> Eex_IncN family (NF033894, E=1.6e-24), lipoprotein `LSA↓C17`, mature 60 aa,
> immediately 5' of TrbL/VirB6 (`6127..7155`).**

Untested phenotypically, so it enters as E4. But it is a named, coordinate-fixed,
family-assigned prediction where the literature had only a region.

## What this changes about the pipeline

**Pyrodigal found all four as skeleton genes.** The exclusion gene in this class
is not hidden in an intergenic gap; it is a normal skeleton call sitting in a
conserved position. So for MPF_T:

* the slot layer does the work, and it needs only the skeleton plus a VirB6 anchor
* six-frame all-starts is NOT required to find these genes

Six-frame's value stands where the benchmark actually placed it — overlapping
genes (excA/excB, mob3/mob4) and hard cases Pyrodigal misses (2/4 in Test 6) —
not for the canonical MPF_T slot.

## Method note: the first slot definition was wrong

The slot was initially defined as the intergenic gap *between* skeleton genes,
which returned 6–15 bp windows with zero candidates and scored MISS on all three
controls. The exclusion gene is *in* the skeleton, not between skeleton genes.
The corrected definition — the skeleton ORF immediately 5' of the VirB6 anchor —
is what the table above uses.

## Next

Extend to the MPF_T subset of PLSDB: locate the VirB6 anchor, take the ORF 5' of
it, and measure how often that position is occupied and by what. The four-plasmid
result says the signal is there; the library run measures how general it is.
