# The slot tracks VirB5, not VirB6

## The result

The invariant is **eex immediately 3' of VirB5 (offset +1)**, holding across three
architectures. VirB6 is the element that relocates.

| group | eex at +1 from VirB5 | records | clusters | unique seqs |
|---|---|---|---|---|
| canonical | 83.9% | 380/453 | 177/200 (88.5%) | 87/105 (82.9%) |
| inverted IncI2 | 69.5% | 346/498 | 62/69 (89.9%) | 9/11 (81.8%) |
| inverted non-IncI2 | 79.6% | 43/54 | 21/29 (72.4%) | 8/11 (72.7%) |

Modal skeleton order, frame = VirB6 strand:

| group | modal order |
|---|---|
| canonical | VirB1 VirB2 VirB3 VirB4 **VirB5 eex VirB6** VirB8 VirB9 VirB10 VirB11 VirD4 |
| inverted IncI2 | **VirB6 VirB5 eex** VirB3 VirB4 VirB8 VirB9 VirB10 VirD4 VirB11 VirB1 |
| inverted non-IncI2 | **VirB6 VirB5 eex** |

`VirB5 -> eex` adjacency is preserved in all three; VirB6 moves from 3' of eex to
5' of VirB5. **A VirB6-based rule cannot be repaired by a constant offset** -- it
needs -1 in canonical and +2 in both inverted groups. The VirB5 rule is +1
everywhere.

### Two corrections to earlier versions of this file

1. **The 99.8% figure was inflated by nearest-hit selection.** It chose, among all
   Eex hits on a plasmid, the one NEAREST the anchor, then asked how often that
   one was adjacent -- close to circular. Using the best-SCORING hit gives 69.5%.
2. **The transcription frame was taken from VirB4's strand** while the
   canonical/inverted classification and the slot measurement both use the
   VirB5/VirB6 strand. Where VirB4 faces the other way the order string came out
   reversed, which made the non-IncI2 inverted group appear to put eex on the
   opposite side of VirB5. Corrected: both inverted groups agree.

### The n=1 concern, resolved

The VirB5-vs-VirB6 discrimination came entirely from IncI2, i.e. from one
architecture, so it was a single natural experiment however many records it
spanned. The 87 non-IncI2 inverted plasmids are an independent replicate and they
show the same rearranged order and the same +1 rule. Both are small at the level
that matters -- 11 unique Eex sequences each -- so the cluster-level agreement
(89.9% and 72.4%) is the figure to quote, never the record-level one.

## Why it was invisible until now

Canonical MPF_T order is `VirB5 -> [Eex] -> VirB6`. The exclusion gene is
simultaneously immediately 3' of VirB5 AND immediately 5' of VirB6, so the two
rules make identical predictions and no canonical plasmid can separate them. 85%
of the data is canonical.

IncI2 separates them. Its order is `[85aa] -> VirB6 -> VirB5 -> [Eex]`: VirB5 has
moved to the 3' side of VirB6, and **the Eex gene moved with VirB5**.

## How it was found, which was not by looking for it

The chain started from an artefact. Removing the VirB5-before-VirB6 assumption
admitted 554 previously-excluded plasmids; they were 97.8% "empty" between the
anchors, which dropped record occupancy 82.4% -> 73.8%. That drop was a proxy
artefact, and chasing it produced, in order:

1. the inverted class is real and strand-aware (10.8%, not ~50%)
2. it is 85% IncI2 and 97.6% MPF_T -- not mistyped
3. its 5'-of-VirB6 occupant fails every architectural test for an exclusion
   protein (no lipobox, hydrophobicity +0.83, no Pfam-A family, no ExcA homology)
4. yet 99.4% of those plasmids carry NF033894 (Eex_IncN) SOMEWHERE -- a higher
   rate than canonical plasmids at 80.7%
5. it sits at a fixed offset: +2 from VirB6, which is +1 from VirB5

Step 4 is what turned it. "The gene is missing here" and "the gene is elsewhere
here" look identical from inside the slot.

## The uncomfortable consequence

**The correct anchor is the worst-performing one.** From the completed ladder,
VirB5 fails on 27.1% of VirB4+ plasmids; VirB6 fails on 4.9%.

An earlier argument in `slot_definitions.md` treated the move to a VirB6-anchored
definition as a clear gain: "+1,521 plasmids, a 28% larger denominator, bought by
depending on a 4.9%-failure anchor instead of a 27.1%-failure one." That trade is
wrong as stated. A larger denominator measuring the wrong position is worse than a
smaller one measuring the right position. The VirB6 rule is a **proxy that is
correct only under canonical order**, and its error is not random: it fails
systematically on rearranged skeletons, which are the subpopulation most likely to
carry divergent families.

## What replaces it: determine architecture first, then apply the offset

The dilemma was false. It is not necessary to choose one anchor globally:

    1. use the CORE anchors (VirB4/8/9/10/11, VirD4 -- all ~4.1-5.0% failure)
       to establish whether the skeleton is canonical or rearranged
    2. canonical  -> eex at +1 from VirB5, equivalently -1 from VirB6 (equivalent
                     here, 83.9% vs 82.7%, so the weak anchor is not required)
    3. rearranged -> eex at +1 from VirB5, which is +2 from VirB6

The 1,286 `no_VirB5_only` plasmids need not be discarded, and the VirB6 proxy's
error rate becomes bounded rather than unknown: inverted plasmids are 589/5,669 =
**10.4%** of the pairable set, so ~10% is an UPPER bound on the proxy's error on
the VirB5-negative subset, and lower to the extent architecture is recoverable
from the core anchors.

**The original error was not "used VirB6". It was assuming canonical order without
checking.** The `VirB5 -> eex` block being intact is what makes the architecture-
first approach work: the skeleton rearranges as blocks, so the core anchors carry
information about where VirB5 should be even when PF07996 fails to detect it.

## Tension with the project's own adjacency rule

The position table concluded that exclusion genes sit **adjacent to their
partner**, not adjacent to VirB6 as such: IncC `sfx` abuts `traN` (its target),
pCF10 `prgA` abuts `prgB`, pAD1 `sea1` abuts `asa1`.

Entry-exclusion targets recorded so far are VirB6-class: SXT TraG aa 606-608,
ICEBs1 ConG loop 276-295, IncC TraG_C. So if the IncI2 exclusion gene tracks
VirB5 rather than VirB6, one of two things must give:

1. it violates the adjacent-to-partner rule, or
2. **its target is not VirB6.**

Option 2 is the interesting one. VirB5 is the pilus-tip adhesin. An exclusion
protein targeting VirB5 rather than the VirB6 channel component would sit closer
to surface exclusion than to entry exclusion mechanistically -- which connects
directly to the scope question in `pipeline_v2.md` 7, since the entry-only scope
was what removed the TraDIS-class targets in the first place.

**The published rule carries the same ambiguity as our pipeline did.** The
statement that IncP/IncW/IncN exclusion genes abut a VirB6 homolog is
underdetermined for the canonical case: RP4 `trbK` lies between `trbJ` (VirB5) and
`trbL` (VirB6) and is adjacent to BOTH. The literature is not wrong; it is
unidentifiable from canonical plasmids, for exactly the reason our own two rules
were unidentifiable.

**IncI2 is the natural experiment that separates them.** If it holds, this is a
substantive correction to a published rule rather than an internal pipeline fix.

Caveat before any of that is used: the wording attributed to the literature here
comes from project notes and must be checked against the source before it is
repeated as a quotation, and the mechanism claim requires the target to be
identified, not inferred from position.
