# MPF_F: entry criterion (in progress)

## The seed baseline caught a broken run

The first pass reported **0/5 seeds for all 13 families**. The scan was working --
each F-specific family admitted ~10,200-11,000 plasmids that MOB-suite
independently calls MPF_F, against its own count of 11,714. Only the seed check
was broken: PLSDB is RefSeq-based and the seed records are INSDC, so none of
AP001918.1, AP000342.1, CP033514.1, AF250878.1, KJ817376.1 is in PLSDB.

Four have RefSeq equivalents (NC_002483.1, NC_002134.1, NZ_CP033514.1,
NC_002305.1). SXT is an ICE, not a plasmid, so PLSDB legitimately lacks it and it
was injected into the cache exactly as RP4's TPA record was injected into the mash
sketch. Now sealed as **A12** in `scripts/lib/assertions.py`.

Without the seed column, "TraU admits 11,189 plasmids" would have been taken as a
clean result.

## Two groups, and the second is disqualified

| group | families | admits | MOB-suite MPF_T among admissions |
|---|---|---|---|
| F-specific | TraU, TraH, TraC_F_IV, TraG_N, TraE, TrbC_Ftype, F_T4SS_TraN, TraV, TraF | 15-17% of DB | 93-376 |
| cross-class | TrbI, P-loop_TraG, TrwB_AAD_bind, TraG-D_C | 25-35% | **6,019-9,034** |

The second group admits the MPF_T set wholesale and cannot serve as an entry
criterion.

**Partial convergence with the cross-class map -- 2 of 4, not 4 of 4.** An
earlier draft here claimed the disqualified families were "exactly" the mapped
ones. Checked: the map resolved PF03743, PF10412, PF11130; the disqualified set is
PF03743, PF10412, PF12696, PF19044. PF03743 and PF10412 are in both. PF11130 was
mapped but is the CHOSEN entry criterion, not disqualified. PF12696 and PF19044
were disqualified but had been dropped from the map by its one-to-one RBH collapse.

The mismatch is more informative than the convergence would have been:
**PF19044 (P-loop_TraG) and PF11130 (TraC_F_IV) are domains of the SAME protein**,
F TraC -- yet PF11130 admits at 93% MPF_F purity while PF19044 admits at 57%. TraC
carries a generic P-loop ATPase domain that also hits MPF_T VirB4 proteins, and an
F-specific domain that does not. The same holds for PF10412/PF12696 on TraD.

So the entry criterion is not "the right protein" but "the SPECIFIC DOMAIN of the
right protein". Choosing PF19044 instead of PF11130 would admit 18,223 plasmids
including 6,019 MPF_T -- the same protein, the wrong domain.

## Choosing among the F-specific families

**Ranking by "most conservative" (fewest admissions) is the wrong criterion**, and
was the first choice made here. MPF_T did not pick VirB4 because it admitted
fewest -- it picked VirB4 because it is the most CONSERVED component and stays
detectable in divergent systems. The explicit argument at the time was that using
VirB6 as entry would systematically lose divergent VirB6, which is where novel
families are.

By admission count the two smallest are F_T4SS_TraN (10,865) and TrbC_Ftype
(10,883). Both are wrong choices:

* **F_T4SS_TraN** -- PF06986 is a DOMAIN model, median 39.8% target coverage even
  on short hits. An entry criterion decides what can never be seen.
* **TrbC_Ftype** -- a pilin, VirB2-class. The measured two-tier structure puts
  pilus components in the high-failure tier (VirB2 GA failure 20.9%), so a pilin
  entry filter fails systematically on pilus-divergent plasmids.

**DECIDED: TraC_F_IV (PF11130), 11,120 admissions (15.4%), 5/5 seeds, 93% of
admissions independently MPF_F.** The VirB4 homolog, structural counterpart of the
MPF_T criterion, and specifically its F-specific domain rather than its generic
P-loop one.

Jaccard 0.950 against TrbC_Ftype gives an intersection of ~10,720: TraC_F_IV
admits ~400 that TrbC_Ftype misses and misses ~164 that it catches. The 400 are
plasmids whose pilin is undetectable, which is the failure mode the two-tier
result predicts for a pilin filter.

Recorded against the choice: **TraE recovers 94.3% of MOB-suite's MPF_F set
against TraC_F_IV's 88.7%.** If third-party recall were the criterion TraE wins.
It is not -- MOB-suite is coarser, gives no component coordinates, and MPF_T used
CONJscan -- but the gap is real and is not being quietly dropped.

Three classes now converge on the VirB4-family ATPase as entry criterion:
MPF_T VirB4 (PF03135), MPF_F TraC (PF11130), MPF_FA ConE (annotated
"VirB4-like ATPase ConE", 831 aa, in the newly extracted ICEBs1 seed).

## Open: the MPF_T overlap is not necessarily error

93-376 admissions per F-specific family are called MPF_T by MOB-suite. MOB-suite
is not ground truth -- MPF_T used CONJscan, and MOB-suite is coarser with no
component coordinates. These are either miscalls or genuine dual-system plasmids
carrying two T4SS. If the latter they belong in BOTH catalogues, and summing rows
across catalogues would double-count them. Needs a `dual_system` flag, not silent
overlap.

## Not yet done

* MPF_F architecture determination. MPF_T's was the sign of (VirB6 - VirB4),
  validated at 98.1% on plasmids with independent ground truth before being
  applied. MPF_F needs its own, validated the same way. **Do not presuppose the
  number of architectures** -- MPF_T ended with three, and R27 splits Tra1/Tra2
  into separate regions, so MPF_F may not be binary either.
* MPF_F slots, per class, not transferred from VirB5|VirB6.
