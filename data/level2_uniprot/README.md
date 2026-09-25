# Level 2 — UniProt / Swiss-Prot layer

`scripts/70_fetch_uniprot.py` → `uniprot_entries.tsv` + `uniprot_raw.json`

## Lookup route: identity by sequence, not by name

Searching UniProt by my own EMBL protein accession found **4 of 18** proteins,
because UniProt cross-references whichever EMBL record it curated from, not the
one I pinned. Resolving the frozen protein's **CRC64 checksum** through UniParc,
then querying UniProtKB for active entries on that UniParc id, found **15 of 18**.
A checksum hit is guaranteed to be the same protein.

**All 15 matched the level-1 frozen SHA. No cross-database sequence disagreement.**

| | count |
|---|---|
| reviewed (Swiss-Prot) | 4 — F TraS `P09129`, F TraG `P33790`, ICEBs1 YddJ `P96647`, ICEBs1 ConG `P96644` |
| unreviewed (TrEMBL) | 11 |
| no active UniProtKB entry | 3 — R621a ExcA, R621a TraY, pKPC TrbL-like (UniParc ids exist) |

Notably **RP4 TrbK is TrEMBL only** (`Q79AR9`) despite being the best-characterised
system in the set, while reviewed `trbK` entries exist for *Agrobacterium*
(`P54912`) and *Sinorhizobium* (`P55401`).

## FINDING 1 — the sequence conflict, verified end to end

`P09129` (F TraS) carries the only conflict in the set:

> `Sequence conflict  144..173  in Ref. 1 and 2`

Ref. 1 = J Mol Biol 198:1 (1987), PMID 3323526 · Ref. 2 = Frost, Microbiol Rev
58:162 (1994), PMID 7915817.

Verified directly against the frozen records, since one of the conflicting
submissions is already in the repo:

| record | protein | length |
|---|---|---|
| `AP001918.1` (complete F, 2000) | `BAA97970.1` | **173 aa** |
| `X06915.1` (1987 traS/traT submission) | `CAA30011.1` | **149 aa** |

```
residues 1-143     IDENTICAL
first divergence   residue 144   <- exactly UniProt's conflict start
1987 (149 aa)  ...FGAILVTVIMAFL DSI          RYR
2000 (173 aa)  ...FGAILVTVIMAFL IQLDIGRYQFVGVIDAINSYVKNKKLSRVK
```

A frameshift-class discrepancy truncating the protein 24 residues early — in the
**C-terminus of the entry-exclusion factor**. RP4 TrbK loses activity to an 8-aa
C-terminal truncation, so any functional work built on the 1987 construct was
missing the true C-terminal 30 aa. This is the single most consequential
provenance item found so far.

**My extractor had silently missed `X06915.1`'s TraS**: the CDS has no `/gene`
qualifier and its note reads "TraSp protein", which `TraS\b` cannot match. The
conflicting record was in the repo and invisible to the pipeline.

## FINDING 2 — UniProt's RP4 TrbK signal boundary is off by one

```
Q79AR9      Signal 1..21, Chain 22..69  -> mature 48 aa, +1 = G
haase1996   signal 1..22, mature 23..69 -> mature 47 aa, +1 = C
position     20 21 22 23 24  =  L  A  G  C  D
```

A lipoprotein's mature +1 must be the lipidated cysteine. Residue 22 is Gly and
cannot be lipidated; residue 23 is the Cys of the C23G mutant. **haase1996 is
correct and the UniProt feature is wrong** — it looks like an automatic SignalP
call that ignores the lipobox. Do not take level 2 over level 1 here.

## FINDING 3 — ConG has 6 annotated TM helices, not 7

`P96644` annotates six: 127-147, 157-177, 311-331, 333-353, 372-392, 401-421.
The curated context says "7 TM". The specificity loop **276-295** falls in the
long loop between TM2 and TM3, consistent with it being surface-exposed, but the
count needs reconciling against davis2021 before the topology claim is used.

## FINDING 4 — independent support for the TraG variable window, from a different data type

`P33790` (F TraG) annotates `Region 654..677 Disordered`, overlapping the
610-673 divergence window, inside `Topological domain 434..938 Periplasmic`.
R100 TraG (`Q6SIZ9`) carries the same feature at `656..675` — the 2-residue
offset expected from the 940 vs 938 length difference.

Disorder prediction is independent of sequence divergence, so this is corroborating
evidence of a genuine variable segment. It does **not** supply the 610-673
boundaries — UniProt's region is 654-677 — so the provenance of the exact range
is still open and still needs Audette 2007's body text.

`P33790` also carries `Chain None..938  Protein TraG*` — the periplasmic fragment
bragagnolo2022 characterises — but with an **undefined start**, so it cannot be
used to bound the domain either.

---

# Second pass — mined from the same JSON, no new requests

`scripts/71_mine_uniprot.py` → `uniprot_mined.tsv` · `scripts/72_signal_audit.py` → `signal_audit.tsv`

## FINDING 5 — 610-673 is NOT a structural slice

Every feature boundary in `P33790`: `1, 53, 54, 73, 74, 329, 330, 350, 351, 412,
413, 433, 434, 654, 677, 915, 923, 938`. **Neither 610 nor 673 is among them.**
The whole 434-938 stretch is a *single* periplasmic topological domain, so the
window is not a topology cut. Nearest boundary to 673 is 677 (4 off); nearest to
610 is 654 (44 off). `Chain Protein TraG*` has `modifier: UNKNOWN` for its start,
so it cannot bound anything either.

**The "domain slice" explanation is eliminated.** The window must come from
experimental mapping. Provenance still open; still needs Audette 2007 body text.

## FINDING 6 — exclusion proteins DO have dedicated domain families

§4 asserts these genes "most have no Pfam family" and that homology search
"saturates fast". **All 7 exclusion proteins with a UniProt entry carry domain
hits**, several of them exclusion-specific:

| element | gene | domain anchors |
|---|---|---|
| F | traS | `InterPro:IPR060798` (TraS) · `NCBIfam:NF010304` |
| R100 | traS | `InterPro:IPR018898` (**Eex_TraS**) · `Pfam:PF10624` (TraS) |
| RP4 | trbK | `InterPro:IPR027584` (TrbK_RP4) · `NCBIfam:TIGR04359` · `PROSITE:PS51257` |
| pKM101 | eex | `InterPro:IPR047937` (**Eex_IncN-like**) · `NCBIfam:NF033894` · `PROSITE:PS51257` |
| R64 | excA | `InterPro:IPR016389` (Exclusion-determining_protein) · `NCBIfam:NF033891` |
| SXT | eex | `NCBIfam:NF041429` (**EexR**) |
| ICEBs1 | yddJ | `InterPro:IPR028075` (DUF4467) · `Pfam:PF14729` |

These are HMM-searchable. The pipeline can anchor on them directly instead of
relying on position alone. Note F TraS and R100 TraS carry **different** InterPro
ids for the same gene (`IPR060798` vs `IPR018898`) — Swiss-Prot vs TrEMBL
annotation depth, not biology.

## FINDING 7 — "VirB6" is three different sequence families

`partner_family: VirB6` is a structural analogy, not a family. The actual anchors
split by MPF class:

| MPF class | elements | domain family |
|---|---|---|
| MPF_T | RP4 trbL, pKM101 traD | `Pfam:PF04610` (TrbL) · `IPR007688` |
| MPF_F | F/R100/IncC/SXT traG | `Pfam:PF07916` (TraG_N) · `IPR012931` |
| MPF_I | R64 traY | `IPR027628` (**DotA_TraY**) · `TIGR04346` |
| MPF_FA | ICEBs1 conG | `CDD:cd06261` (TM_PBP2) · `IPR000515` (**MetI-like**) |

ICEBs1 ConG is the starkest: **MetI-like is an ABC-transporter permease domain**,
not a conjugation family at all. A single VirB6 HMM will not find these; the scan
needs four anchors, chosen by MPF class.

## FINDING 8 — the TrbK signal error generalises to two more lipoproteins

`scripts/72_signal_audit.py` checks every entry mechanically.

| protein | UniProt Signal | lipobox Cys | verdict |
|---|---|---|---|
| RP4 TrbK `Q79AR9` | 1..21 | **23** | `SUSPECT_OFF_BY_ONE` |
| ICEBs1 YddJ `P96647` | 1..27 | **19** | `SUSPECT_MISMATCH` — 8 residues past the Cys |
| pKM101 Eex `Q79SE7` | **none annotated** | **15** | `LIPOBOX_UNANNOTATED` |
| ICEBs1 ConG `P96644` | 1..25 | none | correct — SPase I substrate |

All three lipoboxes are textbook:

```
pKM101 Eex   MKKLLLLIPFFLVA|C ...   mature 61 aa
ICEBs1 YddJ  MKNLFIFLSLMMMFVLTA|C ...   mature 108 aa
RP4 TrbK     MKKSNFIAVAALAAVMAASLAG|C ...   mature 47 aa
```

pKM101 Eex is the worst case: **UniProt annotates no signal peptide at all**, for
a protein whose defining paper (`pohlman1994`) names a lipid attachment motif in
its title. ICEBs1 YddJ is annotated as an SPase I substrate when it is a
lipoprotein.

Shared architecture — MKK/MKN start, hydrophobic stretch, lipobox cysteine at
15/19/23 — is a candidate anchor in its own right, and one that sequence identity
would never reveal (RP4 TrbK vs pKM101 Eex are at background identity).

## FINDING 9 — two seed sequences are flagged preliminary by UniProt

`CAUTION` comments: IncC TraG_C (`A0ABN7GQU1`) and SXT Eex (`A0ABT6U7Q2`) are both
*"derived from a whole genome shotgun (WGS) entry which is preliminary data"*.
Both are also the only two proteins in the set with **no AlphaFold model**.

## FINDING 10 — no experimental structures at all

13 of 15 have `AlphaFoldDB` models; 4 have `SMR`. **PDB cross-references: zero.**
Confirms §8's "no high-resolution structure for TraS, any Eex, TrbK, ExcA, EexC or
YddJ" — and means the roadmap's AF2 step is already largely done, not pending.

## FINDING 11 — `eexA`/`eexB` is a further name collision

`gene:eexB` in UniProt returns **`Q79AR9` — RP4 TrbK**, named "Entry exclusion
protein B"; `gene:eexA` returns `Q79AS0`, which is **RP4 TrbJ** (confirmed by
CRC64 against `BN000925.1`). Unrelated to R27's eexA/eexB (Gunton 2008). Two
independent systems use the same two gene names.

## R27 CLOSED by homolog anchoring — no paper needed

`AF250878.1` names neither `trhZ` nor `trhG`. Searching the record with UniProt
homologs from a related plasmid resolved both:

| gene | ORF | protein | length | coords | identity to homolog |
|---|---|---|---|---|---|
| eexA (trhZ) | `R0018` | `AAF69857.1` | 130 aa | 22505..22897 (+) | 78.1% to `A0A0H2XLF7` |
| trhG | `R0128` | `AAF69966.1` | 1329 aa | 113294..117283 (+) | 76.1% to `A0A0H2XKG9` |

Corroborated positionally, not just by score: `R0015/R0016/R0017/R0018` map exactly
onto the curated Z operon `trhO(015)-orf016(eexB)-orf017-trhZ(eexA)`, and `R0128`
sits 9 bp downstream of `trhH`, matching the Tra1 F operon `trhF-trhH-trhG`.
As a byproduct, **eexB = `R0016` / `AAF69855.1`, 279 aa** (out of scope, surface).

**R27's counterexample is now MEASURED: 90,396 bp apart, 15x outside the ±6 kb
window — and on the SAME strand.** That corrects the earlier `same_strand: false`:
the rule fails on distance and operon membership, not orientation.
