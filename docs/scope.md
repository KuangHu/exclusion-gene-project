# Scope: ENTRY EXCLUSION ONLY

Set 2026-09-01. The database covers **entry exclusion** — the recipient-side
block on DNA translocation after mating-pair formation. Surface exclusion
(reduced stable mating-pair formation) is **out of scope**.

## Dropped by this cut

| gene | element(s) | why |
|---|---|---|
| `traT` | F, R100, pED208 | surface exclusion; target disputed; the TraT cryo-EM work argues against a specific partner interaction at all |
| `eexB` | R27 | surface exclusion, outer membrane (EexA is the entry factor, inner membrane) |
| `sfx` | IncC (pVCR94, R16a) | surface exclusion; targets TraN adhesin |
| `prgA` | pCF10 | surface exclusion, steric, no partner |
| `sea1` | pAD1 | surface exclusion |
| `ses` | pLS20 | surface exclusion |
| `mbeD` | ColE1 | mobilizable, not a conjugative T4SS; no VirB6 axis |

Elements retained only for a surface-exclusion gene are dropped entirely:
**R100, pED208, R16a, pCF10, pAD1, pLS20, ColE1** — and the negative controls
**pAM373, RSF1010** lose their comparator, so they move to `docs/negatives.md`
pending an entry-exclusion-appropriate negative set.

This also retires the whole `sfx` name-collision problem (§7) from the active
schema, since neither meaning of `sfx` is entry exclusion. The synonym record is
kept in `config/synonyms.tsv` so the collision stays documented.

## Retained: 10 elements, 10 entry-exclusion pairs

F · R64 · R621a · RP4 · pKM101 · R27 · IncC · SXT/R391 · ICEBs1 · pKPC_UVA01

The unit of record remains the PAIR (§0). Every retained row has an identified or
predicted **donor-side VirB6-family target**, which is the axis entry exclusion
actually operates on.
