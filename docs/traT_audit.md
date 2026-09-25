# Audit: F vs R100 TraT — derived from pinned accessions

**Run:** 2026-08-31 · **Method:** `scripts/45_traT_audit.py` (Biopython global
alignment, BLOSUM62, gap −11/−1) on CDS translations pulled from the pinned
GenBank records. Nothing here is transcribed from either literature claim.

## Sources

| | accession | protein_id | length |
|---|---|---|---|
| F | `AP001918.1` (99,159 bp, *E. coli* K-12) | `BAA97971.1` | 244 aa |
| R100 | `AP000342.1` (94,281 bp, *S. flexneri* 2b) | `BAA78882.1` | 243 aa |

## Result

**99.2% identical (242/244 aligned columns). Exactly two differences, both in the
first four residues of the signal peptide:**

```
F      1 MMKTKKLMMVALVSSTLALSGCGAMSTAIKKRNLEVKTQMSETIWLEPASERTVFLQIKN
          ||.|||||||||||||||||||||||||||||||||||||||||||||||||||||||||
R100   1 -MKMKKLMMVALVSSTLALSGCGAMSTAIKKRNLEVKTQMSETIWLEPASERTVFLQIKN
```

- alignment col 1 — F `M`, R100 gap
- alignment col 4 — F `T`, R100 `M`

**Residues 61–244 are 100% identical.** TraT is a lipoprotein; both differences
fall N-terminal of the lipobox, so the *mature* proteins are identical.

## What this does and does not refute

**Corrected 2026-08-31.** An earlier version of this note read the result as
bearing on the "residues 116-120" specificity window. That conflated two
independent claims and demoted the wrong one.

| Claim | Derived from | Status after this audit |
|---|---|---|
| A five-residue window at TraT **116-120** determines surface-exclusion specificity (Harrison et al. 1992) | **R6-5 (Sfx IV) vs ColB2-K98 (Sfx II)** | **Untouched.** This audit does not test it -- neither element is F or R100. |
| A single **G120A** substitution switches **F<->R100** specificity (Frost lab) | F vs R100 | **Dead.** The mature proteins are identical; there is no G120A difference to do the switching. |

F residues 116-120 are `RGYEG` and so are R100's (at its own 115-119). That is a
statement about F and R100 only. The 116-120 window was never an F/R100 claim,
and testing it requires the pair it was actually derived from.

## Extending to the pair that does test the window

**Half available, half blocked.**

* **R6-5 (Sfx IV)** -- `X52553.1`, TraT `CAA36788.1`. Available.
* **ColB2-K98 (Sfx II)** -- **not in NCBI.** `ColB2[All Fields] AND traT[All Fields]`
  returns nothing; the only ColB2 transfer-region submission, `U51860.1`, stops at
  TraL/TraE/TraK/TraB/TraP and never reaches *traT*. `ColB2-K98` as a term returns
  nothing at all.

So the window test cannot be completed from public sequence as things stand. The
ColB2-K98 *traT* sequence would have to come from the 1992 paper itself or from
the authors.

### A numbering hazard that must be settled before the test is run

`X52553.1` is annotated in **mature** numbering -- `sig_peptide` AA -20..-1,
`mat_peptide` "TraT lipoprotein (AA 1-223)". The F and R100 GenBank records are
annotated in **precursor** numbering (244 and 243 aa). The two conventions differ
by **20 residues**.

"Residues 116-120" therefore designates two different segments depending on which
convention the 1992 paper used. Precursor 116-120 in F = mature 96-100. Until that
is resolved from the paper, the window cannot be compared across these records at
all -- and this is on top of the separate one-residue F/R100 drift below.

## The numbering trap this exposes

Because of the single-residue N-terminal indel, **F residue *n* = R100 residue
*n*−1** for everything past position 4:

> F 116–120 (`RGYEG`) = **R100 115–119** (`RGYEG`)

A naive positional comparison at 116–120 in both proteins returns F `RGYEG` vs
R100 `GYEGA` and looks like a specificity difference. It is an off-by-one. This
is exactly the failure §5.3 legislates against, and it is why
`specificity_residues` must be stored as a range **against a named reference
accession**. The seed table stores the window twice, once per accession, rather
than once as a bare number.

## Consequence for the schema

`exclusion_group_label` for F and R100 cannot be separated by TraT sequence.
Any F/R100 group distinction must rest on TraS/TraG, which do differ
(F TraS 173 aa `BAA97970.1`; R100 TraS 159 aa `BAA78881.1`).
