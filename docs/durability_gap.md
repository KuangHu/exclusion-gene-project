# Durability gap: ICEBs1's authoritative record is not in the repo

`ICEBs1` has `accession_authoritative = CP171645.1` — the *B. subtilis* 168
chromosome, **10.7 MB**, larger than the entire rest of the project. It lives at
`/global/scratch/users/kh36969/exclusion_gene/hosts/`, not in `data/seed/genbank/`.

`scripts/63_verify_targets.py` caught this as `RECORD_MISSING` on its first run.

## Why it matters

§7 requires a frozen accession set per release. Every other element's
authoritative record is in-repo and version-controlled. ICEBs1's is on a shared
Lustre filesystem with no quota, no backup and no version history — if it is
cleaned up, the ICEBs1 coordinates in the database become unverifiable.

## The fix, not yet done

Extract the ICEBs1 element span from the chromosome and freeze **that** as the
authoritative record, with the host chromosome demoted to `accession_related`.
This is what `accessions.status: CHROMOSOMAL` was always meant to imply.

Blocked on the element boundaries. ICEBs1 integrates at *trnS-leu2*; the span
should be taken from the published annotation rather than inferred from the
`conB..yddM` gene range, since the integrase and attachment sites sit outside
the operon and would be silently truncated.

Until then the verifier reads through to scratch, and the coordinates are
correct but **not durably anchored**.
