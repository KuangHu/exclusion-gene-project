# GA thresholds fail systematically on these proteins

## The claim

For mobile-element proteins of this length and divergence, a pHMM library's
gathering cutoff systematically loses divergent family members, and pairwise
profile search (phmmer against confirmed members) is a required complement rather
than an optimisation.

This is stated as a methods conclusion because it has now been the true cause of
four separate findings that each looked structural at first.

## The four instances

| # | observation | what it looked like | what it was |
|---|---|---|---|
| 1 | `PF07996` scores IncP TrbJ at 17.7 (RP4) and 22.2 (R751) against GA 24.7 | "VirB5 is undetectable in IncP" | GA cutoff; both pass easily by E-value (7.9e-07, 3.4e-08) |
| 2 | MMseqs2 returns zero cross-hits among the four seeds at the proposed parameters | "the seed families are unrelated" | k-mer prefilter fails on 69-76 aa proteins; phmmer recovers both families |
| 3 | `PF07996` scores **0 of 1,221** slot candidates at GA | "these candidates are not VirB5" | phmmer against 228 confirmed VirB5 matched **162 of the 228 in the VirB5 length band** |
| 4 | ~~VirB2 fails on 20.9%, VirB3 on 9.0%~~ | ~~"a two-tier anchor structure"~~ | **RETRACTED -- see below** |

### Instance 4 is RETRACTED. The two-tier anchor structure STANDS.

The recall that appeared to collapse VirB2 (20.9 -> 3.5), VirB3 (9.0 -> 3.2) and
VirB5 (27.1 -> 3.4) was invalid, and its own precision control showed why:

| anchor | conflict rate | what recall actually pulled in | length |
|---|---|---|---|
| VirB5 | **85.6%** | VirB4/VirB4b 47.1%, VirB3+VirB4 23.3%, VirB1 11.3% | median 238 -> **816** |
| TraN_F | **99.9%** | **VirB10 in 7,484 of 7,504** | median 932 -> **403** |

VirB5 members are 200-288 aa; the "recovered" proteins were ~816 aa, i.e. VirB4
ATPases. Recall was **taking proteins from other anchor families**, not finding
divergent members of its own.

Two design errors, both mine:

1. Representatives were the LONGEST member per sequence cluster
   (`max(v, key=len)`), which systematically selects the most domain-rich and
   most promiscuous query.
2. There was no conflict check, even though this same file already argued that
   recall is the unsafe direction because a false positive inflates completeness.
   The argument was written and then not applied.

So VirB5's 27.1% is NOT a threshold artefact on this evidence, the anchor dilemma
does NOT dissolve, and no VirB5-conditional figure needs recomputing on this
basis. `TraN_F` 97.7% -> 5.0% is withdrawn entirely.

**Instances 1-3 survive, and the difference is instructive: each had a control.**
Instance 1 is a direct measurement on named proteins. Instance 2 is a seed test
with a known-positive expectation. Instance 3 was run with a length band AND a
negative control (60-100 aa band: 0/486). Instance 4 had none.

A corrected recall is implemented with medoid representatives instead of longest,
rejection of any protein already assigned to another anchor family at GA, and a
length-range filter. Its numbers replace instance 4 only if they survive the same
precision control.

### A third failure mode, recorded

**Writing down a safeguard and then not applying it.** The asymmetry between
decontamination (safe) and recall (unsafe) is stated earlier in this very file.
The recall was then run without the precision control that asymmetry implies. The
guard existed as prose, not as code -- which is the distinction A1-A10 exist to
enforce.

## Why it matters beyond bookkeeping

It explains a figure that was puzzling on its own: the named exclusion families
cover 46% of slot records but only **22% of unique proteins**. If GA loses
divergent members at this rate, a family's HMM describes its conserved core and
not its membership, and any coverage figure computed at GA understates the family.

It also means "no HMM family hits this protein" is not evidence that the protein
is novel. That inference was used, and it was wrong three times in this project:
the retracted VirB7 blind spot, the unrecognised VirB5 in the candidate pool
(186 unique) and in the recovered set (190 unique).

## What this does NOT license

phmmer against confirmed members is more sensitive, so it also admits more false
positives. Every use of it here is thresholded at E <= 1e-3 or 1e-5 and paired
with a control:

* decontamination uses it to REMOVE anchors, where a false positive costs a
  candidate and is the safe direction
* recall uses it to ADD anchors, where a false positive inflates completeness --
  so recall must be reported alongside GA, never instead of it

Both directions are implemented in `scripts/lib/decontaminate.py` and the recall
table keeps the GA column (`data/anchors/anchor_second_pass_recall.tsv`).

## Two failure modes recorded alongside

**Reading a conclusion off the tail of a count-sorted list.** The MPF_F readout was
first reported as "the >=3/5 aggregation is uninformative -- the families are MOB
machinery, not T4SS structure". That was the last ten lines of a list sorted by
seed count descending; the head was a clean 14-family core hitting all five seeds.

**A Jaccard over the wrong universe.** MPF_F seed overlap was reported as Jaccard
0.222, suggesting two distinct anchor sets. Computed over all 253 families, of
which 167 hit exactly one seed -- transposases, resistance, rep. Over the T4SS
skeleton the overlap is complete. The statistic was right and the universe was
wrong.
