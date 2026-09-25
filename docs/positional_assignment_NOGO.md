# Positional assignment: step 5 does NOT give a go, for VirB5

`scripts/90_positional_assign.py`, steps 1 and 5 only, on 300 VirB4⁺ plasmids.
Outputs: `data/positional/length_priors.tsv`, `loo_validation.tsv`.

## Step 1 — length priors (this part worked)

Built from **GA-level HMM hits only**, no exclusion family involved.

| VirB | n | median | accept band (IQR ×1.5) |
|---|---|---|---|
| VirB1 | 251 | 215 | 141–305 |
| VirB2 | 241 | 99 | 69–141 |
| VirB3 | 197 | 105 | 87–119 |
| VirB4 | 313 | 852 | 610–1098 |
| **VirB5** | **237** | **237** | **200–288** |
| **VirB6** | **300** | **346** | **289–401** |
| VirB8 | 303 | 232 | 199–263 |
| VirB9 | 297 | 294 | 232–344 |
| VirB10 | 307 | 392 | 326–486 |
| VirB11 | 325 | 338 | 307–367 |

VirB7 excluded as designed — no family, and its band would overlap the 69–76 aa
slot occupants.

**Useful byproduct:** the measured VirB5 band (200–288, median 237) independently
supports the `aa > 150` conjunct chosen from n=4. The cut sits 50 aa below the
observed Q1 across 237 instances.

## Step 5 — the reported 100% is close to tautological

| VirB | known | attempted | correct | accuracy |
|---|---|---|---|---|
| VirB2 | 241 | 36 | 36 | 100% |
| VirB3 | 197 | 80 | 80 | 100% |
| VirB4 | 313 | 67 | 67 | 100% |
| **VirB5** | **237** | **6** | 6 | **100% — 6 cases, meaningless** |
| VirB6 | 300 | 2 | 2 | 2 cases |
| VirB8–VirB10 | ~300 ea | 76–133 | all | 100% |

**Do not read these as validation.** The test only fires when the plasmid already
follows canonical VirB order *and* exactly one gene sits in the gap. Under those
conditions the canonical middle position is the answer by construction, so 100%
mostly confirms that the surviving cases follow canonical order — not that
position predicts identity for an unknown gene.

## The rejection breakdown is the real result

| VirB | known | no_flank | strand | **order** | multi | **count** |
|---|---|---|---|---|---|---|
| VirB4 | 313 | 22 | 2 | **131** | 77 | 14 |
| **VirB5** | 237 | 16 | 1 | 87 | 22 | **105** |
| **VirB6** | 300 | 35 | 13 | 115 | 23 | **112** |
| VirB9 | 297 | 9 | 0 | **130** | 0 | 25 |

**1. `count` dominates exactly at VirB5/VirB6 (105 and 112).** That is the slot
case: `VirB4 | VirB5 | <occupant> | VirB6` puts two genes where the canonical
order expects one, so condition 3 refuses. This is correct and conservative
behaviour — but it means **the four-condition rule structurally cannot assign
anything in the slot region of a plasmid that HAS a slot occupant**, which is the
region of interest. Positional assignment cannot help where it is most wanted.

It also explains VirB5's tiny sample: only 6 of 237 instances are testable,
*because* something usually occupies the adjacent slot.

**2. `order` rejections are large everywhere (70–131 per position) — but part of
that is MY BUG, not biology.**

Condition 2 was implemented as **global** canonical VirB numbering
(`opos = {lab: i for i, lab in enumerate(order)}`, then `ri > li`). Pipeline v2
§4.4 explicitly requires the opposite: *"name slots by the locally observed
adjacent pair (X|Y), never by forcing a global VirB order, or rearranged plasmids
get mis-scored."* I made exactly the error that section warns against, so my
attribution of these rejections to Guglielmini 2013's virB5/virB6-after-virB10
observation was premature. Some unknown fraction is the artefact.

Worse, the two rejection classes are **not independent**: an extra gene in the
gap perturbs the global numbering and can trigger an `order` rejection instead of
a `count` one. So `order + count` cannot be partitioned cleanly, and any
denominator built from them (e.g. "111 order-passing VirB5 cases") is not clean
either. Figures derived from this table are indicative, not measurements.

## Verdict

**NO-GO for steps 3 and 4.** Not because accuracy is low, but because it is
unmeasurable where it matters:

* VirB5 has 6 testable cases. The design's own criterion — "if VirB5 accuracy on
  IncN/IncW is 95%+, extrapolating to IncP is justified" — cannot be evaluated at
  n=6.
* Even with perfect accuracy, condition 3 refuses the slot region by
  construction, so steps 3–4 would add skeleton annotation away from the slot and
  nothing where the project needs it.
* Canonical order fails often enough (`order` rejections) that condition 2 is
  weaker than assumed.

**What was gained:** the length-prior table, which is independently useful and
corroborates the `aa > 150` VirB5 conjunct at n=237 rather than n=4.

**Superseded by a direct measurement.** Rather than infer occupancy from these
contaminated rejection counts, `scripts/91_gap_occupant_lengths.py` measures the
VirB5|VirB6 gap directly — both anchors applied with their own thresholds, gap
contents enumerated, occupant lengths binned. That is the number to use; see
`docs/slot_occupancy.md`.

**Reconsider if** step 1a shows the slot is frequently EMPTY. Empty slots have no
extra gene, so `count` would match and VirB5's testable sample would grow. The
occupancy number therefore decides whether this module is worth revisiting — one
more reason to read it first.
