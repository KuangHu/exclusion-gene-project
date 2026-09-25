# Syntenic slot pipeline: method, candidates, and the case for them

Job 26311281 (`scripts/181_slot_pipeline.py`), MPF_T / F / I / FA on PLSDB.
118 slots scored over 22,574 plasmids and 3,669,337 proteins.

---

# 1. The method, step by step

## Step 0 — define the machinery block

Recognised machine genes = the class's CONJScan accessory profiles + the five
generics (`virb4`, `traU`, `t4cp1`, `t4cp2`, `tcpA`) + the 14 anchor families, all
at GA. Consecutive hits separated by **≤3 unrecognised genes** are one block;
a wider gap is a block break, not a slot.

## Step 1 — define a slot, syntenically, without choosing an anchor

    slot = the unrecognised gene(s) between a PAIR of recognised machine genes
    slot key = sorted(profileA, profileB)          -- orientation-normalised

This is what makes the method transferable: it needs to know where the machinery
is, never where the exclusion gene is. `occupancy` is computed over elements where
the slot is **defined** (both flanks present), not over all class members.

Position inside the block is **recorded, never gated**.

## Step 2 — the discriminator

Per slot, over its occupants:

| metric | machine part | IS/transposon scar | random ORF | **EEx slot** |
|---|---|---|---|---|
| occupancy | high | low | low | **high** |
| mean pairwise identity | **high** | **high** | low | **low** |
| length CV | low | low | high | **low** |

Only an exclusion slot is *repeatedly occupied, by a different sequence each time,
at a consistent length*. High identity means one gene; low identity with stable
length means one FUNCTION under diversifying selection.

## Step 3 — sub-GA contamination check

Occupants are re-scanned against every CONJScan profile at **E≤1e-3**, not GA. A
pilin below GA is invisible to the recogniser and would pass as an unannotated
occupant — VirB2/VirB5 are short, hypervariable and membrane-embedded, i.e. the
same layer-2 signature as an EEx. Slots >30% sub-GA are flagged
`machine_contaminated` and kept as negative controls, not deleted.

## Step 4 — rank by distance to the known slots

Occupancy, mean identity and length CV are z-scored within class; the slots that
carry a known EEx define a centroid; every slot is ranked by Euclidean distance to
it. **No hand-chosen weights.** The ranking never sees the EEx labels.

## Step 5 — Rules A/B/C as a filter, not a ranker

`length ≤250 AND (lipobox_strict OR tm_count ≥1)`, assigning
`lipid` / `single_pass` / `multi_pass`. Applied INSIDE a surviving slot.
Recall on the reference panel is **6/7**, not clean — see §4.

---

# 2. Validation

Blind gate, pre-registered: ≥3 of 4 known EEx slots in their class's top 10.

| class | known EEx | slot | rank | percentile |
|---|---|---|---|---|
| MPF_I | ExcA | `traY \| trbA` | **1 of 31** | 3.2% |
| MPF_T | Eex_IncN | `virB6 \| virB8` | **2 of 38** | 5.3% |
| MPF_FA | DUF4467 | `orf15 \| virb4` | **3 of 14** | 21.4% |
| MPF_F | TraS_R100 | `traG \| t4cp1` | **5 of 35** | 14.3% |

**4 of 4, gate passed.** Achieved on occupancy + identity + length CV alone — no
topology, no length prior, no positional term.

Robust to the pre-filter: at `MIN_DEFINED` = 10 / 30 / 60 the ranks are 5/6/6/5,
4/1/2/3 and 3/1/2/1 respectively. All stay in the top 10 at every setting; the
"top 4" headline was a `MIN_DEFINED=30` artefact.

---

# 3. The candidates

Ranked by distance to the known-EEx centroid. None carries any hit from the six
exclusion models, or any CONJScan profile even at E≤1e-3.

| class | rank | slot | occ% | mean id | med aa | len CV | dist |
|---|---|---|---|---|---|---|---|
| MPF_F | 1 | `traF \| trbC` | 94.7 | 16.9 | 154 | 0.342 | **0.316** |
| MPF_I | 2 | `traE \| traI` | 99.5 | 25.6 | 194 | 0.468 | 0.352 |
| MPF_I | 3 | `traK \| trbA` | 100.0 | 19.1 | 144 | 0.372 | 0.359 |
| MPF_FA | 1 | `orf17b \| virb4` | 81.0 | 28.8 | 160 | 0.508 | 0.430 |
| MPF_F | 2 | `traL \| tfc7` | 95.6 | 20.9 | 133 | 0.501 | 0.508 |
| MPF_F | 3 | `traB \| traE` | 75.0 | 15.7 | 173 | 0.416 | 0.524 |

For comparison, the known EEx slots span occupancy 0.9–99.4%, mean identity
12.6–41.2%, length CV 0.19–1.32. **The candidates sit inside that envelope on all
three axes** — which is what "distance to centroid" is measuring.

## Architecture

Every candidate is **interior to the machinery block**, flanked by recognised
machine genes on both sides. That matches the knowns: measured on the prototypes,
TrbK sits between virB5 and virB6 (gap 1), ExcA between trbA and traY (gap 1–2),
and only DUF4467 is at a block edge. The earlier claim that all four sit at block
edges was wrong.

**One flank recurs across each class's known slots**, and this is the real
positional signal:

```
MPF_T    virB6  in 6 of 8 known-EEx slots
MPF_F    traG   in 2 of 2
MPF_FA   orf15  in 2 of 2
```

So the biological position is "adjacent to virB6" / "adjacent to traG", and the
pair-key **fragments one position into several slot IDs** depending on which other
machine gene was recognised in that element. See §4.

---

# 4. Why these are plausible entry-exclusion genes

Six independent properties, none of which was used to select them except where noted:

1. **Conserved position, in the machine.** Present at a fixed syntenic position in
   75–100% of elements that define it, flanked by conjugation genes on both sides.
   Entry exclusion acts at the mating pair, so it must be delivered with the
   machine; a gene that travels with the T4SS in most members of a class is doing
   something for the T4SS.
2. **The occupant is not conserved.** Mean pairwise identity 15.7–28.8% among
   occupants of the SAME position. A structural component cannot do this — it has
   to keep interacting with fixed partners.
3. **Length is constrained anyway.** CV 0.34–0.51 despite that divergence. Free
   sequence with constrained size is the signature of a maintained function, not
   of decay or of a mobile-element scar.
4. **This combination is what the four known EEx look like**, and it is the ONLY
   thing used to rank them. The gate shows the combination recovers all four
   knowns blind, so the candidates are being scored by a criterion demonstrated
   to work on the positives.
5. **Not the machine, by two independent tests.** No CONJScan profile hits them at
   GA (by construction) *or* at E≤1e-3 (the sub-GA check, added after that check
   caught a real machine slot at 52% TraT).
6. **Mechanistically coherent size.** 133–194 aa median, against known EEx at
   69/126/147/159/173/220 aa.

**Why this is exclusion rather than some other accessory function** is the weakest
link and should be stated as such. Properties 1–3 identify a *hypervariable
position inside the conjugation machinery under length constraint*. Entry
exclusion is the best-characterised function with that profile — an arms race over
recognition specificity, which is exactly what drives sequence divergence at a
fixed functional site. But nothing here excludes other specificity-driven
functions at the same position: surface exclusion, a partner-recognition subunit,
or an immunity function unrelated to exclusion.

---

# 5. What would raise or kill these

| test | decides |
|---|---|
| Rules A/B/C inside each slot | membrane anchoring — required by the mechanism |
| 3Di/Foldseek across occupants | do non-homologous occupants share a FOLD? convergent architecture is the Rule B signature |
| occupant vs partner allele | does the occupant covary with the machine's specificity determinant? |
| experimental | mating assay; nothing short of it makes a candidate a finding |

---

# 6. Known limitations

**Slot fragmentation.** The pair-key splits one biological position into several
slot IDs. MPF_T's Eex_IncN occupies 8 slots, 6 of which flank virB6. This
fragments the reference set (occupancy 0.9% to 99.4% among "known" slots) and
almost certainly fragments candidates the same way. The fix is to key a slot by a
single recurring flank rather than a pair, and it is not implemented.

**Blindness is partial.** The RANKING is blind. The SLOT CONSTRUCTION is not —
which profiles count as machine, the 3-gene break rule, `MIN_DEFINED=30` were all
chosen by someone who already knew where the four EEx sit. On MPF_G/B/C there is
nothing to tune against, so any implicit tuning here will not transfer. **If
G/B/C rank quality drops sharply, suspect this first, not the absence of EEx.**

**MPF_FA's rank 3 of 14 is the 21st percentile** and near-uninformative on its own.
The 4/4 is one strong result (MPF_I, 1 of 31), two solid, and one weak.

**"Not machine" is only as good as CONJScan.** Sub-GA clean means no CONJScan
profile reaches them, not that they are not structural. A component with no model
at all is indistinguishable from a novel EEx slot here.

**Discovery bias.** An EEx sitting outside the block is invisible to a method built
on in-block synteny. All four knowns are in-block, but they were found by cloning
the tra region in the 1980s.

**Single-occupant slots are invisible.** Step 2 needs several DIFFERENT occupants.
A conserved EEx not under arms-race selection would score as a machine component.

**Evidence tier: E4.** No experimental support for any candidate. The ranking is
a prediction from three summary statistics.
