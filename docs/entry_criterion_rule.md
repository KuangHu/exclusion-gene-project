# Entry criterion selection: the rule, and the threshold that does not transfer

## The rule

Rank candidates by **cross-class contamination**, not by MOB-suite purity.
MOB-suite purity is a reference column only.

    exclude   cross-class admissions (plasmids already admitted by another
              class's entry criterion) above ~1% of admissions
    prefer    highest recall on an independent known-positive set, where one
              exists
    prefer    the VirB4-family ATPase, which has been the most conserved
              component in every class measured
    never     an anchor that also defines the slot -- that makes every admitted
              plasmid slot-ready by construction

## Why the >=90% MOB-suite purity threshold is class-dependent

It was correct for MPF_F and wrong for MPF_I, and the difference is what the
impurity is made of.

**MPF_F**: the four excluded families were eating 6,000-9,000 plasmids of a
DIFFERENT CLASS -- TrbI 6,344 MPF_T, TrwB_AAD_bind 8,885 MPF_T. Real cross-class
contamination. Excluding them was right.

**MPF_I**: I_traU scored 88.7% and would have been rejected by the same
threshold. Its impurity breaks down as:

| | n | % of impurity |
|---|---|---|
| MOB-suite untyped `-` | 404 | **96.0%** |
| MPF_F | 12 | 2.9% |
| MPF_T | 5 | 1.2% |

Cross-class contamination is **17 of 3,732 = 0.46%**, the lowest of any
candidate. The 88.7% was measuring MOB-suite's failure to TYPE, not
contamination. Transferring the threshold would have rejected the VirB4-family
ATPase for a reason that does not apply.

**This will matter again for MPF_FA**, where MOB-suite coverage of Gram-positive
elements is likely worse still. A purity threshold there would reject nearly
everything for the same wrong reason.

## The independent positive set, which MPF_T and MPF_F never had

The Eex census fixed a target: NF033891/ExcA on 2,294 plasmids, of which 2,266
fall outside both existing entry criteria. MOB-suite calls 2,134 of those MPF_I.

That allowed a recall test no purity figure could substitute for:

| group | admits | ExcA recall |
|---|---|---|
| narrow (traT/traV/traE) | ~1,850 | **77.8-80.8%** |
| broad (traY/trbA/trbB/traW/traU...) | ~3,500-3,900 | **96.6-99.8%** |

The 2,022 plasmids broad admits and narrow does not are 82.4% MOB-suite MPF_I,
1,384 of them carry 14 of 17 T4SS_I_* profiles, and 391 are ExcA-positive. Those
are complete systems. **Narrow's 99.0% purity was specificity bought by
discarding real members** -- the exact failure the recall test exists to catch.

## MPF_I decision

**`T4SS_I_traU`** -- 3,732 admitted, 96.6% ExcA recall, cross-class 17 (0.46%),
and the VirB4-family ATPase, preserving convergence with MPF_T VirB4, MPF_F TraC
and MPF_FA ConE.

`I_traY` recalls higher (99.8%) and is REJECTED: TraY is the slot anchor, so
entering on it would make every admitted plasmid slot-ready by construction --
the circularity avoided when MPF_T entered on VirB4 rather than VirB5.

`I_traR` is excluded outright: 179 cross-class admissions (4.8%), of which 166
MPF_T. That is the only genuine contamination among the seventeen.
