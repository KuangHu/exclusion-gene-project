# Incidental observations worth following up

## R16a carries an IncP-type `trbK`

`NZ_KX156773.1` (R16a, IncA, 173,094 bp) annotates
`/gene="trbK"`, *"entry exclusion lipoprotein TrbK"*, **78 aa**, `WP_000718002.1`,
at **35,684..35,920** — far from its TraG at 164,442..168,056.

Found only because it false-positived as IncC `eexC`. An IncA/IncC element
carrying an IncP-family entry-exclusion lipoprotein is not something the seed
schema predicts: IncC exclusion is supposed to run through EexC/TraG_C. Either
R16a has acquired a second, IncP-type exclusion system, or the annotation is a
homology-transfer artefact.

Belongs in the **candidates** table (E4), never the seed. Worth checking whether
a `trbL`/VirB6 homolog sits near 35.7 kb — if one does, this is a second
independent instance of the RP4 `trbK`–`trbL` module in a completely different
Inc group.
