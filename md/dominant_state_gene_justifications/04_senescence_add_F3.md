# NOT implemented: adding F3 to Senescence up-genes

**Status: proposed, then reverted the same session after finding
stronger contrary evidence. Left in place as a record of the reasoning
and the self-correction, not as an active recommendation.**

## The original proposal

`F3` (tissue factor) was the best-correlated of the 14 currently-
unassigned genes against `Senescence / SASP` (r=0.41), ahead of its
second-best match (`EMT`, r=0.27), and stronger than several of
`Senescence`'s own existing up-genes. The original reasoning noted `F3`
fits the same biological family as an already-present, validated module
member (`SERPINE1`, also a coagulation-pathway gene), and is not obviously
a cell-identity marker the way complement/MHC-II genes are.

## Why this was reverted

The same check that reverted the `Interferon` additions (see
`03_interferon_add_ICAM1_CD74_VCAM1.md`) applies directly to `F3`: it is
one of the 14 genes already present in
`results/curated_inputs/postviral_fatigue_dominant_state_genes.csv` at
`Tier = exploratory`, under the `complement_endothelial` sub-module, with
`Rationale = "Coagulation marker; tissue-composition risk"`. The atlas's
own curators had already flagged `F3` specifically -- not just the
complement/MHC-II genes -- as a composition-risk gene, for a documented
reason directly analogous to the endothelial/vascular-composition concern
already raised for `ICAM1`/`VCAM1`. This directly undercuts the original
proposal's claim that `F3` was meaningfully different in risk profile
from the rest of the 14-gene exploratory set.

## Change made

**None.** `Senescence / SASP`'s up-gene list is unchanged from its
original 10 genes.

## What would revisit this

The same standard as the `Interferon` reversal: a direct check of `F3`'s
correlation with `Senescence` restricted to samples/cell lines without
substantial vascular/endothelial content, rather than inferring safety
from a whole-compendium correlation number alone.
