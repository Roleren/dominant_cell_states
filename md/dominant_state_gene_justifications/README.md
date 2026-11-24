# Dominant-state gene-set change justifications

Each file here documents one proposed change to
`scripts/dominant_state_module_definitions.R`'s `base_dominant_state_modules()`:
the empirical evidence that motivated it, why it isn't simply explained by
this atlas measuring translation (ribosome occupancy / clean-CDS usage)
rather than transcription, and what would overturn the decision. All
evidence comes from the 129-gene x 3,857-sample Spearman correlation
matrix computed in `../dominant_state_pca_umap_correlation_patterns.md`,
refined in `../dominant_state_gene_set_revision_plan_2026-09-25.md`.

**Two of the four originally-proposed changes were implemented; two were
proposed and then reverted the same session after a follow-up check found
better evidence than the correlation numbers alone.** Both outcomes are
kept here because the reversal is itself informative -- see
`03_interferon_add_ICAM1_CD74_VCAM1.md` and `04_senescence_add_F3.md`.

**End-to-end verification (2026-09-25):** a full `--scope all` pipeline
rebuild completed with 54/55 steps `ok` (the sole error being the
pre-existing, unrelated `manual_translons` gap), total runtime 81.2
minutes (see `../history/dominant_state_pipeline_performance_plan_2026-09-25.md`,
written from timing this same run). The regenerated marker panel has 127
genes (129 - `B2M` - `UBC`), confirmed absent from
`human_dominant_cell_states_clean_cds_gene_expression.csv`. The ISR/
Senescence sub-score script (item 5 below) was run once immediately after
the edit and again after the full rebuild, producing numerically
identical results both times.

| file | change | status |
|---|---|---|
| [01_baseline_remove_B2M_UBC.md](01_baseline_remove_B2M_UBC.md) | Remove `B2M`, `UBC` from `Baseline (control)` up-genes | **Implemented** |
| [02_oxphos_remove_LDHA_PDK1_downgenes.md](02_oxphos_remove_LDHA_PDK1_downgenes.md) | Remove `LDHA`, `PDK1` from `OXPHOS`'s down-genes only (kept as `Hypoxia` up-genes) | **Implemented** |
| [03_interferon_add_ICAM1_CD74_VCAM1.md](03_interferon_add_ICAM1_CD74_VCAM1.md) | Add `ICAM1`, `CD74`, `VCAM1` to `Interferon` up-genes | Proposed, then **reverted** |
| [04_senescence_add_F3.md](04_senescence_add_F3.md) | Add `F3` to `Senescence` up-genes | Proposed, then **reverted** |
| [05_isr_and_senescence_subscore_split.md](05_isr_and_senescence_subscore_split.md) | New standalone sub-score script: `isr_atf4_direct`/`isr_terminal_upr`, `senescence_acute_sasp`/`senescence_commitment` | **Implemented** (standalone, not yet pipeline-wired) |

## Why the two additions were reverted: a correction worth reading in full

The 14 genes already present in the atlas's 129-gene clean-CDS panel but
not assigned to any of the 9 modules were initially assumed to be
orphaned candidates from an earlier, still-blocked "antigen presentation"
idea. Checking their actual source
(`results/curated_inputs/postviral_fatigue_dominant_state_genes.csv`)
found they are not orphaned at all: they are **all 14 of the file's
`Tier = exploratory` genes** (as opposed to the 43 `Tier = conservative`
genes that are actually scored), organized into two named sub-modules
(`complement_endothelial`, `ebv_bcell`), and every one of them carries a
`Rationale` ending in **"tissue-composition risk"** -- including `ICAM1`,
`VCAM1`, `CD74`, and `F3`, the four genes this session's correlation
analysis had initially proposed promoting as the "safer" subset. The
atlas's own prior curators had already reached the opposite conclusion,
with better information (presumably direct biological knowledge of what
each gene marks) than a whole-compendium correlation coefficient can
provide on its own. The additions were reverted on finding this.

**The general lesson, stated plainly**: before proposing a gene-set change
from a correlation number alone, check whether the gene already has
curator annotations, tier assignments, or rationale notes elsewhere in the
codebase -- a correlation cannot distinguish "real pathway co-regulation"
from "cell-composition confound," but a domain expert's documented
reasoning about what a specific gene marks biologically often can, and
takes priority when the two disagree.

## Not proposed as changes, and why

`CDKN2A` stays in `Senescence` (weak cohesion is plausibly a fast-SASP
vs. slow-chromatin-commitment timescale mismatch, not an error -- see the
revision plan); `SNAI1` stays in `EMT` (weak but not wrong-signed);
`HSPA5`/`ASNS` stay in `ISR` (weak cohesion, not wrong-signed -- part of
the sub-score discussion in the revision plan rather than a removal
candidate).
