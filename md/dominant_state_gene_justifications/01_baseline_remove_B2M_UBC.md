# Remove B2M and UBC from the Baseline (housekeeping/normalization) gene set

## The evidence

`Baseline (control)`'s 16 up-genes are meant to be stable, condition-
independent reference genes: `dominant_cell_states_fst_clean_cds.R`
computes `gene_z <- sweep(expr_norm, 2, baseline, "-")`, i.e. every other
state's score is expressed *relative to* the mean of this gene set in that
same sample. If a "baseline" gene is not actually stable, it injects
condition-driven noise directly into the reference level every other
state is measured against -- not just its own score.

Leave-one-out cohesion (each gene's mean correlation with the other 15
baseline genes, across all 3,857 samples):

| gene | cohesion |
|---|---|
| RPS27A | 0.389 |
| RPS18 | 0.377 |
| RPL27 | 0.365 |
| RPL30 | 0.352 |
| EEF1A1 | 0.330 |
| RPL13A | 0.319 |
| PPIA | 0.315 |
| PSMB4 | 0.283 |
| TPT1 | 0.269 |
| GAPDH | 0.269 |
| HPRT1 | 0.219 |
| YWHAZ | 0.143 |
| TBP | 0.131 |
| ACTB | 0.101 |
| **UBC** | **-0.014** |
| **B2M** | **-0.040** |

`B2M` and `UBC` are the only two genes with *negative* cohesion -- they
move independently of, not together with, the other 14 genes meant to
represent a stable reference level.

## Why this isn't a transcription-vs-translation artifact

Both genes have a straightforward, literature-established reason to be
unstable that has nothing to do with translational control specifically:

- **`B2M`** (beta-2-microglobulin) is a core structural component of MHC
  class I and a well-documented interferon-stimulated gene at the
  transcript level. There is no reason to expect its translation to be
  *more* stable than its transcription here -- if anything, a gene this
  central to a complex whose whole point is to respond to interferon
  should be expected to move at both levels.
- **`UBC`** (ubiquitin C, one of the four human ubiquitin-coding loci) is
  itself part of the proteotoxic/heat-shock stress response --
  free-ubiquitin demand increases under many stress conditions this
  compendium profiles.

So this is not a case where GO/literature is "wrong because it's
transcription-focused" -- both genes were never good universal reference
genes in the first place, and this atlas's own translation data
independently confirms that.

## Change made

`scripts/dominant_state_module_definitions.R`,
`base_dominant_state_modules()`, `Baseline (control)` up-gene list:

```diff
- c("ACTB", "GAPDH", "HPRT1", "TBP", "PPIA", "UBC", "YWHAZ",
-   "EEF1A1", "PSMB4", "RPS18", "RPL27", "RPL30", "RPS27A", "TPT1",
-   "RPL13A", "B2M")
+ c("ACTB", "GAPDH", "HPRT1", "TBP", "PPIA", "YWHAZ",
+   "EEF1A1", "PSMB4", "RPS18", "RPL27", "RPL30", "RPS27A", "TPT1",
+   "RPL13A")
```

16 genes -> 14 genes.

## Blast radius

This changes the per-sample normalization reference used by every state's
score, not just `Baseline`'s own reporting -- the entire pipeline
downstream of `clean_cds` needs re-running to pick this up (`fst_pages`,
`clean_cds`, and everything after it in the dependency graph). This is the
highest-impact single change in this batch.

## What would overturn this

If a future re-check shows `B2M`/`UBC` were only "weak" because of a
Interferon/UBC-response subset of samples pulling the aggregate cohesion
down, while the majority of samples show stable co-expression, that would
argue for keeping them and instead handling the Interferon/stress subset
as an outlier-robust weighting problem rather than a removal. Not checked
here -- the leave-one-out number is a compendium-wide average, not
per-subgroup.
