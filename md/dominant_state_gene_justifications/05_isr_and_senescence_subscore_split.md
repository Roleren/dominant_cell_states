# Implemented: ISR and Senescence sub-score split

**Status: implemented as a standalone, additive script**
(`scripts/dominant_state_isr_senescence_subscores.R`), **not yet wired
into the main pipeline's dependency graph.** It reads the same per-sample
gene expression matrix and Baseline-relative normalization the main
`clean_cds` scoring uses, and does not modify
`dominant_state_module_definitions.R`'s ISR or Senescence module
definitions -- both states keep their original, single combined score;
this adds two extra sub-score pairs alongside them.

## What it computes

- `isr_atf4_direct` = mean(`ATF4`, `ATF3`, `ASNS`) relative to baseline
- `isr_terminal_upr` = mean(`DDIT3`, `PPP1R15A`, `DDIT4`, `HERPUD1`,
  `HSPA5`, `XBP1`, `CEBPB`) relative to baseline
- `senescence_acute_sasp` = mean(`SERPINE1`, `IGFBP7`, `TGFBI`, `TIMP1`,
  `CCND1`, `IL6`, `CDKN1A`, `CXCL8`) relative to baseline
- `senescence_commitment` = mean(`CDKN2A`) minus `LMNB1`, relative to
  baseline (mirrors the existing module's up/down convention: `CDKN2A` up,
  `LMNB1` down, both real senescence-commitment markers)

Output: `results/dominant_state_isr_senescence_subscores.csv` (one row per
`Run`) plus `..._subscores_summary.csv` (means, SDs, and the
cross-sub-score correlation for each pair).

## Result (confirmed identical before and after the full pipeline rebuild)

| pair | correlation |
|---|---|
| `isr_atf4_direct` vs. `isr_terminal_upr` | **r = 0.58** |
| `senescence_acute_sasp` vs. `senescence_commitment` | **r = 0.066** |

The ISR sub-scores correlate moderately (0.58) -- consistent with the
revision plan's characterization that ISR is a real internal gradient,
not a hard binary split (all 10 original ISR genes showed positive
leave-one-out cohesion with each other; this confirms that holds at the
sub-score level too, just with room to diverge under specific conditions
as proposed). The Senescence sub-scores are **essentially uncorrelated**
(0.066) -- a much cleaner and stronger confirmation of the fast-SASP
vs. slow-chromatin-commitment timescale-mismatch hypothesis than the
single-gene cohesion number alone suggested (`CDKN2A`'s cohesion with the
SASP core was -0.124; the sub-score-level correlation across all 3,857
samples shows the two programs are close to statistically independent
compendium-wide).

## Verification

Run once immediately after the `Baseline`/`OXPHOS` gene-set edits (against
a gene expression matrix that still had the now-removed baseline genes
present as extra, unused rows), then re-run after the full pipeline
rebuild completed and regenerated
`human_dominant_cell_states_clean_cds_gene_expression.csv` with the
corrected 127-gene panel (confirmed via
`dominant_pipeline_last_run_summary.csv`: 54/55 steps `ok`, the sole
error being the pre-existing, unrelated `manual_translons` gap). **Both
runs produced numerically identical sub-score correlations** (ISR: r =
0.581 both times; Senescence: r = 0.066 both times) -- expected, since
none of the genes this script reads (`ATF4`, `ATF3`, `ASNS`, `DDIT3`,
`PPP1R15A`, `DDIT4`, `HERPUD1`, `HSPA5`, `XBP1`, `CEBPB`, `SERPINE1`,
`IGFBP7`, `TGFBI`, `TIMP1`, `CCND1`, `IL6`, `CDKN1A`, `CXCL8`, `CDKN2A`,
`LMNB1`) were touched by the panel edit, and the corrected 14-gene
`Baseline` set was already in effect for both runs.

## Not yet done

Testing the sub-scores' actual value: checking whether
`isr_atf4_direct`/`isr_terminal_upr` diverge specifically under prolonged
vs. acute ER stress within the same study (the cleanest direct test of the
translation-efficiency-switch hypothesis for `ATF4`), and whether
`senescence_commitment` specifically tracks replicative/oncogene-induced
senescence conditions while `senescence_acute_sasp` tracks broader acute
stress conditions (the direct test proposed alongside the new
metadata-field recommendation in the revision plan). Both require the
condition-metadata cross-referencing this note's script does not yet do.
