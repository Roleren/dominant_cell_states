# Remove LDHA and PDK1 from OXPHOS's down-genes (keep them as Hypoxia up-genes)

## The evidence

`OXPHOS (mitochondrial)`'s down-genes are meant to be the Warburg-effect
counter-markers: glycolysis/hypoxia genes expected to *decrease* when
OXPHOS capacity is high. Their correlation with the OXPHOS up-gene set
should therefore be negative.

| gene | role | mean correlation with OXPHOS up-genes |
|---|---|---|
| SLC2A1 | down | -0.036 (weak, correctly signed) |
| **LDHA** | down | **+0.183** (wrong sign) |
| **PDK1** | down | **+0.134** (wrong sign) |

Both `LDHA` and `PDK1` correlate *positively*, not negatively, with
OXPHOS's up-gene set in this compendium's clean-CDS usage data.

Critically, this is not a case of `LDHA`/`PDK1` being bad markers in
general: as `Hypoxia`'s own up-genes, both cohere correctly and
reasonably well with the rest of that module (`LDHA` 0.175, `PDK1` 0.171,
mid-pack among `CA9` 0.250, `BNIP3` 0.246, `VEGFA` 0.170, `EGLN3` 0.166,
`SLC2A1` 0.155, `NDRG1` 0.148). The problem is specific to their *reuse*
as an OXPHOS-down marker, not their quality as hypoxia markers.

## Why this is very plausibly a transcription-vs-translation effect, not an error in the underlying biology

The Warburg-effect logic (glycolysis up, OXPHOS down, under HIF1A) is
well-established **at the mRNA level**. This atlas's own
`dominant_rdg_consensus_candidates.csv` independently ranks `LDHA` as its
#1 consensus-tier gene, with a documented **-13.9% clean-CDS usage change
under ISR/ER-stress conditions** -- a *decrease* in translation despite
HIF-driven transcriptional induction being extremely well established.
The working hypothesis already recorded for that finding
(`dominant_state_open_analysis_2026-09-24.md`) is a **translational
throttle layered on top of real transcriptional induction**: `LDHA` mRNA
plausibly still goes up under hypoxia/stress (literature is not wrong at
the transcript level), but its *ribosome allocation specifically to the
clean CDS*, which is what this atlas measures, does not track that
increase 1:1 -- and may instead track a sample's general permissiveness
for cap-dependent/uORF-gated translation, a variable that also happens to
correlate with how well OXPHOS genes translate. That produces exactly the
observed pattern (positive correlation with OXPHOS at the translation
level, despite being a real, opposite-direction Warburg gene at the
transcript level) without requiring the underlying transcriptional
biology to be wrong.

## Change made

`scripts/dominant_state_module_definitions.R`,
`base_dominant_state_modules()`, `OXPHOS (mitochondrial)` down-gene list:

```diff
- c("SLC2A1", "LDHA", "PDK1")
+ c("SLC2A1")
```

`LDHA` and `PDK1` remain unchanged as `Hypoxia / HIF1A program` up-genes
(no edit needed there -- they were already correctly placed).

## Blast radius

Changes `OXPHOS`'s down-gene list only. `LDHA`/`PDK1` are not removed
from the 129-gene marker panel (they are still required via `Hypoxia`),
so this does not change `fst_pages`' gene-selection set or the total
panel size -- it only changes how `OXPHOS`'s own score is computed
(one fewer down-gene subtracted).

## What would overturn this

A direct comparison of `LDHA`/`PDK1` clean-CDS usage against matched
bulk-RNA-seq (or non-clean-CDS-normalized total footprint count) for the
same samples/studies, which would directly test the "transcript up,
translation-allocation flat/down" hypothesis instead of inferring it from
correlation structure alone. Not done yet -- flagged as the concrete
follow-up in the revision plan.
