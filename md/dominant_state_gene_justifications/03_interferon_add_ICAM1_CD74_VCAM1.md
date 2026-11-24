# NOT implemented: adding ICAM1, CD74, VCAM1 to Interferon up-genes

**Status: proposed, then reverted the same session after finding
stronger contrary evidence. Left in place as a record of the reasoning
and the self-correction, not as an active recommendation.**

## The original proposal

14 genes are present in the atlas's 129-gene clean-CDS panel but not
assigned to any of the 9 dominant-state modules. Checking each against
all 9 modules' up-gene sets by correlation, `ICAM1` (r=0.45), `CD74`
(r=0.35), and `VCAM1` (r=0.34) correlated with `Interferon`'s existing
up-genes comparably to or better than some of that module's own weaker
existing members, and `Interferon` has only 8 up-genes -- the fewest of
any module besides `Hypoxia`. The initial reasoning treated these three as
a "safer" subset of the 14, distinct from the complement/MHC-II/B-cell
genes (`C1QA/B/C`, `CFB`, `C3`, `HLA-DRA/DRB1`, `CD79A`, `MS4A1`) already
suspected of driving the LCL cell-composition confound, on the grounds
that adhesion molecules and MHC-II invariant chain are less
identity-locked than complement or B-cell lineage markers.

## Why this was reverted

Checking where these 14 genes actually come from (rather than assuming
they were simply orphaned candidates) found
`results/curated_inputs/postviral_fatigue_dominant_state_genes.csv` --
the atlas's own curated fatigue candidate gene table. **All 14 of the
"unassigned" genes are already in this file, explicitly tagged
`Tier = exploratory`** (as opposed to the 43 genes tagged
`Tier = conservative` that actually get scored), split into two named
sub-modules: `complement_endothelial` (`C3, CFB, C1QA, C1QB, C1QC, VWF,
F3, ICAM1, VCAM1`) and `ebv_bcell` (`CD74, HLA-DRA, HLA-DRB1, MS4A1,
CD79A`). Every single row's own `Rationale` column ends with the same
phrase: **"tissue-composition risk."** `ICAM1` and `VCAM1` are tagged
specifically as "Endothelial adhesion marker; tissue-composition risk,"
and `CD74` as "EBV/B-cell context marker; tissue-composition risk" --
i.e. the atlas's own original curators had *already* identified and
deliberately excluded exactly these three genes from scoring, for
exactly the reason (composition confound) this session's correlation
analysis was trying to screen for, and had already grouped them with the
other 11 exploratory-tier genes rather than treating them as a distinct,
safer group.

This means the original proposal's attempt to distinguish "safer" (`ICAM1,
CD74, VCAM1`) from "confound-risk" (the other 11) genes using correlation
strength alone was not just insufficiently supported -- it was answering
a question the atlas's curators had already answered with better
information (presumably direct knowledge of what each gene marks
biologically) than a correlation coefficient can provide. A gene
correlating with `Interferon` is exactly what would be expected whether
that correlation reflects (a) real interferon-pathway co-regulation, or
(b) immune/vascular cell content correlating generally with expression of
immune/vascular genes -- correlation strength alone cannot distinguish
these two explanations, and the curators' domain-specific tagging is the
better available evidence here.

## Change made

**None.** `Interferon / antiviral`'s up-gene list is unchanged from its
original 8 genes (`STAT1, ISG15, MX1, OAS1, IFIT1, IFIT3, IFI44,
CXCL10`).

## What would revisit this

If `ICAM1`/`CD74`/`VCAM1` specifically (not the full 14-gene exploratory
set) are ever re-scored at `conservative` tier and checked against
per-cell-line eta-squared / composition metrics directly (rather than
inferred from raw correlation with `Interferon`), and shown to hold up
across many cell lines rather than being concentrated in
vascular/B-cell-rich contexts, that would be a legitimate reason to
reconsider -- but that check has not been done, and the curators' own
tagging is the more informative prior until it is.
