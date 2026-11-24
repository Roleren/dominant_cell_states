# Is the LCL outlier signature EBV reprogramming, or just B-cell identity?

Updated: 2026-09-25

## The question

`dominant_state_pca_umap_correlation_patterns.md` showed LCL (lymphoblastoid
cell line) samples forming a completely separate island on PC3 of the
gene-level PCA, isolated from the rest of the 3,857-sample compendium.
LCLs are not an arbitrary cell line: they are **EBV-transformed B cells**
-- Epstein-Barr virus immortalizes primary B cells into LCLs by driving a
continuously active viral gene-expression program (the "latency III"
program), not a one-time hit-and-run transformation event. So the
question is real and answerable with this atlas's own data: is the PC3
outlier signature evidence of ongoing viral reprogramming, or is it just
"these are B cells, and B cells look different"?

**Short answer: the evidence favors real, EBV-driven reprogramming, not
generic B-cell/blood-lineage identity.** Three independent checks all
point the same way, and the specific mechanism has direct, named
literature support.

## Check 1: the signature is not explained by "these are blood/immune cells"

If the PC3 signature were just B-lymphocyte identity, other blood-lineage
cell lines in the compendium should show something similar. They don't.
Comparing mean log2 expression of the 8 canonical interferon-stimulated
genes (ISGs) already in this atlas's `Interferon / antiviral` module
(`STAT1, ISG15, MX1, OAS1, IFIT1, IFIT3, IFI44, CXCL10`) across the three
blood-lineage cell lines in the compendium:

| cell line | lineage | EBV-transformed? | mean IFN-gene score |
|---|---|---|---|
| **LCL** | B lymphocyte | **yes** | **4.74** |
| K562 | myeloid (CML) | no | 0.76 |
| HSPC | hematopoietic stem/progenitor | no | 0.77 |
| whole compendium (all 3,857 samples) | -- | -- | 1.50 |

LCL's score is **~6x higher than the other two blood-lineage lines**
(which are statistically indistinguishable from each other, 0.76 vs 0.77)
and ~3.2x the whole-compendium average. Since K562 and HSPC are equally
"blood cells" but not EBV-transformed, and neither shows any elevation,
generic hematopoietic/immune lineage identity is ruled out as the
explanation. (B-cell identity markers themselves, `MS4A1, CD79A, HLA-DRA,
HLA-DRB1, CD74`, are of course sky-high in LCL specifically -- mean score
8.68 vs 1.99 in HSPC and 0.15 in K562 -- but that's expected and
uninteresting on its own; B-cell markers being high in a B-cell line is
not evidence of anything beyond correct cell identity. The ISG elevation
is the separate, interesting signal.)

Per-gene, every single one of the 8 ISGs is elevated, not just one or two
driving an average:

| gene | LCL mean log2 | whole-compendium mean | delta (log2, i.e. fold-change) |
|---|---|---|---|
| MX1 | 6.26 | 0.82 | **+5.43** (~43x) |
| ISG15 | 7.03 | 2.96 | +4.07 (~17x) |
| OAS1 | 4.80 | 1.05 | +3.75 (~13x) |
| IFIT3 | 4.41 | 1.28 | +3.13 (~9x) |
| IFI44 | 3.92 | 0.85 | +3.07 (~8x) |
| STAT1 | 6.29 | 3.27 | +3.02 (~8x) |
| IFIT1 | 3.30 | 1.34 | +1.96 (~4x) |
| CXCL10 | 1.94 | 0.41 | +1.54 (~3x) |

This is about as complete a canonical ISG signature as could be asked
for.

## Check 2: it's constitutive, not an artifact of an experimental infection condition

185 LCL samples come from 6 studies with genuinely different designs --
critically, most of them **never deliberately infected or interferon-
treated their cells**:

| study | author | n | condition | mean IFN score |
|---|---|---|---|---|
| PRJNA262009 | Battle | 72 | (none -- population-genetics baseline samples) | 4.76 |
| PRJNA275386 | Cenik | 50 | WT | 4.85 |
| PRJNA266341 | Cenik | 50 | infected | 4.82 |
| PRJNA599422 | Ouspenskaia | 7 | mutant | 3.41 |
| PRJNA292112 | Wang | 4 | WT | 4.77 |
| PRJNA281481 | Reid | 2 | WT | 4.32 |

The Battle et al. samples (72 of 185, the single largest contributor) are
a classic population-genetics ribosome-profiling panel -- genetically
diverse LCLs used as a tractable, immortalized proxy system, with no
infection or interferon challenge in the experiment at all. The Cenik
"WT" samples (baseline controls in a separate study) score **4.85 -- if
anything slightly higher than the same study's own "infected" condition
(4.82)**. In other words, deliberately infecting these cells in that one
study added essentially nothing on top of the LCL line's own baseline --
consistent with the baseline already sitting near a ceiling set by
constitutive EBV latency, not by whatever the "infected" treatment was
in that specific experiment. Across every study and every condition
label, the score stays in a tight 4.3-4.85 band (excluding the small,
7-sample Ouspenskaia set at 3.41) regardless of what was actually done to
the cells. That is the signature of a property of the cell line itself,
not of a specific treatment.

Consistent with this, **94.6% of all 185 LCL samples (175/185) call
`Interferon / antiviral` as their calibrated dominant state** -- by far
the most concentrated single-state call rate of any major cell-line
cluster examined this session (compare A549's 73.1% Senescence rate,
itself already called out as an extreme case in
`dominant_state_pca_umap_correlation_patterns.md`).

## Check 3: there is a specific, named literature mechanism, not just a correlation

Epstein-Barr virus's "latency III" program -- the specific viral
gene-expression state that defines LCLs, distinct from the more
restricted latency programs seen in some EBV-associated cancers --
continuously expresses **LMP1** (latent membrane protein 1), a
constitutively active mimic of CD40 signaling that drives **constitutive
NF-kB activation**. This has been directly linked to autocrine interferon
secretion and STAT1 activation in EBV-transformed B cells specifically:
LMP1-induced STAT1 tyrosine phosphorylation was shown to be
"almost exclusively due to the NF-kB-dependent secretion of interferons,"
and EBV latency-III gene products were shown to drive interferon
signaling machinery including IFN-gamma receptor alpha, IRF-4, and STAT1
([Journal of Virology 2005](https://journals.asm.org/doi/full/10.1128/jvi.79.8.4936-4943.2005);
[Journal of Virology 2004](https://journals.asm.org/doi/abs/10.1128/jvi.78.8.4108-4119.2004);
[PMC3717723, NF-kB and IRF7 pathway activation by LMP1](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3717723/)).
This is a precise mechanistic match to what shows up here: LMP1 -> NF-kB
-> autocrine interferon -> STAT1 -> the exact downstream ISG panel
(`ISG15, MX1, OAS1, IFIT1, IFIT3, IFI44, CXCL10`) this atlas's own
Interferon module already tracks.

**A second, translation-specific mechanism exists in the literature but
is more contested, and is worth flagging precisely because this atlas
could help resolve it.** EBV also expresses two extremely abundant,
EBV-encoded small noncoding RNAs, EBER1 and EBER2 (~10^7 copies/cell in
latently infected cells), which bind PKR (the double-stranded-RNA-
activated kinase that phosphorylates eIF2-alpha to shut down translation
during antiviral defense) in vitro
([EMBO Reports, EBER1-PKR binding site](https://link.springer.com/article/10.1093/embo-reports/kvf137)).
However, at least one study found EBERs did **not** measurably change
eIF2-alpha phosphorylation levels between EBV-positive and EBV-negative
cells in living cells, despite the in vitro PKR-binding activity being
real ([Journal of Virology 2005, "Protection from Interferon-Induced
Apoptosis... Is Not Mediated by Inhibition of PKR"](https://journals.asm.org/doi/10.1128/jvi.79.23.14562-14569.2005)).
So the transcriptional/NF-kB-driven ISG story (Check 3's main mechanism)
is well-supported; a direct EBER-PKR-eIF2alpha translational-control
layer on top of it is plausible but not settled in the literature.

## What this atlas specifically adds, and what's still open

Because this atlas's per-gene values are derived from ribosome profiling
(clean-CDS usage / ribosome allocation), not bulk mRNA-seq, **the ISG
elevation documented here is already a translation-level observation**,
not an assumption carried over from transcriptomics literature -- these
genes are more heavily translated in LCLs, which is one level more
specific than "more mRNA is present," though it does not by itself
distinguish "more translation because there is proportionally more mRNA"
from "more translation because of an EBER/PKR-style translational-control
layer on top." The concrete next step to separate these two would be
checking this atlas's own uORF-regulation and translation-efficiency
outputs (`dominant_uorf_regulation_edges.csv`,
`dominant_uorf_structure_rules` outputs) specifically for these 8 ISGs
in LCL samples versus the rest of the compendium, to look for an
LCL-specific translational-control signature (e.g. altered uORF
engagement) beyond a simple proportional increase -- not done in this
pass, and a genuinely novel test this dataset is well-positioned to run,
since the existing literature checked eIF2-alpha phosphorylation directly
rather than genome-wide ribosome allocation.

## Conclusion

Yes -- the evidence here supports treating the LCL PC3/PC4 outlier
signature as real, mechanistically specific evidence of EBV latency-III
viral reprogramming (via LMP1-driven constitutive NF-kB and autocrine
interferon signaling), not an artifact of B-cell lineage identity or of
any specific experimental treatment. All three checks (lineage
specificity against other blood cell lines, constitutiveness across
conditions and studies, and a named, literature-documented mechanism)
point the same direction. The one part of the mechanism that remains
genuinely open -- whether there is an EBER/PKR-style translation-specific
layer beyond the NF-kB/STAT1 transcriptional story -- is exactly the kind
of question this atlas's ribosome-profiling-specific outputs could
address directly, and is worth a dedicated follow-up if this line of
analysis continues.
