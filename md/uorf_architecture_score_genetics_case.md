# A genetics/bioinformatics case: uORF architecture density as a generalizable, genetically-relevant score

Updated: 2026-09-24

## Framing

`uorf_to_cds_effects_and_distance_patterns.md` found, from this atlas's own
82-gene marker panel, that cross-gene variation in clean-CDS output is
predicted almost entirely by uORF *architecture density* -- count, overlap
status, and leader-occupancy fraction -- and essentially not by raw
distance from the last uORF to the CDS. This note asks the harder question:
is that a real, generalizable, previously-uncharacterized pattern, or is it
already known? And if it is a real gap, what is the specific, falsifiable,
purely-computational study (no new wet-lab work required) that would turn
it into a genetics/bioinformatics paper rather than a mechanistic one like
EIF4G2 -- which is the explicit brief here: find something in this
compendium's own genetic and genomic reach, not something that needs a new
reporter assay.

The short version: **the pattern is a real, non-trivial reconciliation of
two literatures that have not been directly cross-compared this way, and it
makes a specific, testable prediction about human genetic variant
interpretation that nobody appears to have tested yet.** The full case, the
literature that both supports and complicates it, and exactly what a
genome-wide follow-up would need to do, are below.

## What is already established (so the paper does not reinvent it)

Three independent literatures already converge on parts of this:

**Population genetics: overlapping and count-increasing uORF variants are
under strong purifying selection.** Whiffin et al. (*Nature Communications*
2019, and a related 2021 companion paper) used gnomAD (15,708+ genomes) to
show that variants creating or disrupting uORFs are depleted relative to
neutral expectation -- i.e., under purifying selection -- and that variants
creating a *new overlapping uORF* (one that extends past the annotated
start codon) are under selection comparable to protein-coding missense
variants, the strongest class they tested. A parallel comparative-genomics
study (Ye et al., *PMC3486843*) found overlapping uORFs are genome-wide
under-represented relative to the null expectation in human, mouse, and
fly alike -- the same "overlap is the worst case" signal, independently,
from evolutionary divergence rather than population variation.

**Reporter-based MPRAs: uORF number, distance, overlap, and Kozak context
all matter when experimentally varied.** Johnstone, Bazzini & Giraldez
(*EMBO J* 2016) and the same group's follow-up comparative study
(*Nat Commun* 2016, "Conservation of uORF repressiveness and sequence
features in mouse, human and zebrafish") used massively parallel reporter
assays to show that uORF count, intercistronic distance, CDS overlap, and
initiation context *jointly* predict repression when each is systematically
varied on a synthetic construct.

**Functional genomics: uORF translation itself, not distance, is what
CRISPR screens find essential.** A very recent preprint (2025,
*bioRxiv* 10.1101/2025.07.07.663426) combined ribosome-profiling-validated
uORFs (978 confirmed translated, from 13 published Ribo-seq datasets) with
CRISPR-Cas9 loss-of-function screens across three cell lines. Only 16% of
tested uORFs were essential for proliferation, essentiality was highly
cell-type-specific (6 of 155 universal across all three lines), and --
directly consistent with this atlas's own architecture finding --
frameshift indels that *extend a uORF into an overlap with the CDS* were
the change most consistently linked to translational repression. This
paper intersected its uORF set with ClinVar/dbSNP for disease variants, but
did not connect to gnomAD population constraint scoring, and did not build
a continuous, cross-condition redirection-strength score.

## The gap: nobody has connected compendium-scale, ribosome-occupancy-validated architecture density to population constraint as a continuous score

Each literature above uses a different currency. Whiffin-style constraint
work uses **sequence-predicted** uORFs (annotated computationally, not
necessarily confirmed to be translated) scored against **population allele
frequency**. The 2025 essentiality preprint uses **ribosome-profiling-confirmed**
uORFs scored against **binary CRISPR essentiality** in three cell lines.
The MPRA literature uses **synthetic reporters** scored against
**engineered variation** in one feature at a time. None of them build what
this atlas already computes for free: a **continuous architecture-density
score, measured from real ribosome occupancy across thousands of natural
human samples and dozens of biological conditions**, that could be used to
ask whether the population-genetic constraint on a uORF-disrupting variant
scales with how functionally active that uORF's architecture actually is in
real cells -- not just whether a uORF is predicted to exist.

That is a specific, checkable claim: **a uORF sitting in a dense,
overlapping, compendium-confirmed-active architecture should show stronger
purifying selection at its disrupting/creating variants than a uORF that is
sequence-predicted but rarely used, occupies little of the leader, and does
not overlap the CDS** -- even after controlling for the presence/absence
signal that existing gnomAD-based work already captures. If true, this
would mean current uORF-variant pathogenicity classifiers (e.g.,
UTRannotator, used clinically to flag uORF-creating variants) are
under-using a signal -- real translational architecture density -- that is
sitting unused in public Ribo-seq data and could sharpen variant
interpretation without any new sequencing.

## Why the distance finding specifically matters for this claim, not just as a footnote

The reconciliation from the companion note is directly load-bearing here,
not just a curiosity. The MPRA literature that finds distance matters does
so by *engineering* distance while holding a uORF's other features fixed --
a controlled experiment. Real human genetic variants overwhelmingly do the
opposite: a single-nucleotide variant that creates a new upstream stop
codon, new start codon, or extends an existing uORF changes *architecture*
(count, overlap status, effective uORF length) far more often than it
cleanly changes *distance alone*. If distance is not an independent,
generalizable predictor of output once architecture is accounted for across
real endogenous genes (as this atlas's 82-gene panel shows), then **variant
interpretation tools that lean on positional/distance heuristics are
leaning on the wrong axis for the kind of variant most patients actually
carry**, and should weight overlap-creation and density-increase more
heavily than fine positional detail. This is the concrete, actionable
version of "novel generic pattern": not just an academic reconciliation, but
a specific recommendation for how uORF-variant pathogenicity scoring should
be re-weighted, backed by both this atlas's regression and the independent
CRISPR-essentiality preprint's convergent "overlap-creating indels are what
repress" finding.

## What would make this a real paper, stated honestly

This atlas's current panel is 82 genes, curated for an unrelated purpose
(dominant cell-state marker genes), not a genome-wide uORF survey. That is
enough for a strong pilot signal and a precise hypothesis; it is not enough
statistical power for a standalone genetics paper. The concrete, scoped
program to get there, using infrastructure this atlas already has and
without requiring new wet-lab work:

1. **Score genome-wide, not just the 82-gene panel.** This atlas's own
   `dominant_manual_translons.R` / translon categorization pipeline and
   `coverageByTranscriptFST()` machinery already do this per-gene; the
   blocker is the same one flagged in
   `dominant_cell_state_iteration_review_2026-09-24.md` -- the marker panel
   needs the lab's improved translon predictor before overlapping uORFs in
   particular can be called reliably at scale (currently a hand-curated
   bridge list of 27 features across 18 genes). This is the single
   rate-limiting step for turning this from a pilot into a genome-scale
   result.
2. **Compute the same architecture-density score (leader-occupancy
   fraction, uORF count, overlap status, inter-uORF spacing) for every
   gene with a confidently-called 5' UTR**, using the existing
   3,857-sample compendium this atlas already draws on -- no new
   sequencing.
3. **Join that score, gene-by-gene, to public gnomAD constraint metrics for
   each gene's uORF-disrupting and uORF-creating variant classes**
   (Whiffin/UTRannotator-style variant calls are already public and
   downloadable; this is a data-integration step, not new lab work).
4. **Test whether architecture density predicts constraint strength beyond
   presence/absence**, i.e. whether a uORF's compendium-measured activity
   level adds predictive power over the existing binary "has a uORF"
   signal already used in variant classifiers. This is the single
   regression that either supports or kills the central claim, and it
   is directly falsifiable with data that already exists publicly.
5. **Separately confirm the distance-reconciliation point genome-wide**:
   refit the distance-to-CDS effect within the leader-only-uORF class and
   the overlap class separately (this atlas's 82-gene panel does not yet
   have enough overlap-class genes, 17, to do this with real power) to
   check whether distance re-emerges as significant once it is measured
   from the correct, class-appropriate reference coordinate.

## What would falsify this, stated in advance

If step 4 above finds that architecture density adds *no* predictive power
over simple uORF presence/absence for population constraint, the central
genetics claim is wrong and the finding reduces to a translation-biology
observation only (still worth the companion note, not worth a genetics
paper). If step 5 finds that distance remains non-significant even within
class, restricted to the correct reference coordinate, and at genome scale,
that would upgrade "architecture density, not distance, is the generic
lever" from a pilot-scale observation to a much stronger, more surprising
claim worth leading a paper with on its own.

## Bottom line

This is not yet a finished discovery; it is a sharp, falsifiable hypothesis
sitting at a real gap between three literatures that have not been directly
cross-compared this way, backed by a real (if small) pilot signal from this
atlas's own compendium, and testable entirely with public data plus this
atlas's existing, already-built compendium-scoring infrastructure --
exactly the kind of genetics/bioinformatics result that does not require
new lab reagents, in contrast to EIF4G2's causal-mechanism story, which
does.
