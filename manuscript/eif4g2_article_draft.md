> **Working-draft status (remove before submission).** This text reports only
> completed public-data analyses. It does not yet contain the causal reporter,
> endogenous protein or immunopeptidomic experiments required for the proposed
> Nature claim. Square-bracketed author metadata remains to be completed.

# A tRNA-selective switch in an immunogenic EIF4G2 upstream ORF

[Author list]

[Affiliations]

*Correspondence to: [corresponding author and email]*

## Summary

**Transfer RNA abundance can reshape translation according to codon demand,
whereas upstream open reading frames partition scanning ribosomes between short
peptides and downstream proteins. Whether these regulatory layers intersect at
specific human transcripts is largely unknown. Here we identify a selective
redistribution of ribosome occupancy into an ATG-initiated upstream open
reading frame in the leader of EIF4G2, which encodes the stress-responsive
translation factor DAP5. Overexpression of tRNA-Glu(UUC), but not the matched
metastasis-promoting tRNA-Arg(CCG), shifted occupancy from the EIF4G2 coding
sequence into this upstream reading frame without a commensurate change in the
dominant messenger RNA isoform. The response was the strongest
tRNA-Glu-specific outlier among count-supported, non-overlapping leader reading
frames and survived biological-replicate, read-length, frame, coordinate and
peak-deletion challenges. Independently processed ultraviolet-stress data
reproduced the switch in both short and canonical ribosome footprints, and loss
of the recycling factor ABCE1 produced a second recurrence. The same upstream
reading frame has independently been shown to yield an HLA-presented peptide
[1–4]. These observations nominate a tRNA-selective, recycling-sensitive
translation branch with a defined immunogenic output, while reporter and
immunopeptidomic tests are required to establish causality.**

## A discrete translation branch in the EIF4G2 leader

Changes in transfer RNA pools can selectively alter the stability and ribosome
occupancy of codon-biased transcripts [1]. Upstream open reading frames (uORFs)
provide a second layer of selectivity by changing whether scanning ribosomes
translate a leader peptide, reach a downstream start, or reinitiate [5,6]. The
intersection of these layers is potentially consequential in cancer, where
tRNA pools, non-canonical initiation and antigen presentation are all remodelled,
but it is difficult to resolve from conventional gene-level translation
measurements.

We searched a compendium of public human ribosome-profiling data for
condition-dependent redistribution between annotated leader uORFs and a
disjoint, uniquely assigned main coding sequence (CDS). Candidate discovery was
followed by an exact-union validator that collapsed sequencing lanes into true
biological units, required count and reading-frame support, and challenged each
effect by coordinate perturbation and deletion of high-occupancy positions.
This progression was deliberately conservative: 95 review rows collapsed to 60
candidate-direction rows, and none was accepted directly from the discovery
queue.

The strongest surviving candidate was EIF4G2 transcript ENST00000339995
(RefSeq NM_001418). Its 308-nucleotide leader contains two non-overlapping
annotated uORFs. uORF1 is an ATG-initiated, 51-nucleotide reading frame at
transcript positions 23–73; uORF2 is a TTG-initiated, 93-nucleotide reading
frame at positions 128–220. The 2,724-nucleotide main CDS begins at position 309
with the conserved non-AUG GUG start of EIF4G2 [7]. The uORF1 peptide is
`MEVAAGTEWRLQQRLL*`, and an independent study has directly detected an
EIF4G2-uORF-derived HLA ligand and demonstrated leader- and uORF-start-dependent
presentation after mitotic arrest [4]. Thus, the biological output of this
branch is defined independently; the question addressed here is how ribosome
occupancy is routed into it.

## tRNA-Glu selectively shifts EIF4G2 occupancy into uORF1

The discovery contrast derives from an isogenic MDA-MB-231 experiment in which
tRNA-Glu(UUC) was stably overexpressed [1]. The SRA submission split each
biological replicate across three sequencing lanes, with separate accession
identifiers. Collapsing those lanes yielded two tRNA-Glu and two control
biological units. We defined clean-CDS allocation as clean-CDS counts divided by
the sum of clean-CDS and union-uORF counts. The tRNA-Glu samples contained 1,118
union-uORF and 2,653 clean-CDS footprints, compared with 2,033 and 10,720 in
controls. Clean-CDS allocation fell from 0.841 to 0.704, a 0.137 absolute shift
towards the uORF union. Total-coverage-normalized uORF occupancy increased
(log2 ratio 0.748) while clean-CDS occupancy decreased (log2 ratio −0.388).

Both biological units supported the effect. Clean-CDS allocations were 0.700
and 0.706 after tRNA-Glu overexpression, compared with 0.826 and 0.851 in
controls; the replicate ranges did not overlap. Union-uORF reads were strongly
phase-0 enriched in both arms (0.811 and 0.748 of phase-informative reads), and
the shift retained an absolute magnitude of at least 0.120 when both feature
boundaries were moved from −2 to +2 nucleotides. uORF1 carried most of the
leader signal and showed a +0.136 uORF-allocation shift; uORF2 was not used as a
mechanistically resolved feature.

A matched experiment from the same source study and parental cell background
provided a stringent specificity control. tRNA-Arg(CCG) also promoted
metastasis in the original work [1], but its overexpression produced no
corresponding EIF4G2 redistribution: clean-CDS allocation changed by +0.018,
with overlapping biological-unit ranges and coordinate-shift support below the
predefined 0.03 gate. Stable-line generation, cell background, laboratory,
translation-inhibitor exposure and the broad metastatic phenotype are therefore
insufficient to explain the tRNA-Glu result.

We next ranked all non-overlapping leader uORFs shared between the tRNA-Glu and
tRNA-Arg datasets under identical count requirements. EIF4G2 uORF1 had a
tRNA-Glu-minus-tRNA-Arg allocation difference of +0.159 and ranked first at
each of three minimum-count thresholds: first of 35, 27 and 18 eligible
features. Depending on the threshold, it lay 13.3–16.0 unscaled median absolute
deviations above the eligible-background median. These values are descriptive
outlier measures rather than z-scores or inferential probabilities. They show
that the effect is not merely the most appealing member of an unreported large
class: within this matched feature set, EIF4G2 is exceptional.

## The tRNA contrast survives technical and RNA-level alternatives

The tRNA-Glu libraries contain a real technical complication. Lane 3 is
30-nucleotide dominant in both conditions and replicates; lanes 1 and 2 are
28-nucleotide dominant in both arms of replicate 1 and in control replicate 2,
but 30-nucleotide dominant after tRNA-Glu overexpression in replicate 2. We
therefore reconstructed all 12 tRNA-Glu/control and four tRNA-Arg/control raw
libraries with biological insert length retained. P-site offsets were inferred
separately for each run and length against 4,670 external CDS-start windows.
A study-wide pool admitted a length only when it passed the phase gate in every
run of that perturbation; a complementary lane-matched pool admitted a length
only when it passed all four case and control runs from the same technical lane
and required all three tRNA-Glu lanes to contribute. The exact-uORF1 count,
phase, direction and replicate-range gates were fixed before the final lane-3
case and tRNA-Arg target profiles were examined.

In the decisive lane-matched raw pool, tRNA-Glu samples contained 1,028 exact
uORF1, 28 uORF2 and 2,568 clean-CDS P sites, compared with 1,870, 65 and 10,488
in controls. Union clean-CDS allocation was 0.709 versus 0.844 (delta −0.136),
and exact-uORF1 allocation was 0.714 versus 0.849 (delta −0.135). Both true
biological units were separated for the union (0.703 and 0.713 versus 0.829 and
0.856) and exact uORF1 (0.708 and 0.719 versus 0.834 and 0.859). Exact uORF1 was
phase-0 dominant in every raw tRNA-Glu/control run. Each independently matched
lane also reproduced the union shift (−0.121, −0.169 and −0.130).

The matched raw tRNA-Arg pool was nearly null (+0.012), yielding a
tRNA-Glu-minus-tRNA-Arg contrast of −0.147 and passing the predefined
specificity gate. The tRNA-Glu effect also passed using the eight study-wide
shared lengths (21, 24 and 26–31 nt; delta −0.141) and canonical 28–32-nt
footprints (−0.145). The 20–23-nt stratum was sparse and failed its exact-uORF1
phase and measurement gates, so it was not used as evidence. Shifting feature
coordinates by −2 to +2 nucleotides left a minimum absolute effect of 0.116;
deleting any single position left at least 0.100, and greedy deletion of the
three most influential positions left −0.060. The old processed tracks showed
the same direction in every technical lane and after deletion of their large
phase-0 position-35 peak. Thus neither footprint composition, an independently
calibrated P-site choice nor one narrow peak accounts for the tRNA-selective
redistribution.

Parallel total-RNA measurements also argue against transcript abundance or
isoform replacement. For NM_001418, mean total-transcript abundance was 200.0
FPKM in controls and 194.9 after tRNA-Glu overexpression (log2 ratio −0.038).
Its fraction of summed EIF4G2 isoform-level total RNA changed by −0.025, and its
ribosome-footprint isoform fraction changed by −0.028, both far smaller than
the 0.137 allocation shift. Short-read isoform estimates are not definitive,
but the dominant transcript remains stable and dominant at the RNA level.

EIF4G2 uORF1 contains two GAG codons, compared with 72 glutamate codons across
the 908-codon main CDS. This suggested a direct decoding model, but a
transcriptome-wide test did not support it. Among 18 strictly count-supported
features, relative glutamate-codon enrichment correlated neither with the
tRNA-Glu allocation response (Spearman rho = 0.129, two-sided P = 0.610) nor
with the tRNA-Glu-specific response (rho = 0.117, P = 0.645). The matched
tRNA-Arg analysis was also null. Direct action at the two GAG codons therefore
remains a reporter hypothesis, not an inference from codon content.

## Orthogonal stress confirms EIF4G2 branch responsiveness

We used unrelated perturbations of ribosome traffic as orthogonal tests of the
EIF4G2 branch, not as evidence that they share the tRNA mechanism. UV
irradiation in HeLa cells is an especially demanding confirmation because the
source experiment intentionally retained 15–34-nucleotide footprints and
reported UV-induced short footprints and ribosome collisions [2]. In the
compendium tracks, two UV and two control samples showed a −0.306 clean-CDS
allocation shift, but pooled processing could not exclude a read-length-
composition or P-site error.

We therefore reprocessed all four raw libraries while retaining read length.
Adapters and random nucleotides were removed according to the source library
design, inserts were uniquely aligned to GRCh38, and P-site offsets were
inferred independently for every run and length. A read length entered the
decisive analysis only when it passed an external CDS-start phase gate in all
four libraries. The resulting shared set comprised 17, 19–22, 24–31 and 33
nucleotides. A physically impossible offset inferred for 16-nucleotide reads in
one control was rejected before target quantification.

Across all shared, well-phased lengths, clean-CDS allocation was 0.511 after UV
and 0.854 in controls (delta −0.343). Both UV replicates (0.495 and 0.524) were
below both controls (0.809 and 0.870). Crucially, the effect recurred separately
in canonical 28–32-nucleotide footprints (0.568 versus 0.909; delta −0.341) and
short 20–23-nucleotide footprints (0.432 versus 0.706; delta −0.274), with
non-overlapping replicate ranges in each stratum. Phase 0 was the dominant
uORF phase in both arms. Coordinate perturbation left minimum absolute effects
of 0.303 and 0.191 for canonical and short reads, respectively. Greedy deletion
of the three positions that most weakened each signal still left shifts of
−0.278 and −0.238. The UV recurrence is therefore not generated by combining
short and canonical footprints or by a narrow position-specific artifact.

Loss of the ribosome-recycling factor ABCE1 in HCT116 cells provided a second
independent recurrence [3]. Across two targeting and two control units,
clean-CDS allocation shifted by −0.197, with non-overlapping replicate ranges,
phase support and coordinate robustness. uORF occupancy rose strongly (log2
ratio 1.737), whereas clean-CDS occupancy also rose modestly (log2 ratio
0.185). UV and ABCE1 loss therefore converge on relative routing into uORF1 but
do not produce identical downstream CDS behaviour. This distinction argues
against treating uORF1 as a simple binary repressor and does not establish that
tRNA-Glu, UV and ABCE1 act through one molecular intermediate.

## A defined immunogenic output and a falsifiable mechanism

The convergence is notable because EIF4G2 encodes DAP5, a translation factor
that supports main-ORF translation from structured and uORF-rich leaders
[5,6], while its own main ORF initiates at GUG and can respond to ribosome
queuing [7,8]. The 51-nucleotide uORF1 is therefore positioned above an
unusually regulated initiation event. Independent immunopeptidomic and reporter
experiments have already shown that this uORF can generate an HLA-presented
peptide and that its start codon is required for presentation after mitotic
arrest [4]. Our analysis adds a tRNA-selective branch response and two
recycling-related recurrences to that established output.

Several interpretations remain open. Increased mature tRNA-Glu could change
elongation through uORF1, initiation upstream of it, or ribosome traffic near
the downstream GUG. Alternatively, overexpressed tRNA-Glu can generate
regulatory fragments, and tRNA-derived fragments have independent functions in
breast cancer [9]. The null codon-dose result makes an indirect, locus-specific
route at least as plausible as a simple two-codon sensor. Moreover, an
independent later MDA-LM2 dataset showed a weak opposite allocation change and
a condition-linked footprint-length/frame imbalance; the data do not support a
general metastatic-state switch.

The present result is therefore a sharply defined causal target rather than a
finished mechanism. The economical first decision is two small reporter
screens: vector, tRNA-Glu and tRNA-Arg crossed with native EIF4G2-leader WT and
uORF1-start mutant in MDA-MB-231 cells; and untreated versus UV crossed with
the same reporter pair in HeLa cells. The tRNA-by-start interaction determines
whether to invest in the main story. UV validates branch responsiveness and
the assay but cannot rescue a failed tRNA result. Only a positive first screen
justifies GAG-to-GAA and pause-site mutants, mature-tRNA-versus-fragment tests,
endogenous DAP5 measurements and targeted HLA immunopeptidomics. If those
experiments reproduce the branch and connect it to peptide presentation, they
would establish a direct route by which a selected tRNA state exposes a
cryptic immunogenic product. If they do not, the public-data convergence
remains a robust translational phenotype without the proposed molecular or
immune consequence.

## Main references

1. Goodarzi, H. *et al.* Modulated expression of specific tRNAs drives gene
   expression and cancer progression. *Cell* **165**, 1416–1427 (2016).
   https://doi.org/10.1016/j.cell.2016.05.046
2. Wu, C. C.-C., Peterson, A., Zinshteyn, B., Regot, S. & Green, R. Ribosome
   collisions trigger general stress responses to regulate cell fate. *Cell*
   **182**, 404–416.e14 (2020). https://doi.org/10.1016/j.cell.2020.06.006
3. Zhu, X., Zhang, H. & Mendell, J. T. Ribosome recycling by ABCE1 links
   lysosomal function and iron homeostasis to 3′ UTR-directed regulation and
   nonsense-mediated decay. *Cell Reports* **32**, 107895 (2020).
   https://doi.org/10.1016/j.celrep.2020.107895
4. Kowar, A. *et al.* Upstream open reading frame translation enhances
   immunogenic peptide presentation in mitotically arrested cancer cells.
   *Nature Communications* **16**, 8008 (2025).
   https://doi.org/10.1038/s41467-025-63405-2
5. Smirnova, V. V. *et al.* Ribosomal leaky scanning through a translated uORF
   requires eIF4G2. *Nucleic Acids Research* **50**, 1111–1127 (2022).
   https://doi.org/10.1093/nar/gkab1286
6. Weber, R. *et al.* DAP5 enables main ORF translation on mRNAs with
   structured and uORF-containing 5′ leaders. *Nature Communications* **13**,
   7510 (2022). https://doi.org/10.1038/s41467-022-35019-5
7. Takahashi, K. *et al.* Evolutionarily conserved non-AUG translation
   initiation in NAT1/p97/DAP5 (EIF4G2). *Genomics* **85**, 360–371 (2005).
   https://doi.org/10.1016/j.ygeno.2004.11.012
8. Kearse, M. G. *et al.* Ribosome queuing enables non-AUG translation to be
   resistant to multiple protein synthesis inhibitors. *Genes & Development*
   **33**, 871–885 (2019). https://doi.org/10.1101/gad.324715.119
9. Goodarzi, H. *et al.* Endogenous tRNA-derived fragments suppress breast
   cancer progression via YBX1 displacement. *Cell* **161**, 790–802 (2015).
   https://doi.org/10.1016/j.cell.2015.02.053
10. Alard, A. *et al.* Breast cancer cell mesenchymal transition and metastasis
    directed by DAP5/eIF3d-mediated selective mRNA translation. *Cell Reports*
    **42**, 112646 (2023). https://doi.org/10.1016/j.celrep.2023.112646
11. Ingolia, N. T., Ghaemmaghami, S., Newman, J. R. S. & Weissman, J. S.
    Genome-wide analysis in vivo of translation with nucleotide resolution
    using ribosome profiling. *Science* **324**, 218–223 (2009).
    https://doi.org/10.1126/science.1168978
12. Tjeldnes, H. *et al.* ORFik: a comprehensive R toolkit for the analysis of
    translation. *BMC Bioinformatics* **22**, 336 (2021).
    https://doi.org/10.1186/s12859-021-04254-w
13. Jiang, H., Lei, R., Ding, S.-W. & Zhu, S. Skewer: a fast and accurate
    adapter trimmer for next-generation sequencing paired-end reads. *BMC
    Bioinformatics* **15**, 182 (2014).
    https://doi.org/10.1186/1471-2105-15-182
14. Dobin, A. *et al.* STAR: ultrafast universal RNA-seq aligner.
    *Bioinformatics* **29**, 15–21 (2013).
    https://doi.org/10.1093/bioinformatics/bts635
15. Qian, Y. *et al.* ALKBH8-mediated codon-specific translation promotes
    colorectal tumorigenesis. *Nature Communications* **16**, 9075 (2025).
    https://doi.org/10.1038/s41467-025-64144-0
16. Rendleman, J. *et al.* Elongationless start-stop elements are
    stress-resilient translation gates that are more repressive than
    uTranslons. *Nucleic Acids Research* **54**, gkag627 (2026).
    https://doi.org/10.1093/nar/gkag627
