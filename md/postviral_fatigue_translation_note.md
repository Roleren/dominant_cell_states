# Post-viral fatigue, viral infection, and ribosome decision reprogramming

Last updated: 2026-05-21

## Status update: 2026-09-24

The "Immediate implementation plan" below (steps 1-7) has been substantially
executed since this note was first written. The dominant-state atlas now has a
real `Post-viral fatigue / ribosome stress` module, a `fatigue_summary`
pipeline step (`postviral_fatigue_pipeline_summary.R`), and iORF/NTT
diagnostics, exactly as proposed. Concrete quantitative results now exist:
DDIT3, PPP1R15A, ATF4, and IFIH1 are current viral-infection clean-CDS
candidates (+87%, +128%, +81%, +38% respectively), and a therapeutic-
perturbation layer ranks "ISR/eIF2B-ATF4 attenuation" as the top hypothesis
with quantified half-normalization targets. See
`dominant_cell_state_scientific_note.md` for the full current results and
`dominant_cell_state_iteration_review_2026-09-24.md` for what is missing next.

The single biggest remaining gap is that, unlike EIF4G2, this story still has
**no lab go/no-go plan**. `long_covid_me_cfs_medications_overview.md`
confirms no existing approved therapy targets the ISR/ATF4/uORF axis directly,
which makes this a genuine white-space therapeutic hypothesis rather than a
crowded one -- and therefore a strong candidate for the same cheap-screen-first
discipline used in `eif4g2_lab_go_no_go_plan.md`. The rest of this document
(literature summary, gene modules, practical model tests) remains a useful
historical record of the original hypothesis and is not superseded, but its
"Immediate implementation plan" and "Bottom line" sections below should be read
as largely completed rather than pending.

## Working hypothesis

Viral infection may leave a persistent translational state where host ribosomes are not simply "globally lower", but are rerouted across transcript architecture. The most interesting version of the hypothesis is:

1. Acute infection activates host shutoff, interferon, PKR/eIF2-alpha, ER stress, stress granules, and mTOR/TOP remodeling.
2. Those programs change start-codon choice and reinitiation behavior.
3. For genes with upstream overlapping ORFs (uoORFs/ouORFs), dense uORF leaders, N-terminal extensions/truncations, or internal ORFs (iORFs), the measured CDS output can shift without a matching RNA abundance change.
4. In a chronic/post-viral state, this can create a durable loss of normal protein-output phenotype, especially in energy metabolism, immune regulation, endothelial/coagulation biology, and stress-response genes.

This is plausible enough to justify a new analysis track. It is not yet a proven biological mechanism for Long COVID or ME/CFS. The most defensible next step is to build a "post-viral fatigue / ribosome programming" state score and test whether infection-associated metadata branches show systematic uORF, uoORF, NTT/NTE, and iORF rerouting.

The important distinction is that this is not just an "interferon-high" score. The proposed state should capture a sample where antiviral/ISR/inflammatory pressure is high while mitochondrial and normal host translation capacity are low or rerouted. That makes it closer to a ribosome-decision phenotype than to a standard transcriptomic inflammation module.

## Role in the final RDG atlas paper

This module should be the main biological test case for the final paper. The paper should first introduce the general RDG atlas and CDS-buffering framework across known dominant-state effects such as ISR/ER stress, interferon, hypoxia, OXPHOS, mTOR/TOP, proliferation, EMT, and senescence. Those states establish that the atlas can turn thousands of Ribo-seq samples into interpretable ribosome-routing and protein-output proxy tables.

The second main point should then use `Post-viral fatigue / ribosome stress` to test something more novel: whether viral and post-viral contexts induce a host ribosome-decision program where overlapping uORFs, dense leader uORFs, NTE/NTT nodes, and iORFs redirect translation away from clean canonical CDS. The biological claim should be treated as a falsifiable hypothesis. Strong evidence would require matched-control shifts, multiple replicated designs, and feature-level support showing that clean-CDS loss or buffering is explained by increased ouORF/iORF-like routing rather than by low coverage, isoform mismatch, or generic interferon expression.

## Literature summary

### Viral infection strongly remodels host translation

Viruses commonly redirect translation by attacking initiation, mRNA availability, ribosome access, stress granules, and innate immune sensors. SARS-CoV-2 is a clear example: nsp1 binds the 40S ribosomal mRNA entry channel and suppresses host mRNA translation while viral RNAs partly evade shutoff. Other viruses use host mRNA degradation, eIF4F remodeling, PKR antagonism, IRES-like entry, leaky scanning, reinitiation, ribosome frameshifting, and stress-granule manipulation.

Important point for this project: these mechanisms are exactly the mechanisms that should change uORF and iORF usage. If initiation is limiting, scanning is altered, eIF2-alpha is phosphorylated, or eIF1/eIF5 stringency changes, then AUG and near-cognate start choice can move. That can turn a transcript from "canonical CDS translated" into "leader ORF, overlapping ORF, N-terminal truncation, or internal ORF translated".

Relevant sources:

- Viral translation is a host-virus "tug of war" over the translation machinery: https://pubmed.ncbi.nlm.nih.gov/36334591/
- Ribosome profiling exposes viral uORFs, short ORFs, and hidden coding capacity: https://www.annualreviews.org/doi/10.1146/annurev-virology-100114-054854
- Viral short ORFs are common and can be functional: https://pmc.ncbi.nlm.nih.gov/articles/PMC7167739/
- HCMV ribosome profiling revealed extensive unannotated viral ORFs and alternative transcript/start usage: https://doi.org/10.1126/science.1227919
- A 2025 pan-viral massively parallel ribosome profiling study identified thousands of unannotated viral ORFs across hundreds of human-associated viruses and found many viral uORFs likely regulate downstream viral protein synthesis: https://pubmed.ncbi.nlm.nih.gov/40504907/
- SARS-CoV-2 nsp1 blocks host translation by occupying the ribosomal mRNA entry channel: https://www.nature.com/articles/s41594-020-0511-8
- SARS-CoV-2 nsp1 host shutoff and viral escape mechanisms are reviewed here: https://pmc.ncbi.nlm.nih.gov/articles/PMC8481097/

The 2025 pan-viral ORF result is especially relevant for our model. It argues that uORF-like regulatory architecture is not a rare curiosity; viruses broadly use short and upstream ORFs as part of translation control. A host cell fighting infection is therefore a setting where both sides are under strong selection around initiation, scanning, and ORF choice.

### ISR/eIF2-alpha is the bridge from viral infection to uORF biology

The integrated stress response is a central antiviral and stress pathway. Viral RNA activates PKR; ER stress activates PERK; nutrient and mitochondrial stress can activate GCN2/HRI-like branches. These converge on eIF2-alpha phosphorylation, lowering general initiation but selectively increasing translation of specific uORF-controlled transcripts such as ATF4, ATF5, DDIT3/CHOP, PPP1R15A/GADD34, and related feedback genes.

This is directly connected to our RDG work. ATF4 is not just a marker gene; it is the canonical proof that a stress state can alter ribosome flow through uORFs to change CDS output. If post-viral fatigue contains a chronic ISR-like substate, the most sensitive signal may be altered feature usage rather than total CDS expression.

Relevant sources:

- ISR translation inhibition and ATF4/ATF5 delayed reinitiation: https://pmc.ncbi.nlm.nih.gov/articles/PMC5720466/
- uORFs differentially regulate gene-specific translation during ISR: https://pmc.ncbi.nlm.nih.gov/articles/PMC5016099/
- eIF2-alpha kinases, including PKR, connect viral infection to translational control: https://pmc.ncbi.nlm.nih.gov/articles/PMC6028073/
- PKR activation and viral countermeasures: https://pmc.ncbi.nlm.nih.gov/articles/PMC3185532/
- Stress granules and virus replication/translation shutoff: https://pmc.ncbi.nlm.nih.gov/articles/PMC4574952/
- ATF4 delayed reinitiation benchmark paper: https://doi.org/10.1073/pnas.0400541101

### Long COVID and ME/CFS converge on immune persistence, viral reactivation, metabolism, and vascular biology

Long COVID and ME/CFS are heterogeneous syndromes, but several recurring themes are relevant to translation:

1. Immune activation and inflammatory signaling can persist after acute infection.
2. EBV reactivation or herpesvirus immune memory is repeatedly implicated in subsets of patients.
3. Mitochondrial, redox, and energy-metabolism abnormalities are common hypotheses.
4. Endothelial/coagulation/microvascular abnormalities are reported in Long COVID literature.
5. Exercise intolerance and post-exertional malaise imply dynamic failure of stress-response recovery, not only a static baseline defect.

These themes fit a translational-control model because immune and metabolic stress pathways converge on eIF2, mTOR, AMPK, HIF1A, mitochondrial stress signaling, and inflammatory translation programs.

Relevant sources:

- Multi-omic Long COVID study identifying early factors including viral load, autoantibodies, diabetes, and EBV viremia/reactivation as predictors of PASC: https://pubmed.ncbi.nlm.nih.gov/35216672/
- Long COVID major findings, mechanisms, and recommendations: https://www.nature.com/articles/s41579-022-00846-2
- Long COVID clinical characteristics and proposed mechanisms: https://pmc.ncbi.nlm.nih.gov/articles/PMC10171433/
- Immune activation and antigen persistence in Long COVID: https://pmc.ncbi.nlm.nih.gov/articles/PMC9996119/
- ME/CFS metabolomics and energy-metabolism abnormalities: https://pmc.ncbi.nlm.nih.gov/articles/PMC6030047/
- Long COVID coagulation/endothelial abnormalities: https://pmc.ncbi.nlm.nih.gov/articles/PMC10113134/

### DecodeME genetics is a risk-prior layer, not a reason to abandon the RDG test

The August 2025 DecodeME GWAS preprint should be part of the ME/CFS framing. It is the first large ME/CFS genetics result with genome-wide significant signals in a carefully recruited case cohort: the reported analyses used up to 15,579 ME/CFS cases and 259,909 population controls of European genetic ancestry and identified eight significant loci across the main, female-stratified, and infection-at-onset analyses. This matters because germline DNA can point toward upstream biology rather than downstream illness consequences.

It does not currently read as "ME/CFS is primarily a variant-only disease". The DecodeME authors explicitly frame ME/CFS as multifactorial rather than monogenic. The lead common-variant effects in their Table 2 are modest: reported odds ratios span roughly `0.927` to `1.100` per effect allele across the eight significant lead variants. Their liability-scale SNP heritability estimate is `9.5%` with standard error `0.6%`, and the preprint notes that this estimate is lower than some twin-based estimates while still supporting a genetic contribution. These values fit a polygenic susceptibility model where variants shift probability and biology, not a model where the current RDG atlas is fruitless unless every sample comes from the DecodeME genetic-risk population.

For this project the separation should be explicit:

| Layer | What it can answer now | What it cannot answer alone |
| --- | --- | --- |
| DecodeME germline risk | Which biological entry points are genetically prioritized for ME/CFS follow-up | Whether a public Ribo-seq infection/stress sample has ME/CFS-like translational rewiring |
| RDG atlas over public Ribo-seq | Whether infection, ISR, interferon, and fatigue-adjacent metadata branches reproducibly reroute ribosomes through uORF, ouORF, NTE/NTT, iORF, and clean-CDS nodes | Whether the same rerouting has altered penetrance in ME/CFS risk-allele carriers |
| Future patient/genotype-aware RDGs | Whether risk genotype, trigger, cell context, and translational state interact | Not available from generic public Ribo-seq alone |

The practical consequence is to continue the post-viral/fatigue RDG analysis, but keep the claim precise. Public Ribo-seq can test an acquired or trigger-responsive ribosome-programming mechanism that might sit downstream of infection and stress in many genetic backgrounds. A later ME/CFS-specific study should add patient-derived or genotype-aware material to test whether DecodeME-prioritized loci or polygenic risk modify the same RDG buffering phenotypes.

The current fatigue state score should not automatically absorb DecodeME-nearby genes as expression markers. The genes highlighted in early DecodeME interpretation, including `BTN2A2`, `OLFM4`, `RABGAP1L`, `ZNFX1`, `FBXL4`, `CA10`, `ARFGEF2`, and `CSE1L`, are better added first as a separate GWAS-prior candidate layer for transcript-architecture, translon, and FST-readiness review. They should enter the dominant-state panel only if their expression or RDG behavior improves the state definition.

Relevant source:

- DecodeME Collaboration preprint, `Initial findings from the DecodeME genome-wide association study of myalgic encephalomyelitis/chronic fatigue syndrome`: https://www.medrxiv.org/content/10.1101/2025.08.06.25333109v1

## What this suggests for RiboCrypt

The current dominant states already contain parts of the biology:

| Existing state | Relevant piece |
| --- | --- |
| Interferon / antiviral | Acute and persistent antiviral signaling |
| ISR / ER stress | eIF2-alpha/uORF-mediated stress translation |
| Hypoxia / HIF1A | Tissue stress, inflammation, mitochondrial compensation |
| OXPHOS | Energy production and mitochondrial state |
| mTOR / TOP | Global translation capacity and ribosome biogenesis |
| Senescence / SASP | Chronic inflammatory/stress phenotype |

The missing state is not simply "infection". It should be a signed post-viral fatigue state that asks whether samples combine antiviral/ISR/inflammatory activation with reduced mitochondrial/translation capacity. This could capture chronic post-viral biology better than an acute interferon-only state.

## Proposed new dominant state: Post-viral fatigue / ribosome stress

This should probably be a signed state, not a simple up-only state.

### Up genes

These genes mark antiviral persistence, ISR, inflammatory stress, and post-viral damage signals:

```text
IFI27, IFIT1, IFIT3, IFITM1, IFITM3, ISG15, MX1, OAS1, OAS2,
STAT1, STAT2, IRF7, DDX58/RIGI, IFIH1, CXCL10,
ATF4, ATF3, DDIT3, DDIT4, PPP1R15A, ASNS, TRIB3, SESN2, GDF15,
HSPA5, XBP1, HERPUD1,
NFKBIA, TNFAIP3, IL6, CXCL8, CCL2,
C3, CFB, C1QA, C1QB,
SERPINE1, F3, ICAM1, VCAM1
```

### Down genes

These genes mark mitochondrial output, mitochondrial biogenesis, and host translation capacity:

```text
NDUFS1, NDUFA9, NDUFA1, SDHB, UQCRC1, COX5A, COX6C,
ATP5F1A, ATP5MC1, CS, TFAM, PPARGC1A,
RPS6, RPLP0, RPL32, RPL13A, RPS3,
EIF4B, EIF3E, EEF2, PABPC1
```

### Initial scoring rule

Use cleaned CDS expression first:

```text
fatigue_score =
  mean(z[antiviral + ISR + inflammation + vascular/complement up genes])
  - mean(z[OXPHOS + mitochondrial biogenesis + host translation capacity down genes])
```

Then compare this state to feature-level RDG usage:

```text
RDG fatigue branch =
  samples with high fatigue_score
  versus all other samples,
  and where possible versus matched controls within infection/stress studies
```

## Candidate RDG attack points

### 1. uORF / uoORF rerouting

Use the existing RDG feature-expression framework to test whether fatigue-high or infection-high samples increase leader uORF or overlapping uoORF relative use while clean_CDS falls.

Priority marker transcripts:

```text
ATF4, PPP1R15A, DDIT4, DDIT3, HSPA5, XBP1, IFIT1, IFIT3,
ISG15, IFI27, CXCL10, OAS1, MX1, STAT1, IRF7, GDF15,
SLC2A1, LDHA, EGLN3, MKI67, SNAI2, MMP2, LMNB1
```

### 2. iORF / internal-initiation rerouting

The current pipeline intentionally ignored internal ORFs for clean-CDS masking. For fatigue/infection biology, that should become a separate feature class.

Proposed iORF classes:

| Class | Meaning | Expected phenotype |
| --- | --- | --- |
| iORF_out_of_frame | Internal ORF inside CDS, different frame | Competes with canonical protein; may create peptide fragments |
| NTT_in_frame | Internal in-frame start downstream of canonical start | N-terminal truncation; canonical protein domain lost |
| internal_near_cognate | CTG/GTG/TTG/ACG internal start support | Stress/infection-sensitive alternate initiation |
| downstream_rescue_CDS | Internal start after uoORF bypass | Similar to ATF4-like altered reinitiation |

Important caveat: same-frame internal initiation cannot be proven from total CDS coverage alone, because downstream RPFs are shared with canonical translation. We need P-site frame, start-codon peak support, 5-prime-CDS versus 3-prime-CDS imbalance, and ideally harringtonine/lactimidomycin-like TIS evidence if available.

For out-of-frame internal ORFs, the evidence can be stronger from standard Ribo-seq because frame and stop-codon behavior can separate the iORF from the canonical CDS. For in-frame internal truncations, the evidence must come from coverage geometry and start-site enrichment:

```text
strong iORF/NTT evidence =
  low canonical-start-region coverage
  preserved downstream CDS coverage
  local start peak at internal AUG / near-cognate codon
  compatible P-site frame
  reproducible in matched infected/stress branch
```

### 3. Feature ratios to add

For every gene/transcript in fatigue/infection branches:

```text
leader_uORF_relative_use
overlapping_uORF_relative_use
internal_ORF_relative_use
clean_CDS_relative_use
CDS_5prime_to_3prime_coverage_ratio
canonical_start_region_to_downstream_CDS_ratio
iORF_start_peak_score
frame_purity_score
```

The key iORF phenotype is not only "iORF high". It is:

```text
upstream canonical CDS low
downstream CDS high or preserved
internal start/frame signal high
normal protein N-terminus likely lost
```

### 4. Infection/fatigue metadata branches

The current metadata should be mined for:

```text
infected, infect, viral, virus, SARS, COVID, EBV, IFN, interferon,
cytokines, LPS, TNF, stress, ER stress, arsenite, thapsigargin,
hypoxia, starvation, heat shock, NaCl, inhibitor, drug, treatment
```

Then each infection/stress branch should be split by design:

1. Direct infected versus mock/control within study and cell line.
2. Cytokine/interferon versus control.
3. Stress/inhibitor versus vehicle/control.
4. Chronic-like branches, if metadata has repeated or persistent infection labels.

## Candidate fatigue-state gene modules

### Conservative v1 module

This avoids making the state too broad:

```text
Up:
IFI27, IFIT1, IFIT3, IFITM1, ISG15, MX1, OAS1, STAT1, IRF7,
ATF4, ATF3, DDIT3, DDIT4, PPP1R15A, ASNS, TRIB3, GDF15,
HSPA5, XBP1, HERPUD1,
CXCL10, CCL2, IL6, NFKBIA, TNFAIP3, SERPINE1

Down:
NDUFS1, NDUFA9, SDHB, UQCRC1, COX5A, ATP5F1A, ATP5MC1,
CS, TFAM, PPARGC1A,
RPS6, RPLP0, RPL32, RPS3, EIF4B, EEF2, PABPC1
```

### More exploratory v2 module

Add complement/coagulation/endothelial and EBV/B-cell context:

```text
Extra up:
C3, CFB, C1QA, C1QB, C1QC, VWF, F3, ICAM1, VCAM1,
CD74, HLA-DRA, HLA-DRB1, MS4A1, CD79A
```

Use this only if the samples include relevant tissue/blood/immune contexts. In cell-line Ribo-seq, these may become tissue-composition markers rather than fatigue biology.

The same gene set is saved in machine-readable form as:

```text
dominant_cell_states/results/curated_inputs/postviral_fatigue_dominant_state_genes.csv
```

Use the conservative rows first for scoring. The exploratory rows should be reported separately or used only as sensitivity analysis.

## Practical model tests

The first model should not try to prove chronic fatigue directly. It should test a narrower translation hypothesis:

```text
Do infection-like and fatigue-like sample states show reproducible redistribution of
ribosome allocation from clean canonical CDS toward uORF, overlapping uORF, NTT/NTE,
or internal ORF features?
```

Recommended outputs:

| Output | Purpose |
| --- | --- |
| fatigue_state_scores.csv | Sample-level signed score and component scores |
| fatigue_metadata_enrichment.csv | Which metadata terms are enriched for the state |
| fatigue_rdg_feature_usage.csv | Feature-level usage for high-score branches |
| fatigue_iORF_candidates.csv | Genes with internal ORF / NTT-like evidence |
| fatigue_matched_control_effects.csv | Direct contrast against matched WT/mock/vehicle samples |

Minimum columns for candidate ranking:

```text
Gene, transcript_id, feature_class, state, branch_term, n_in, n_out,
clean_CDS_delta, uORF_delta, ouORF_delta, iORF_delta,
canonical_start_to_downstream_ratio_delta,
frame_score_delta, start_peak_score_delta,
matched_control_delta, replicated_designs, q_value
```

The most useful biological output is not a list of high scores. It is an interpretable effect size:

```text
Predicted canonical protein-output buffering/loss for each gene and condition.
```

Example target statement:

```text
In interferon/post-viral-like samples, gene X shows 35% lower clean-CDS allocation,
mostly explained by increased overlapping uORF and internal ORF usage.
```

## What would make the hypothesis stronger

Evidence tiers:

| Tier | Evidence |
| --- | --- |
| 1 | Fatigue/post-viral score enriches infected, cytokine, interferon, or viral-reactivation metadata branches |
| 2 | Fatigue-high branches show lower clean_CDS relative use in multiple RDG genes |
| 3 | The drop is explained by higher uORF/uoORF/iORF feature use, not low read count |
| 4 | Effects replicate across multiple studies/designs, not only one infection experiment |
| 5 | iORF/NTT effects show frame/TIS support and 5-prime-CDS depletion |
| 6 | Candidate genes connect to phenotype: mitochondrial output, immune exhaustion, endothelial/coagulation, neuronal/metabolic function |

## Immediate implementation plan

1. Add `Post-viral fatigue / ribosome stress` as a signed dominant state.
2. Save both conservative and exploratory gene-set definitions, but score conservative first.
3. Rerun clean-CDS dominant-state scoring and metadata enrichment.
4. Add infection/fatigue-specific branch sets to RDG generation.
5. Extend RDG feature annotation to include internal ORFs as a separate feature class, not as CDS masks.
6. Add iORF/NTT diagnostics:
   - iORF feature coverage
   - canonical start-to-downstream CDS imbalance
   - frame support
   - start-codon support
7. Rank genes by:
   - fatigue-state enrichment
   - clean_CDS loss
   - uORF/uoORF/iORF gain
   - matched-control support
   - design replication
   - biological relevance to post-viral fatigue

## Main risk

The major risk is confounding. Infection, interferon, stress, cell line, tissue, fraction, author, and study are heavily nested in public Ribo-seq metadata. A "fatigue" state can easily become an "infection experiment" state or a "blood/immune/tissue composition" state.

The solution is not to avoid the state, but to treat it as a hypothesis generator and require matched-control RDG evidence plus design replication before interpreting gene-specific results.

## Bottom line

The idea is plausible and worth testing. Viral infection is already known to reprogram host translation through exactly the mechanisms that should alter uORF, uoORF, and iORF usage. Long COVID and ME/CFS literature repeatedly points toward immune persistence, EBV/herpesvirus reactivation, mitochondrial stress, vascular/coagulation changes, and post-exertional failure. These are all compatible with a chronic ribosome-decision-state model.

The next concrete RiboCrypt goal should be:

```text
Build a post-viral fatigue / ribosome stress dominant state,
then test whether fatigue-high and infection-high branches show replicated
loss of clean canonical CDS allocation through increased uORF/uoORF/iORF routing.
```
