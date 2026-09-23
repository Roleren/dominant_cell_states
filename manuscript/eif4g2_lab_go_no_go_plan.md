# EIF4G2 rapid lab go/no-go plan

Updated: 2026-08-19

## Decision to make

Before investing in new Ribo-seq, endogenous editing or immunopeptidomics, ask
one narrow question:

> Does tRNA-Glu selectively change translation through the native EIF4G2
> leader in a manner that requires the uORF1 start codon?

UV is an orthogonal branch-control, not the central mechanism. It asks whether
the same EIF4G2 uORF1 can be switched by a defined ribosome-stress condition.
A positive UV result strengthens the locus and assay. It cannot rescue a
negative tRNA-Glu result.

| Decision stage | Smallest experiment | What a pass means | What failure means |
|---|---|---|---|
| Screen A: main story | MDA-MB-231, vector/Glu/Arg × WT/start mutant | tRNA-Glu causally engages uORF1-dependent leader output | Stop or redesign the tRNA mechanism before expensive work |
| Screen B: orthogonal check | HeLa, untreated/UV × WT/start mutant | UV independently engages the same EIF4G2 branch and validates the reporter | UV remains computational support only; it does not kill a positive tRNA result |
| Combined second stage | MDA-MB-231, tRNA × UV × construct | Reveals additivity, saturation or occlusion after both axes work alone | Do not attempt unless Screen A passes |
| Endogenous/immune stage | uORF1 editing, endogenous translation/protein, then HLA peptide | Converts the reporter result into a biological and immune mechanism | Defer until the cheap causal screen is positive |

## Build the smallest informative reporter pair

Use the complete 308-nt leader of EIF4G2 ENST00000339995 and retain the native
GUG main start and immediate EIF4G2 coding context in-frame with a sensitive
luciferase reporter. Make two otherwise identical constructs:

1. **WT leader:** native ATG uORF1, TTG uORF2 and GUG main start.
2. **uORF1-start mutant:** ATG changed to a non-start codon, preferably the
   ATG-to-ATA change already used for the published EIF4G2 uORF reporter.

Sequence the entire leader and junctions. Include a constitutive second
luciferase or equivalent transfection control. Quantify reporter mRNA from the
same wells or matched wells; the primary quantity is translation output
normalized to both transfection and reporter RNA.

This downstream reporter is the cheapest first test of branch function. If it
works, add a separate uORF1-initiation reporter in round 2 rather than heavily
modifying the 51-nt native uORF in the initial construct. The published study
already cloned EIF4G2-leader WT and uORF-start-mutant reporters and used
ATG-to-ATA start mutation, so obtaining those plasmids or sequences could
shorten setup. Its immune reporter replaced the uORF peptide with SIINFEKL and
therefore should be treated as a useful positive-control design, not a native
mechanistic reporter ([Kowar *et al.*, 2025](https://www.nature.com/articles/s41467-025-63405-2)).

## Screen A — does the tRNA story survive a causal reporter?

Use MDA-MB-231, matching the source tRNA experiment. Cross:

- vector, tRNA-Glu(UUC) and tRNA-Arg(CCG);
- WT and uORF1-start-mutant EIF4G2 reporters;
- at least three independently cultured biological experiments, with technical
  wells nested within each experiment rather than counted as replication.

Verify that the tRNA-Glu construct actually increases the intended mature tRNA
species. A Northern blot is preferable for the first decisive experiment
because tRNA modifications can make ordinary reverse transcription
misleading. Also record precursor and prominent fragment species if visible.
If mature tRNA-Glu is not demonstrably increased, the screen is technically
uninformative rather than negative.

Primary analysis: fit the reporter output with tRNA perturbation, leader
construct and their interaction. The key result is the tRNA-Glu × intact-uORF1
interaction. Within-construct significance tests are not a substitute.
tRNA-Arg is the matched specificity control, not merely another negative well.

### Go

Advance when all of the following are true:

- mature tRNA-Glu perturbation is verified in every independent experiment;
- the WT reporter changes under tRNA-Glu in the direction predicted by the
  public occupancy result, and the change is attenuated by uORF1-start
  mutation;
- the interaction has the same direction in all three biological experiments,
  its interval excludes zero, and the effect is large enough to be practical
  (pre-register approximately 20% as a useful first-screen benchmark unless
  pilot precision justifies a different threshold);
- tRNA-Arg does not reproduce the same WT-versus-mutant interaction;
- reporter RNA abundance does not account for the translation change.

### Stop or redesign

- **Stop the tRNA mechanism paper:** verified tRNA-Glu and tRNA-Arg produce the
  same construct interaction, or verified tRNA-Glu has no reproducible
  WT-versus-start-mutant interaction.
- **Redesign the reporter:** WT and start mutant differ at baseline but neither
  responds selectively, the GUG fusion has negligible signal, or reporter RNA
  changes explain the output.
- **Repeat the perturbation, not the conclusion:** mature tRNA-Glu was not
  increased or the construct produced mostly precursor/fragment species.

## Screen B — can UV validate the same branch independently?

Run this first in HeLa, matching the raw UV dataset
([GSE141459](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE141459)).
Cross untreated and UV-treated cells with the same WT and uORF1-start-mutant
reporters in at least three independently cultured experiments.

Do a small UV dose/time pilot rather than copying a dose blindly across lamps
and instruments. Choose the lowest regimen that reproducibly activates a
ribotoxic-stress marker while retaining sufficient viability and reporter RNA.
Record lamp wavelength, calibrated dose, recovery time and viability. Include
phosphorylated p38/JNK and, if already routine in the lab, phospho-eIF2alpha as
engagement controls. The source study links UV-induced ribosome collisions to
ZAK-dependent stress signalling, but the first EIF4G2 test need not establish
that whole pathway ([Wu *et al.*, 2020](https://pubmed.ncbi.nlm.nih.gov/32610081/)).

### UV gate

Pass when UV changes normalized WT reporter output and the response is
materially attenuated by uORF1-start mutation, with stable reporter RNA and
the same interaction direction across biological experiments. Failure of the
stress markers makes the experiment technically uninformative. A stress-
positive but start-mutation-insensitive response argues that UV changes the
reporter through another route and does not validate this uORF1 branch.

If this gate passes, a second-stage ZAK perturbation can ask whether collision
signalling is required. That is mechanistic refinement, not part of the cheap
initial decision.

## Combine tRNA and UV only after the two screens

If Screen A and Screen B both pass, cross vector/tRNA-Glu/tRNA-Arg with
untreated/UV and WT/start-mutant reporters in MDA-MB-231. This experiment asks
whether tRNA-Glu and UV are additive, occlusive or independent. Do not require
additivity: saturation or occlusion could be informative. The decisive model
contains tRNA, UV, construct and their interactions, with biological experiment
as a blocking factor.

Interpretation:

- **Both single-axis screens pass:** strong reason to invest in endogenous
  validation.
- **tRNA passes, UV fails:** retain the main tRNA story; UV remains a public-
  data recurrence rather than an experimental pillar.
- **UV passes, tRNA fails:** EIF4G2 is a real stress-responsive uORF branch, but
  abandon or radically narrow the tRNA paper.
- **Both fail with adequate controls:** do not spend on Ribo-seq or HLA work.

## Spending ladder after a go

1. **Cheap:** repeat/lock the reporter interaction; verify mature tRNA,
   reporter RNA, stress engagement and endogenous DAP5 by immunoblot.
2. **Moderate:** add an uORF1-initiation reporter, GAG-to-GAA and codon-5/pause
   mutants, and mature-tRNA-versus-fragment perturbations.
3. **Expensive:** edit the endogenous uORF1 start, perform focused endogenous
   translation profiling, and quantify DAP5 protein with orthogonal methods.
4. **Most expensive:** targeted HLA immunopeptidomics for the exact EIF4G2-uORF
   peptide, only after the endogenous branch is causal.

The paper should call UV a confirmation/orthogonal recurrence. The title and
central model live or die on Screen A; the immune claim lives or dies later on
the endogenous peptide measurement.
