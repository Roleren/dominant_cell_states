#!/usr/bin/env bash

set -euo pipefail

# Trim and uniquely align the tRNA-Glu(UUC), tRNA-Arg(CCG), and matched-control
# RPF libraries while retaining biological insert length. Runtime products live
# outside the source repository.

input_dir="${GOODARZI_TRNA_FASTQ_DIR:-/media/roler/S/data/Bio_data/runtime_inputs/dominant_cell_states/Goodarzi_trna_fastq}"
output_dir="${GOODARZI_TRNA_OUTPUT_DIR:-/media/roler/S/data/Bio_data/runtime_inputs/dominant_cell_states/Goodarzi_trna_length_resolved}"
skewer_bin="${SKEWER_BIN:-/media/roler/S/data/Bio_data/runtime_inputs/dominant_cell_states/tools/skewer/skewer}"
star_bin="${STAR_BIN:-/home/roler/bin/STAR-2.7.4a/bin/Linux_x86_64/STAR}"
star_index="${STAR_INDEX:-/media/roler/S/data/Bio_data/references/Homo_sapiens_GRCh38_101/STAR_index/genomeDir}"
threads="${GOODARZI_TRNA_THREADS:-8}"
force="${GOODARZI_TRNA_FORCE:-false}"

default_runs="SRR3131003 SRR3131004 SRR3131007 SRR3131008 SRR3131011 SRR3131012 SRR3131015 SRR3131016 SRR3131019 SRR3131021 SRR3131023 SRR3131025 SRR3130940 SRR3130941 SRR3130944 SRR3130945"
read -r -a runs <<<"${GOODARZI_TRNA_RUNS:-$default_runs}"
adapter="AGATCGGAAGAGCACACGTCT"

if [[ ! -x "$skewer_bin" ]]; then
  echo "Skewer executable not found: $skewer_bin" >&2
  exit 1
fi
if [[ ! -x "$star_bin" ]]; then
  echo "STAR executable not found: $star_bin" >&2
  exit 1
fi
if [[ ! -d "$star_index" ]]; then
  echo "STAR index not found: $star_index" >&2
  exit 1
fi

mkdir -p "$output_dir/trimmed" "$output_dir/aligned" "$output_dir/logs"

for run in "${runs[@]}"; do
  if [[ ! "$run" =~ ^SRR313(09(40|41|44|45)|10(03|04|07|08|11|12|15|16|19|21|23|25))$ ]]; then
    echo "Unsupported Goodarzi tRNA RPF run: $run" >&2
    exit 1
  fi

  input_fastq="$input_dir/$run.fastq.gz"
  trimmed_fastq="$output_dir/trimmed/$run.trimmed_15_40nt.fastq.gz"
  trim_log="$output_dir/logs/$run.skewer.log"

  if [[ ! -s "$input_fastq" ]]; then
    echo "Missing FASTQ: $input_fastq" >&2
    exit 1
  fi
  gzip -t "$input_fastq"

  if [[ "$force" == "true" || ! -s "$trimmed_fastq" ]]; then
    temporary_fastq="$trimmed_fastq.incomplete"
    echo "Trimming $run"
    "$skewer_bin" \
      -m tail \
      -x "$adapter" \
      -l 15 \
      -L 40 \
      -r 0.1 \
      -t "$threads" \
      -1 \
      "$input_fastq" \
      2>"$trim_log" | gzip -c >"$temporary_fastq"
    gzip -t "$temporary_fastq"
    mv "$temporary_fastq" "$trimmed_fastq"
  else
    echo "Using existing trimmed FASTQ: $trimmed_fastq"
  fi
  gzip -t "$trimmed_fastq"

  alignment_dir="$output_dir/aligned/$run"
  aligned_bam="$alignment_dir/Aligned.sortedByCoord.out.bam"
  mkdir -p "$alignment_dir"
  alignment_ready="false"
  if [[ -s "$aligned_bam" ]] && samtools quickcheck "$aligned_bam"; then
    if [[ ! -s "$aligned_bam.bai" ]]; then
      samtools index -@ "$threads" "$aligned_bam"
    fi
    if samtools idxstats "$aligned_bam" >/dev/null; then
      alignment_ready="true"
    fi
  fi

  if [[ "$force" == "true" || "$alignment_ready" != "true" ]]; then
    recovery_stamp="$(date -u +%Y%m%dT%H%M%SZ)"
    recovery_reason="incomplete"
    if [[ "$force" == "true" && "$alignment_ready" == "true" ]]; then
      recovery_reason="forced_previous"
    fi
    if [[ -e "$aligned_bam" ]]; then
      mv "$aligned_bam" "$aligned_bam.${recovery_reason}_$recovery_stamp"
    fi
    if [[ -e "$aligned_bam.bai" ]]; then
      mv "$aligned_bam.bai" \
        "$aligned_bam.bai.${recovery_reason}_$recovery_stamp"
    fi
    if [[ -e "$alignment_dir/_STARtmp" ]]; then
      mv "$alignment_dir/_STARtmp" \
        "$alignment_dir/_STARtmp.${recovery_reason}_$recovery_stamp"
    fi
    echo "Aligning $run"
    "$star_bin" \
      --genomeDir "$star_index" \
      --readFilesIn "$trimmed_fastq" \
      --readFilesCommand zcat \
      --runThreadN "$threads" \
      --outFileNamePrefix "$alignment_dir/" \
      --outSAMtype BAM SortedByCoordinate \
      --outSAMattributes NH HI AS nM MD \
      --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
      --outFilterMultimapNmax 1 \
      --outFilterMismatchNoverLmax 0.1 \
      --limitBAMsortRAM 12000000000
    samtools index -@ "$threads" "$aligned_bam"
  else
    echo "Using existing alignment: $aligned_bam"
  fi
  samtools quickcheck "$aligned_bam"
done

echo "Goodarzi tRNA length-preserving preprocessing complete: $output_dir"
