#!/usr/bin/env bash

set -euo pipefail

# Reprocess the four GSE141459 monosome Ribo-seq runs while retaining read
# length. Runtime data and BAMs remain outside the repository.

input_dir="${GSE141459_FASTQ_DIR:-/media/roler/S/data/Bio_data/runtime_inputs/dominant_cell_states/GSE141459_fastq}"
output_dir="${GSE141459_OUTPUT_DIR:-/media/roler/S/data/Bio_data/runtime_inputs/dominant_cell_states/GSE141459_length_resolved}"
skewer_bin="${SKEWER_BIN:-/media/roler/S/data/Bio_data/runtime_inputs/dominant_cell_states/tools/skewer/skewer}"
star_bin="${STAR_BIN:-/home/roler/bin/STAR-2.7.4a/bin/Linux_x86_64/STAR}"
star_index="${STAR_INDEX:-/media/roler/S/data/Bio_data/references/Homo_sapiens_GRCh38_101/STAR_index/genomeDir}"
threads="${GSE141459_THREADS:-8}"
force="${GSE141459_FORCE:-false}"

default_runs="SRR11569080 SRR11569081 SRR11569082 SRR11569083"
read -r -a runs <<<"${GSE141459_RUNS:-$default_runs}"
adapter="NNNNNNCACTCGGGCACCAAGGA"

if [[ "${#runs[@]}" -eq 0 ]]; then
  echo "GSE141459_RUNS resolved to an empty run list." >&2
  exit 1
fi
for run in "${runs[@]}"; do
  if [[ ! "$run" =~ ^SRR1156908[0-3]$ ]]; then
    echo "Unsupported GSE141459 run: $run" >&2
    exit 1
  fi
done

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
  input_fastq="$input_dir/$run.fastq.gz"
  trimmed_fastq="$output_dir/trimmed/$run.trimmed_15_34nt.fastq.gz"
  trim_log="$output_dir/logs/$run.skewer.log"

  if [[ ! -s "$input_fastq" ]]; then
    echo "Missing FASTQ: $input_fastq" >&2
    exit 1
  fi
  if ! gzip -t "$input_fastq"; then
    echo "FASTQ is incomplete or corrupt: $input_fastq" >&2
    exit 1
  fi

  if [[ "$force" == "true" || ! -s "$trimmed_fastq" ]]; then
    temporary_fastq="$output_dir/trimmed/$run.trimmed_15_34nt.fastq.gz.incomplete"
    echo "Trimming $run"
    "$skewer_bin" \
      -m tail \
      -x "$adapter" \
      -l 19 \
      -L 38 \
      -r 0.1 \
      -t "$threads" \
      -1 \
      "$input_fastq" \
      2>"$trim_log" |
      awk '
        NR % 4 == 1 { print; next }
        NR % 4 == 2 { print substr($0, 5); next }
        NR % 4 == 3 { print; next }
        NR % 4 == 0 { print substr($0, 5) }
      ' |
      gzip -c >"$temporary_fastq"
    gzip -t "$temporary_fastq"
    mv "$temporary_fastq" "$trimmed_fastq"
  else
    echo "Using existing trimmed FASTQ: $trimmed_fastq"
  fi

  alignment_dir="$output_dir/aligned/$run"
  aligned_bam="$alignment_dir/Aligned.sortedByCoord.out.bam"
  mkdir -p "$alignment_dir"
  if [[ "$force" == "true" || ! -s "$aligned_bam" ]]; then
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
    if [[ ! -s "$aligned_bam.bai" ]]; then
      samtools index -@ "$threads" "$aligned_bam"
    fi
  fi
done

echo "GSE141459 length-preserving preprocessing complete: $output_dir"
