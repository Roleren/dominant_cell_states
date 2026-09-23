#!/usr/bin/env bash

set -euo pipefail

# Download only the ribosome-protected-fragment libraries required for the
# tRNA-Glu(UUC) and matched tRNA-Arg(CCG) EIF4G2 audit. Total-transcript
# libraries are deliberately excluded.

output_dir="${GOODARZI_TRNA_FASTQ_DIR:-/media/roler/S/data/Bio_data/runtime_inputs/dominant_cell_states/Goodarzi_trna_fastq}"
selected_runs="${GOODARZI_TRNA_RUNS:-all}"
download_tmp_dir="${GOODARZI_TRNA_DOWNLOAD_TMPDIR:-$output_dir}"
mkdir -p "$output_dir" "$download_tmp_dir"

manifest=(
  "SRR3131003|https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR313/003/SRR3131003/SRR3131003.fastq.gz|8c1c8b95e5241f560319c9d1a1ccada0"
  "SRR3131004|https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR313/004/SRR3131004/SRR3131004.fastq.gz|c2c6ec8016e99b8f806836a8bf87be0d"
  "SRR3131007|https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR313/007/SRR3131007/SRR3131007.fastq.gz|225fa19ad5acd5947ac7fbddfd2da558"
  "SRR3131008|https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR313/008/SRR3131008/SRR3131008.fastq.gz|e4e17f2d32af91ca010d486f8a3e4f6d"
  "SRR3131011|https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR313/001/SRR3131011/SRR3131011.fastq.gz|aefab321e1b2d8af6a52372a04342839"
  "SRR3131012|https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR313/002/SRR3131012/SRR3131012.fastq.gz|349353e05f037aa6e81d54f05f78c5d8"
  "SRR3131015|https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR313/005/SRR3131015/SRR3131015.fastq.gz|d9d283f1ebae4951c28600089e275d22"
  "SRR3131016|https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR313/006/SRR3131016/SRR3131016.fastq.gz|ee03ba1eeebbc9d245402597055f43fb"
  "SRR3131019|https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR313/009/SRR3131019/SRR3131019.fastq.gz|72e720b1121bd868696ffe88341d9570"
  "SRR3131021|https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR313/001/SRR3131021/SRR3131021.fastq.gz|d7c263f6f568a8784987dd0a1c2228a1"
  "SRR3131023|https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR313/003/SRR3131023/SRR3131023.fastq.gz|59e2cce4319f8c4dbd0dc5f5ec3c6dbc"
  "SRR3131025|https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR313/005/SRR3131025/SRR3131025.fastq.gz|dcb7cd39211847c9d3232f0c15808c68"
  "SRR3130940|https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR313/000/SRR3130940/SRR3130940.fastq.gz|c6460062262445319b0bd7899a6a547f"
  "SRR3130941|https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR313/001/SRR3130941/SRR3130941.fastq.gz|d7b4c9961e49fea0607d75e24583c612"
  "SRR3130944|https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR313/004/SRR3130944/SRR3130944.fastq.gz|2d63ad5f45cbf8a75d04112e73bf276d"
  "SRR3130945|https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR313/005/SRR3130945/SRR3130945.fastq.gz|90a6df1aeb80473d71ffdcd9229807c4"
)

run_selected() {
  local run="$1"
  if [[ "$selected_runs" == "all" ]]; then
    return 0
  fi
  [[ " $selected_runs " == *" $run "* ]]
}

for entry in "${manifest[@]}"; do
  IFS='|' read -r run url expected_md5 <<<"$entry"
  if ! run_selected "$run"; then
    continue
  fi
  destination="$output_dir/$run.fastq.gz"
  if [[ -s "$destination" ]]; then
    if printf '%s  %s\n' "$expected_md5" "$destination" | \
         md5sum -c --status && gzip -t "$destination"; then
      echo "Using verified FASTQ: $destination"
      continue
    fi
    echo "Existing FASTQ failed integrity checks; move it aside before retrying: $destination" >&2
    exit 1
  fi

  temporary="$download_tmp_dir/$run.fastq.gz.incomplete"
  echo "Downloading or resuming $run"
  wget -c --progress=dot:giga -O "$temporary" "$url"
  printf '%s  %s\n' "$expected_md5" "$temporary" | md5sum -c -
  gzip -t "$temporary"
  if [[ "$download_tmp_dir" == "$output_dir" ]]; then
    mv "$temporary" "$destination"
  else
    destination_temporary="$destination.incomplete"
    cp "$temporary" "$destination_temporary"
    printf '%s  %s\n' "$expected_md5" "$destination_temporary" | md5sum -c -
    gzip -t "$destination_temporary"
    mv "$destination_temporary" "$destination"
    rm -f "$temporary"
  fi
done

echo "Goodarzi tRNA RPF FASTQ download complete: $output_dir"
