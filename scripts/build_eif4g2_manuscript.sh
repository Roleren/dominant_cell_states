#!/usr/bin/env bash

set -euo pipefail

analysis_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
manuscript_dir="$analysis_dir/manuscript"
output_dir="$manuscript_dir/generated"
pandoc_bin="${PANDOC_BIN:-/home/roler/bin/pandoc}"

if [[ ! -x "$pandoc_bin" ]]; then
  pandoc_bin="$(command -v pandoc || true)"
fi
if [[ -z "$pandoc_bin" || ! -x "$pandoc_bin" ]]; then
  echo "pandoc was not found. Set PANDOC_BIN to an executable." >&2
  exit 1
fi

mkdir -p "$output_dir"

env -u LC_ALL Rscript --vanilla \
  "$analysis_dir/scripts/validate_eif4g2_manuscript.R"

"$pandoc_bin" \
  "$manuscript_dir/eif4g2_article_draft.md" \
  "$manuscript_dir/eif4g2_figure_plan_and_legends.md" \
  "$manuscript_dir/eif4g2_methods.md" \
  "$manuscript_dir/eif4g2_extended_data.md" \
  "$manuscript_dir/eif4g2_lab_go_no_go_plan.md" \
  "$manuscript_dir/eif4g2_claims_evidence_and_gaps.md" \
  --from=gfm \
  --to=docx \
  --standalone \
  --toc \
  --metadata title="EIF4G2 working manuscript and evidence package" \
  --output="$output_dir/eif4g2_working_manuscript.docx"

echo "Wrote $output_dir/eif4g2_working_manuscript.docx"
