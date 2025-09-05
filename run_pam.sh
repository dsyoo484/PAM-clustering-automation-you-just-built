#!/bin/bash
# run_pam.sh : Bash wrapper to call pam_clustering.R

INPUT="/BiO/Access/home/user9/Data/dataset1/ratio/ASV_table.blast_NCBI_16S_L7.txt"
OUTDIR="/BiO/Access/home/user9/PAM_results"
CLUSTER_N=10
COUNT=TRUE
STAGE=6   # ✅ 전체 실행

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
RSCRIPT="${SCRIPT_DIR}/pam_clustering.R"

command -v Rscript >/dev/null 2>&1 || { echo "❌ Rscript not found"; exit 1; }
[[ -f "$RSCRIPT" ]] || { echo "❌ R script not found: $RSCRIPT"; exit 1; }

# 입력을 절대경로로
if [[ -f "$INPUT" ]]; then
  INPUT_ABS="$(readlink -f "$INPUT")"
else
  echo "❌ Input not found: $INPUT"
  echo "현재 위치: $(pwd)"
  echo "여기 파일 목록:"; ls -l
  exit 1
fi

mkdir -p "$OUTDIR/plot"
OUTDIR_ABS="$(readlink -f "$OUTDIR")"

echo "▶ INPUT  : $INPUT_ABS"
echo "▶ OUTDIR : $OUTDIR_ABS"
echo "▶ k max  : $CLUSTER_N"
echo "▶ TSS    : $COUNT"
echo "▶ STAGE  : $STAGE (1=read ... 6=all)"

Rscript "$RSCRIPT" \
  --input "$INPUT_ABS" \
  --outdir "$OUTDIR_ABS" \
  --cluster_n "$CLUSTER_N" \
  --count "$COUNT" \
  --stage "$STAGE"

echo "✅ Done → $OUTDIR_ABS"

