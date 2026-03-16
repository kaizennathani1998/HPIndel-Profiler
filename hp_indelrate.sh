#!/usr/bin/env bash
set -euo pipefail

# Usage:
#   ./hp_indelrate.sh /path/to/sample.bam /path/to/summary.txt \
#     [reference.fa] [reference.fa.fai] [ezh2_hp.bed] [ezh2_nonhp.bed]
#
# Required arguments:
#   1) BAM file
#   2) Summary TXT file
#
# Optional arguments:
#   3) Reference FASTA path
#   4) Reference FASTA index (.fai) path
#   5) HP BED path
#   6) non-HP BED path
#
# Defaults for optional arguments:
#   <script_dir>/resources/hg19.fa
#   <script_dir>/resources/hg19.fa.fai
#   <script_dir>/resources/ezh2_hp.bed
#   <script_dir>/resources/ezh2_nonhp.bed
#
# Output:
#   # Not included for workflow runs -  Writes a per-BAM txt file in the BAM directory named:  <bam_prefix>.ezh2HP_indelRate.txt
#   1) Appends one tab-delimited row to the summary file:
#        BAM  Run  Sample  HP_BASES  NONHP_BASES  HP_INDELS  NONHP_INDELS
#        HP_indel_rate  nonHP_indel_rate  HP_over_nonHP_ratio

if [[ $# -lt 2 || $# -gt 6 ]]; then
  echo "Usage: $0 /path/to/sample.bam /path/to/summary.txt [reference.fa] [reference.fa.fai] [ezh2_hp.bed] [ezh2_nonhp.bed]" >&2
  exit 1
fi

BAM="$1"
SUMMARY_TXT="$2"

# --- resolve default resource paths relative to script location ---
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RESOURCE_DIR="${SCRIPT_DIR}/resources"

DEFAULT_REF="${RESOURCE_DIR}/hg19.fa"
DEFAULT_REF_FAI="${RESOURCE_DIR}/hg19.fa.fai"
DEFAULT_HP_BED="${RESOURCE_DIR}/ezh2_hp.bed"
DEFAULT_NONHP_BED="${RESOURCE_DIR}/ezh2_nonhp.bed"

# --- optional overrides ---
REF="${3:-$DEFAULT_REF}"
REF_FAI="${4:-$DEFAULT_REF_FAI}"
HP_BED="${5:-$DEFAULT_HP_BED}"
NONHP_BED="${6:-$DEFAULT_NONHP_BED}"

# --- derive output path (same dir as BAM) ---
BAM_DIR="$(cd "$(dirname "$BAM")" && pwd)"
BAM_BASE="$(basename "$BAM")"
BAM_PREFIX="${BAM_BASE%.bam}"
OUT_TXT="${BAM_DIR}/${BAM_PREFIX}.ezh2HP_indelRate.txt"

# --- derive run name ---
# Change this if your BAM directory structure is different.
Run="$(basename "$(dirname "$BAM")")"
if [[ -z "${Run:-}" ]]; then
  Run="NA"
fi

# --- sample name is BAM prefix ---
Sample="$BAM_PREFIX"

# --- sanity checks ---
if [[ ! -f "$BAM" ]]; then
  echo "ERROR: BAM not found: $BAM" >&2
  exit 1
fi

# Ensure BAM index exists (.bam.bai or .bai)
if [[ -f "${BAM}.bai" ]]; then
  BAM_INDEX="${BAM}.bai"
elif [[ -f "${BAM%.bam}.bai" ]]; then
  BAM_INDEX="${BAM%.bam}.bai"
else
  echo "ERROR: BAM index not found for: $BAM" >&2
  echo "       Expected either: ${BAM}.bai or ${BAM%.bam}.bai" >&2
  echo "       Create it with: samtools index \"$BAM\"" >&2
  exit 1
fi

if [[ ! -f "$REF" ]]; then
  echo "ERROR: Reference FASTA not found: $REF" >&2
  exit 1
fi

if [[ ! -f "$REF_FAI" ]]; then
  echo "ERROR: Reference FASTA index not found: $REF_FAI" >&2
  exit 1
fi

if [[ ! -f "$HP_BED" ]]; then
  echo "ERROR: HP BED not found: $HP_BED" >&2
  exit 1
fi

if [[ ! -f "$NONHP_BED" ]]; then
  echo "ERROR: non-HP BED not found: $NONHP_BED" >&2
  exit 1
fi

# --- consistency check: ensure .fai matches .fa when both defaults are not used ---
# This only checks basename convention, not file content.
EXPECTED_FAI="${REF}.fai"
if [[ "$REF_FAI" != "$EXPECTED_FAI" ]]; then
  echo "WARNING: Provided FASTA index path does not match the usual convention (${EXPECTED_FAI})." >&2
  echo "         Proceeding anyway because an explicit .fai path was provided." >&2
fi

# --- denominators: sum of depth across each region ---
HP_BASES=$(samtools depth -a -b "$HP_BED" "$BAM" | awk '{s+=$3} END{print s+0}')
NONHP_BASES=$(samtools depth -a -b "$NONHP_BED" "$BAM" | awk '{s+=$3} END{print s+0}')

# --- numerator: indel events from mpileup read-bases string ---
count_indels_mpileup () {
  local bed="$1"
  samtools mpileup -f "$REF" -l "$bed" -aa -A -B -Q 0 -q 0 "$BAM" 2>/dev/null \
  | awk '
    function count_indels(s,    n,c,i,len) {
      n=0
      for (i=1; i<=length(s); i++) {
        c=substr(s,i,1)
        if (c=="+" || c=="-") {
          i++
          len=""
          while (i<=length(s) && substr(s,i,1) ~ /[0-9]/) {
            len = len substr(s,i,1)
            i++
          }
          if (len=="") continue
          i += len - 1
          n++
        }
      }
      return n
    }
    { total += count_indels($5) }
    END { print total+0 }'
}

HP_INDELS=$(count_indels_mpileup "$HP_BED")
NONHP_INDELS=$(count_indels_mpileup "$NONHP_BED")

# --- compute rates + ratio ---
read -r HP_RATE NONHP_RATE RATIO < <(
awk -v hpI="$HP_INDELS" -v hpB="$HP_BASES" \
    -v nhI="$NONHP_INDELS" -v nhB="$NONHP_BASES" '
BEGIN {
  hp_rate = (hpB>0)? hpI/hpB : 0
  nh_rate = (nhB>0)? nhI/nhB : 0
  ratio   = (nh_rate>0)? hp_rate/nh_rate : 0
  printf("%.10f %.10f %.6f\n", hp_rate, nh_rate, ratio)
}'
)

# --- write detailed per-BAM output ---

## --- Skipped for workflow run ----


# --- ensure summary directory exists ---
SUMMARY_DIR="$(dirname "$SUMMARY_TXT")"
mkdir -p "$SUMMARY_DIR"

# --- create summary file with header if missing or empty ---
if [[ ! -f "$SUMMARY_TXT" || ! -s "$SUMMARY_TXT" ]]; then
  printf "BAM\tRun\tSample\tHP_BASES\tNONHP_BASES\tHP_INDELS\tNONHP_INDELS\tHP_indel_rate\tnonHP_indel_rate\tHP_over_nonHP_ratio\n" > "$SUMMARY_TXT"
fi

# --- append one row to summary file ---
printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
  "$BAM" "$Run" "$Sample" "$HP_BASES" "$NONHP_BASES" "$HP_INDELS" "$NONHP_INDELS" \
  "$HP_RATE" "$NONHP_RATE" "$RATIO" >> "$SUMMARY_TXT"

echo "Wrote per-BAM output: $OUT_TXT" >&2
echo "Updated summary file: $SUMMARY_TXT" >&2
