#!/usr/bin/env bash
# download_reference.sh — robust reference FASTA fetcher (with retries, checksum, indexing)
# Usage:
#   bash download_reference.sh [REF_DIR] [REF_URL] [OUT_BASENAME] [SHA256_OPT]
#
# Defaults:
#   REF_DIR      = ref
#   REF_URL      = https://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/chr8.fa.gz
#   OUT_BASENAME = hg19_chr8
#
# Examples:
#   bash download_reference.sh
#   bash download_reference.sh ref https://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/chr8.fa.gz hg19_chr8
#   bash download_reference.sh ref <URL> <BASENAME> <SHA256SUM>
#
# Notes:
# - Uses curl with --http1.1 to avoid certain HTTP/2 PROTOCOL_ERRORs (e.g., Dropbox).
# - If URL ends with .gz, file is gunzipped to REF_DIR/OUT_BASENAME.fa.
# - Builds samtools (.fai) and BWA indexes.
# - Requires: curl, gunzip, samtools, bwa

set -euo pipefail
IFS=$'\n\t'

REF_DIR="${1:-ref}"
REF_URL="${2:-https://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/chr8.fa.gz}"
OUT_BASENAME="${3:-hg19_chr8}"
SHA256_OPT="${4:-}"

command -v curl >/dev/null 2>&1 || { echo "ERROR: curl not found in PATH." >&2; exit 127; }
command -v gunzip >/dev/null 2>&1 || { echo "ERROR: gunzip not found." >&2; exit 127; }
command -v samtools >/dev/null 2>&1 || { echo "ERROR: samtools not found in PATH." >&2; exit 127; }
command -v bwa >/dev/null 2>&1 || { echo "ERROR: bwa not found in PATH." >&2; exit 127; }

mkdir -p "${REF_DIR}"

TMP="${REF_DIR}/${OUT_BASENAME}.fa.gz.part"
OUT_GZ="${REF_DIR}/${OUT_BASENAME}.fa.gz"
OUT_FA="${REF_DIR}/${OUT_BASENAME}.fa"

echo "[download_reference] Downloading: ${REF_URL}"
# --fail for HTTP errors, --retry for transient failures, --http1.1 to dodge HTTP/2 issues on some hosts
if ! curl -L --fail --retry 5 --retry-delay 3 --http1.1 -o "${TMP}" "${REF_URL}"; then
  echo "[download_reference] ERROR: download failed from ${REF_URL}" >&2
  echo "  Tips: If using Dropbox, ensure you have a stable 'dl=1' share link (not a time-limited token)." >&2
  exit 1
fi

# If the URL wasn't gzipped, rename TMP appropriately
if [[ "${REF_URL}" != *.gz ]]; then
  mv "${TMP}" "${OUT_FA}.part"
  OUTDOWN="${OUT_FA}.part"
else
  mv "${TMP}" "${OUT_GZ}"
  OUTDOWN="${OUT_GZ}"
fi

# Optional checksum
if [[ -n "${SHA256_OPT}" ]]; then
  echo "${SHA256_OPT}  ${OUTDOWN}" | sha256sum -c -
fi

# Decompress if needed
if [[ "${OUTDOWN}" == *.gz ]]; then
  echo "[download_reference] Decompressing to ${OUT_FA}"
  gunzip -c "${OUTDOWN}" > "${OUT_FA}"
  rm -f "${OUTDOWN}"
else
  mv "${OUTDOWN}" "${OUT_FA}"
fi

# Basic sanity check
if [[ ! -s "${OUT_FA}" ]]; then
  echo "[download_reference] ERROR: reference file is empty: ${OUT_FA}" >&2
  exit 1
fi
if ! head -c 1 "${OUT_FA}" >/dev/null 2>&1; then
  echo "[download_reference] ERROR: cannot read ${OUT_FA}" >&2
  exit 1
fi

echo "[download_reference] Indexing with samtools and BWA"
samtools faidx "${OUT_FA}"
bwa index "${OUT_FA}"

echo "[download_reference] Done. Reference ready at ${OUT_FA}"
