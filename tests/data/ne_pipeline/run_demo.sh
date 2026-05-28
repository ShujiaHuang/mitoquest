#!/usr/bin/env bash
#
# End-to-end smoke test for `mitoquest trans-prep` and `mitoquest ne-estimate`
# using the synthetic mother-child mtDNA cohort in this directory.
#
# Expected outcome (true Ne = 5, seed 20260528):
#   * trans-prep matches all 20 trios and emits 240 PASS rows.
#   * ne-estimate reports Ne ~~ 5 with a tight CI.
#   * Kimura cross-check reports Ne_kimura in roughly the same ballpark.
#
# Usage:
#     bash run_demo.sh                  # uses ../../../bin/mitoquest
#     MITOQUEST=/path/to/mitoquest bash run_demo.sh
#
# Author: Shujia Huang (hshujia@qq.com)
# Date:   2026-05-28

set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${HERE}/../../.." && pwd)"
MITOQUEST="${MITOQUEST:-${REPO_ROOT}/bin/mitoquest}"

if [[ ! -x "${MITOQUEST}" ]]; then
    echo "ERROR: mitoquest binary not found or not executable: ${MITOQUEST}" >&2
    echo "Build it first with: (cd ${REPO_ROOT} && cmake --build build -j8)" >&2
    exit 1
fi

VCF="${HERE}/cohort.vcf"
FAM="${HERE}/cohort.fam"
TSV_OUT="${HERE}/cohort.trans_prep.run.tsv"
JSON_OUT="${HERE}/cohort.ne.run.json"

echo "==> trans-prep:  ${VCF}  +  ${FAM}  ->  ${TSV_OUT}"
"${MITOQUEST}" trans-prep \
    -v "${VCF}" \
    -f "${FAM}" \
    -d 100 \
    -o "${TSV_OUT}"

echo
echo "==> ne-estimate (with --cross-check kimura):  ${TSV_OUT}  ->  ${JSON_OUT}"
"${MITOQUEST}" ne-estimate \
    -i "${TSV_OUT}" \
    --min-vaf 0.10 --max-vaf 0.90 \
    --min-ne  1    --max-ne  50 \
    --cross-check kimura \
    -t 4 \
    -o "${JSON_OUT}"

echo
echo "==> JSON output:"
cat "${JSON_OUT}"
echo

echo "Demo finished.  Generated files:"
echo "  ${TSV_OUT}"
echo "  ${JSON_OUT}"
