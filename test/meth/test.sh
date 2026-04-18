#!/usr/bin/env bash
# Regression test: bwa-mem2 meth-postproc must produce a deterministic, expected
# output against a hand-crafted SAM fixture that exercises:
#   - header rewriting: f/r chrom prefix strip, r-strand @SQ drop
#   - YS:Z tag strip
#   - YD:Z:{f,r} emission from chrom prefix
#   - chimera QC heuristic (longest M < 44% of read length)
#   - QC-fail propagation across QNAME pair groups
#   - MAPQ cap to 1 on chimera detection
#   - unmapped-pair and mate-unmapped passthrough
#   - RNEXT f/r prefix stripping (not '=' or '*')
#
# If bwa-mem2's alignment behavior ever changes (Plan 2), regenerate the fixture:
#     ./bwa-mem2 meth-postproc < fixtures/raw.sam > fixtures/expected.sam
# after verifying correctness of the new output against bwameth.py's as_bam().
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
BWAMEM2="$HERE/../../bwa-mem2"

if [[ ! -x "$BWAMEM2" ]]; then
    echo "ERROR: bwa-mem2 binary not found at $BWAMEM2. Run 'make arm64' first."
    exit 2
fi

cd "$HERE"

"$BWAMEM2" meth-postproc < fixtures/raw.sam > /tmp/meth_actual.sam

# Normalize the dynamic @PG CL field (contains absolute path to binary).
normalize() {
    sed -E 's|CL:[^[:space:]]*meth-postproc|CL:meth-postproc|g' "$1"
}

if ! diff <(normalize fixtures/expected.sam) <(normalize /tmp/meth_actual.sam) > /tmp/meth_diff.txt; then
    echo "REGRESSION: differences found vs fixtures/expected.sam:"
    head -40 /tmp/meth_diff.txt
    exit 1
fi
echo "OK: bwa-mem2 meth-postproc matches fixtures/expected.sam"
