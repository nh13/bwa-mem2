#!/usr/bin/env bash
# Regression test: bwa-mem2 mem --meth end-to-end.
#
# Two layers of assertions:
#
# Layer 1 (always runs):  valid BAM emission.
#   - binary builds and runs with --meth
#   - produces uncompressed BAM readable by samtools
#   - @PG ID:bwa-mem2-meth present
#   - BGZF EOF marker at tail
#   - --set-as-failed / --do-not-penalize-chimeras parse cleanly
#
# Layer 2 (runs if pixi + bwameth.py available):  equivalence to bwameth.py.
#   Builds a bwameth c2t reference, c2t-converts reads, runs BOTH the
#   bwameth.py Python pipeline and `bwa-mem2 mem --meth` on the same
#   converted reads, and diffs structural fields + YD:Z tags + flag
#   distribution. Currently zero diff on the bwa-meth/example fixture.
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
BWAMEM2="$HERE/../../bwa-mem2"
SAMTOOLS="${SAMTOOLS:-samtools}"
BWAMETH_DIR="${BWAMETH_DIR:-$HOME/work/git/bwa-meth}"
BWAMETH_PY="$BWAMETH_DIR/bwameth.py"

if [[ ! -x "$BWAMEM2" ]]; then
    echo "ERROR: bwa-mem2 binary not found at $BWAMEM2. Run 'make arm64' first."
    exit 2
fi
if ! command -v "$SAMTOOLS" >/dev/null 2>&1; then
    echo "ERROR: samtools not found on PATH."
    exit 2
fi

cd "$HERE"

# ---------------------------------------------------------------------------
# Layer 1: BAM emission smoke test (plain index — no c2t required)
# ---------------------------------------------------------------------------

if [[ ! -f ref.fa.bwt.2bit.64 ]]; then
    "$BWAMEM2" index ref.fa >/dev/null 2>&1
fi

"$BWAMEM2" mem --meth -t 2 ref.fa t_R1.fastq.gz 2>/dev/null > /tmp/meth_test.bam

EXPECT_EOF="1f8b08040000000000ff0600424302001b0003000000000000000000"
ACTUAL_EOF="$(tail -c 28 /tmp/meth_test.bam | od -An -v -t x1 | tr -d ' \n')"
if [[ "${ACTUAL_EOF%$'\n'}" != "${EXPECT_EOF}" ]]; then
    echo "FAIL: BGZF EOF marker mismatch (actual=$ACTUAL_EOF)"; exit 1
fi

HDR="$("$SAMTOOLS" view -H /tmp/meth_test.bam 2>&1)"
if echo "$HDR" | grep -qi 'truncated\|EOF marker is absent'; then
    echo "FAIL: samtools reports truncated BAM"; echo "$HDR"; exit 1
fi
if ! echo "$HDR" | grep -q 'ID:bwa-mem2-meth'; then
    echo "FAIL: @PG ID:bwa-mem2-meth missing"; exit 1
fi

TOTAL="$("$SAMTOOLS" view -c /tmp/meth_test.bam 2>/dev/null)"
if [[ "$TOTAL" -lt 1 ]]; then echo "FAIL: zero records in output BAM"; exit 1; fi

"$BWAMEM2" mem --meth --set-as-failed f --do-not-penalize-chimeras \
    ref.fa t_R1.fastq.gz 2>/dev/null > /tmp/meth_test2.bam || true
if [[ ! -s /tmp/meth_test2.bam ]]; then
    echo "FAIL: --set-as-failed + --do-not-penalize-chimeras produced empty output"
    exit 1
fi

echo "OK layer 1: bwa-mem2 mem --meth (records=$TOTAL, BGZF-EOF ok, @PG bwa-mem2-meth ok)"

# ---------------------------------------------------------------------------
# Layer 2: bwa-meth equivalence (requires pixi + bwameth.py + toolshed)
# ---------------------------------------------------------------------------

if ! command -v pixi >/dev/null 2>&1; then
    echo "SKIP layer 2: pixi not on PATH"
    exit 0
fi
if [[ ! -f "$BWAMETH_PY" ]]; then
    echo "SKIP layer 2: bwameth.py not found at $BWAMETH_PY (set BWAMETH_DIR)"
    exit 0
fi
if [[ ! -f "$HERE/pyproject.toml" ]]; then
    echo "SKIP layer 2: pixi env not initialized in $HERE (run 'pixi init && pixi add toolshed')"
    exit 0
fi

export PATH="$(cd "$HERE/../.." && pwd):$PATH"

if [[ ! -f "$HERE/ref.fa.bwameth.c2t.0123" ]]; then
    (cd "$HERE" && pixi run python3 "$BWAMETH_PY" index-mem2 ref.fa >/dev/null 2>&1)
fi

pixi run python3 "$BWAMETH_PY" --reference ref.fa t_R1.fastq.gz t_R2.fastq.gz \
    2>/dev/null > /tmp/meth_oracle.sam

pixi run python3 "$BWAMETH_PY" c2t t_R1.fastq.gz t_R2.fastq.gz 2>/dev/null \
    > /tmp/meth_c2t.fq
"$BWAMEM2" mem --meth -CM -p -T 40 -B 2 -L 10 -U 100 -t 4 \
    ref.fa.bwameth.c2t /tmp/meth_c2t.fq 2>/dev/null > /tmp/meth_mine.bam
"$SAMTOOLS" view /tmp/meth_mine.bam 2>/dev/null > /tmp/meth_mine.sam

grep -v '^@' /tmp/meth_oracle.sam > /tmp/meth_oracle_records.sam

norm() {
    awk -F'\t' 'BEGIN{OFS="\t"} { if ($7 == "=") $7 = $3; print }' "$1" \
        | sort -k1,1 -k2,2n
}

ONE="$(diff <(norm /tmp/meth_mine.sam | cut -f1-9) \
             <(norm /tmp/meth_oracle_records.sam | cut -f1-9) | wc -l | tr -d ' ')"
if [[ "$ONE" != "0" ]]; then
    echo "FAIL layer 2: structural diff in cols 1-9 ($ONE lines)"
    diff <(norm /tmp/meth_mine.sam | cut -f1-9) <(norm /tmp/meth_oracle_records.sam | cut -f1-9) | head -20
    exit 1
fi

TWO="$(diff <(norm /tmp/meth_mine.sam | cut -f6) \
             <(norm /tmp/meth_oracle_records.sam | cut -f6) | wc -l | tr -d ' ')"
if [[ "$TWO" != "0" ]]; then echo "FAIL layer 2: CIGAR diff ($TWO lines)"; exit 1; fi

for d in f r; do
    MINE_YD="$(grep -c "YD:Z:$d" /tmp/meth_mine.sam || true)"
    ORACLE_YD="$(grep -c "YD:Z:$d" /tmp/meth_oracle_records.sam || true)"
    if [[ "$MINE_YD" != "$ORACLE_YD" ]]; then
        echo "FAIL layer 2: YD:Z:$d count mismatch (mine=$MINE_YD oracle=$ORACLE_YD)"
        exit 1
    fi
done

MINE_N="$(wc -l < /tmp/meth_mine.sam | tr -d ' ')"
ORACLE_N="$(wc -l < /tmp/meth_oracle_records.sam | tr -d ' ')"
if [[ "$MINE_N" != "$ORACLE_N" ]]; then
    echo "FAIL layer 2: record count mismatch (mine=$MINE_N oracle=$ORACLE_N)"
    exit 1
fi

MINE_F="$(grep -c YD:Z:f /tmp/meth_mine.sam || true)"
MINE_R="$(grep -c YD:Z:r /tmp/meth_mine.sam || true)"
echo "OK layer 2: bwa-mem2 mem --meth matches bwameth.py (records=$MINE_N, YD:Z:f=$MINE_F YD:Z:r=$MINE_R)"
