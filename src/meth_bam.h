/* SPDX-License-Identifier: MIT */
#ifndef BWAMEM2_METH_BAM_H
#define BWAMEM2_METH_BAM_H

#include <stdint.h>
#include "htslib/sam.h"
#include "bwa.h"
#include "bwamem.h"
#include "bntseq.h"

/*
 * Native BAM emission for `bwa-mem2 mem --meth`.
 *
 * bwa-meth's c2t reference holds each chrom twice (fchr/rchr). This module
 * builds a consolidation map from the bntseq_t, emits a deduped BAM header,
 * and provides mem_aln_t → bam1_t conversion with meth-specific transforms
 * (chrom rewrite, YD:Z tag, chimera QC on CIGAR, pair-level QC propagation).
 *
 * htslib handles all binary BAM encoding, BGZF framing, and I/O.
 */

#ifdef __cplusplus
extern "C" {
#endif

/* --- Chrom map: internal (fchr/rchr) -> output (chr) ---------------- */

typedef struct meth_chrom_map_s {
    int      n_internal;    /* matches bns->n_seqs */
    int      n_output;      /* after f/r collapse */
    int     *out_tid;       /* length n_internal; index into output_names */
    char    *direction;     /* length n_internal; 'f', 'r', or 0 */
    char   **output_names;  /* length n_output; stripped names (NUL-term) */
    int64_t *output_lens;   /* length n_output */
} meth_chrom_map_t;

meth_chrom_map_t *meth_chrom_map_build_from_bns(const bntseq_t *bns);
void              meth_chrom_map_free(meth_chrom_map_t *m);

/* --- BAM writer lifecycle ------------------------------------------- */

typedef struct meth_bam_writer_s {
    htsFile          *fp;       /* htslib output handle (wraps BGZF) */
    sam_hdr_t        *hdr;      /* header with consolidated @SQ list */
    meth_chrom_map_t *cmap;     /* non-owning pointer; caller frees */
} meth_bam_writer_t;

/* Open a BAM writer. `path_or_dash` is "-" for stdout or a filename.
 * Compression: "wb"=deflated, "wu"=uncompressed (BGZF level 0). We default
 * to "wb0" (synonym for wu) for now — `--meth` emits uncompressed BAM as
 * specified. `pg_line` is the CL/VN for the bwa-mem2-meth @PG entry.
 * `bwa_pg` is the existing bwa-mem2 @PG line (may be NULL).
 * Returns NULL on failure. */
meth_bam_writer_t *meth_bam_writer_open(const char *path_or_dash,
                                        meth_chrom_map_t *cmap,
                                        const char *bwa_pg,
                                        const char *meth_pg_cl);

/* Write one bam1_t. Returns 0 on success, -1 on error. */
int meth_bam_writer_write(meth_bam_writer_t *w, bam1_t *b);

/* Close the writer and flush the BGZF EOF marker. Frees internal hdr and
 * htsFile; does NOT free the cmap. Returns 0 on success, -1 on error. */
int meth_bam_writer_close(meth_bam_writer_t *w);

/* --- mem_aln_t -> bam1_t --------------------------------------------- */

/*
 * Convert one alignment from bwa-mem2's internal representation into a
 * bam1_t. The caller owns `b` (allocate via bam_init1 once and reuse, or
 * allocate per-record).
 *
 * This function performs these meth transforms:
 *   - chrom rewrite: p->rid → cmap->out_tid[p->rid]
 *   - strand detection: chrom prefix ('f'/'r') becomes direction[p->rid]
 *   - emits YD:Z:{f,r} aux tag for mapped records when direction != 0
 *   - chimera QC on CIGAR: if longest M/=/X run < 44% of l_seq,
 *     OR 0x200 into flag, clear 0x2, cap mapq at 1 (unless no_chim)
 *   - honors set_as_failed = 'f' or 'r' to force 0x200 on that strand
 *
 * Returns 0 on success, -1 on error. */
int meth_mem_aln_to_bam(bam1_t *b,
                        const mem_opt_t *opt, const bntseq_t *bns,
                        const bseq1_t *s, int n_alns,
                        const mem_aln_t *list, int which,
                        const mem_aln_t *m_,
                        const meth_chrom_map_t *cmap);

/* --- QNAME-group QC-fail propagation ---------------------------------- */

/* Propagate 0x200 across an array of bam1_t*. If any record has 0x200,
 * set 0x200 and clear 0x2 on all others. Leaves mapq untouched (the
 * per-record chimera heuristic handled its own mapq cap). */
void meth_bam_group_propagate_qcfail(bam1_t **group, int n);

#ifdef __cplusplus
}
#endif

#endif
