/* SPDX-License-Identifier: MIT */
#ifndef BWAMEM2_METH_POSTPROC_H
#define BWAMEM2_METH_POSTPROC_H

/* Shared helpers for `bwa-mem2 mem --meth`.
 *
 * Historically this module held a SAM-text post-processor (meth-postproc
 * subcommand). That path has been removed in favor of direct BAM emission
 * from mem_aln_t in bwamem.cpp (see mem_aln2bam). This header retains only
 * the small, unit-testable primitives reused by the direct path. */

#ifdef __cplusplus
extern "C" {
#endif

/* Longest run of M, =, or X CIGAR ops. Returns 0 on "*" or parse error.
 * Accepts a textual CIGAR (e.g. "30M70S"); the binary-CIGAR caller in
 * bwamem.cpp computes the same value directly from the uint32_t array. */
int meth_longest_m(const char *cigar);

#ifdef __cplusplus
}
#endif

#endif
