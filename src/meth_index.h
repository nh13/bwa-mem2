/* SPDX-License-Identifier: MIT */
#ifndef BWAMEM2_METH_INDEX_H
#define BWAMEM2_METH_INDEX_H

#ifdef __cplusplus
extern "C" {
#endif

/* `bwa-mem2 meth-index <in.fasta>` — builds the standard bwa-mem2 index
 * (.pac/.ann/.amb/.0123/.bwt.2bit.64) AND a BS-aware FMI over a C→T
 * projection of the reference (.meth.0123/.meth.bwt.2bit.64). The
 * resulting files are all anchored at <in.fasta>; normal `bwa-mem2 mem`
 * ignores the .meth.* files, and `bwa-mem2 meth` (future PR) loads them. */
int meth_index_main(int argc, char *argv[]);

#ifdef __cplusplus
}
#endif

#endif
