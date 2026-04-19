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
 * ignores the .meth.* files, and `bwa-mem2 meth` loads them. */
int meth_index_main(int argc, char *argv[]);

/* `bwa-mem2 meth <ref.fa> R1.fq [R2.fq]` — thin wrapper that invokes the
 * `mem` pipeline with `--meth` and `--meth-index` auto-injected, so the
 * user gets native BS-aware seeding + post-processing without remembering
 * two long options. Returns the same exit code as `main_mem`. */
int meth_align_main(int argc, char *argv[]);

#ifdef __cplusplus
}
#endif

#endif
