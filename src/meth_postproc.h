/* SPDX-License-Identifier: MIT */
#ifndef BWAMEM2_METH_POSTPROC_H
#define BWAMEM2_METH_POSTPROC_H

#ifdef __cplusplus
extern "C" {
#endif

/* Longest run of M, =, or X CIGAR ops. Returns 0 on "*" or parse error. */
int meth_longest_m(const char *cigar);

/* Main entry. argc/argv start at "meth-postproc" (argv[0] == "meth-postproc"). */
int meth_postproc_main(int argc, char *argv[]);

#ifdef __cplusplus
}
#endif

#endif
