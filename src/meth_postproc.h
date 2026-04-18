/* SPDX-License-Identifier: MIT */
#ifndef BWAMEM2_METH_POSTPROC_H
#define BWAMEM2_METH_POSTPROC_H

#include "kstring.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Longest run of M, =, or X CIGAR ops. Returns 0 on "*" or parse error. */
int meth_longest_m(const char *cigar);

/* Append rewritten header line (with trailing newline) to *out.
 * Returns 1 if emitted, 0 if dropped (r-strand @SQ), -1 on error.
 * `line` must be null-terminated and exclude any trailing newline. */
int meth_rewrite_header_line(const char *line, kstring_t *out);

/* Main entry. argc/argv start at "meth-postproc" (argv[0] == "meth-postproc"). */
int meth_postproc_main(int argc, char *argv[]);

#ifdef __cplusplus
}
#endif

#endif
