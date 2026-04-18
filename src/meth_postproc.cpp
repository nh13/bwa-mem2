/* SPDX-License-Identifier: MIT */
#include "meth_postproc.h"
#include <cstdio>
#include <cstring>

static int usage(FILE *fp)
{
    fprintf(fp,
        "Usage: bwa-mem2 meth-postproc [options] < in.sam > out.sam\n"
        "\n"
        "Options:\n"
        "  --set-as-failed {f,r}       Flag reads aligned to this strand as QC-fail (0x200)\n"
        "  --do-not-penalize-chimeras  Skip the longest-match <44%% chimera heuristic\n"
        "  -h, --help                  Show this message\n");
    return 1;
}

int meth_postproc_main(int argc, char *argv[])
{
    (void)argc; (void)argv;
    return usage(stderr);
}
