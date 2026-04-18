/* SPDX-License-Identifier: MIT */
#include "meth_postproc.h"
#include <cstdio>
#include <cstring>
#include <cstdlib>

int meth_longest_m(const char *cigar)
{
    if (cigar == NULL || cigar[0] == '\0' || (cigar[0] == '*' && cigar[1] == '\0')) return 0;
    int longest = 0;
    const char *p = cigar;
    while (*p) {
        char *end;
        long n = strtol(p, &end, 10);
        if (end == p) return 0; /* malformed */
        char op = *end;
        if (op == '\0') return 0;
        if (op == 'M' || op == '=' || op == 'X') {
            if ((int)n > longest) longest = (int)n;
        }
        p = end + 1;
    }
    return longest;
}

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
