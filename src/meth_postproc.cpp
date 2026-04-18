/* SPDX-License-Identifier: MIT */
#include "meth_postproc.h"
#include <cstdio>
#include <cstring>
#include <cstdlib>

int meth_rewrite_header_line(const char *line, kstring_t *out)
{
    if (line == NULL || out == NULL) return -1;
    if (line[0] != '@') return -1;
    /* Only @SQ needs rewriting; everything else passes through verbatim. */
    if (strncmp(line, "@SQ\t", 4) != 0) {
        kputs(line, out);
        kputc('\n', out);
        return 1;
    }
    /* Find SN: tag */
    const char *sn = strstr(line, "\tSN:");
    if (sn == NULL) { kputs(line, out); kputc('\n', out); return 1; }
    const char *name = sn + 4;            /* after "\tSN:" */
    char prefix = name[0];
    if (prefix == 'r') return 0;          /* drop reverse-strand SQ */
    if (prefix != 'f') {                  /* unexpected; pass through unchanged */
        kputs(line, out);
        kputc('\n', out);
        return 1;
    }
    /* Emit prefix + "\tSN:" + (name without leading 'f') + suffix */
    kputsn(line, (int)(name - line), out); /* includes the "\tSN:" */
    kputs(name + 1, out);                  /* skip the 'f' */
    kputc('\n', out);
    return 1;
}

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
