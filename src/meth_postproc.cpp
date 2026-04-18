/* SPDX-License-Identifier: MIT */
#include "meth_postproc.h"
#include <cstdio>
#include <cstring>
#include <cstdlib>

/* Split the line into fields in-place on a mutable buffer. Returns number of
 * fields, or -1 if max_fields would overflow. */
static int meth_split_fields(char *line, char **fields, int max_fields)
{
    int n = 0;
    if (max_fields <= 0) return -1;
    fields[n++] = line;
    for (char *p = line; *p; ++p) {
        if (*p == '\t') {
            *p = '\0';
            if (n >= max_fields) return -1;
            fields[n++] = p + 1;
        }
    }
    return n;
}

int meth_rewrite_record(const char *line, kstring_t *out,
                        char set_as_failed, int no_chim)
{
    if (line == NULL || out == NULL) return -1;
    size_t len = strlen(line);
    char *buf = (char *)malloc(len + 1);
    if (buf == NULL) return -1;
    memcpy(buf, line, len + 1);

    char *fields[4096];
    int n = meth_split_fields(buf, fields, 4096);
    if (n < 11) { free(buf); return -1; }

    int flag     = atoi(fields[1]);
    char *rname  =       fields[2];
    int mapq     = atoi(fields[4]);
    char *cigar  =       fields[5];
    char *rnext  =       fields[6];
    const char *seq = fields[9];
    int seq_len  = (seq[0] == '*' && seq[1] == '\0') ? 0 : (int)strlen(seq);

    char direction = 0;
    if (rname[0] == 'f' || rname[0] == 'r') {
        direction = rname[0];
        rname = rname + 1;
    }
    if (rnext[0] == 'f' || rnext[0] == 'r') {
        rnext = rnext + 1;
    } /* '=' and '*' stay as-is */

    int mapped = !(flag & 0x4) && direction != 0;
    if (mapped && set_as_failed != 0 && set_as_failed == direction) flag |= 0x200;
    if (mapped && !no_chim) {
        int lm = meth_longest_m(cigar);
        /* threshold: lm < 0.44 * seq_len, integer-safe: 100*lm < 44*seq_len */
        if (seq_len > 0 && 100 * lm < 44 * seq_len) {
            flag |= 0x200;
            flag &= ~0x2;
            if (mapq > 1) mapq = 1;
        }
    }

    /* Emit: qname flag rname pos mapq cigar rnext pnext tlen seq qual */
    ksprintf(out, "%s\t%d\t%s\t%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s",
             fields[0], flag, rname, fields[3], mapq, cigar,
             rnext, fields[7], fields[8], fields[9], fields[10]);

    for (int i = 11; i < n; ++i) {
        if (strncmp(fields[i], "YS:Z:", 5) == 0) continue; /* drop incoming YS */
        kputc('\t', out);
        kputs(fields[i], out);
    }
    if (mapped) {
        kputs("\tYD:Z:", out);
        kputc(direction, out);
    }
    kputc('\n', out);

    free(buf);
    return 1;
}

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
