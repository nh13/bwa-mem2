/* SPDX-License-Identifier: MIT */
#include "meth_postproc.h"
#include <cstdlib>
#include <cstring>

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
