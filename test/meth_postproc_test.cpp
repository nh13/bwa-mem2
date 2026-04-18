/* SPDX-License-Identifier: MIT */
#include "meth_postproc.h"
#include <cassert>
#include <cstdio>
#include <cstring>

static int failures = 0;
#define CHECK(cond, msg) do { if (!(cond)) { fprintf(stderr, "FAIL: %s\n", msg); ++failures; } } while (0)

static void test_longest_m(void)
{
    CHECK(meth_longest_m("100M")      == 100, "100M -> 100");
    CHECK(meth_longest_m("50M2I50M")  == 50,  "50M2I50M -> 50 (broken by I)");
    CHECK(meth_longest_m("10S80M10S") == 80,  "soft clip ignored");
    CHECK(meth_longest_m("5H80M5H")   == 80,  "hard clip ignored");
    CHECK(meth_longest_m("30M10D70M") == 70,  "D breaks M run");
    CHECK(meth_longest_m("40=5X40=")  == 40,  "= and X count separately from M");
    CHECK(meth_longest_m("*")         == 0,   "unaligned");
    CHECK(meth_longest_m("")          == 0,   "empty");
}

static void test_rewrite_sq(void)
{
    kstring_t out = {0, 0, NULL};
    const char *keep_f = "@SQ\tSN:fchr1\tLN:249250621";
    const char *drop_r = "@SQ\tSN:rchr1\tLN:249250621";
    CHECK(meth_rewrite_header_line(keep_f, &out) == 1, "f-SQ keep");
    CHECK(strcmp(out.s, "@SQ\tSN:chr1\tLN:249250621\n") == 0, "f-SQ rewritten");
    out.l = 0; if (out.s) out.s[0] = '\0';
    CHECK(meth_rewrite_header_line(drop_r, &out) == 0, "r-SQ drop");
    free(out.s);
}

static void test_rewrite_passthrough(void)
{
    kstring_t out = {0, 0, NULL};
    const char *hd = "@HD\tVN:1.6\tSO:unsorted";
    CHECK(meth_rewrite_header_line(hd, &out) == 1, "@HD passthrough");
    CHECK(strcmp(out.s, "@HD\tVN:1.6\tSO:unsorted\n") == 0, "@HD unchanged");
    free(out.s);
}

int main(void)
{
    test_longest_m();
    test_rewrite_sq();
    test_rewrite_passthrough();
    if (failures > 0) { fprintf(stderr, "%d test(s) failed\n", failures); return 1; }
    fprintf(stderr, "OK\n");
    return 0;
}
