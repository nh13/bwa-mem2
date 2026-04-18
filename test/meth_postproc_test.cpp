/* SPDX-License-Identifier: MIT */
#include "meth_postproc.h"
#include <cassert>
#include <cstdio>
#include <cstring>
#include <string>

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

static void test_rewrite_record_forward(void)
{
    kstring_t out = {0, 0, NULL};
    /* Aligned to f-strand; no chimera. */
    const char *in =
        "r1\t99\tfchr1\t100\t60\t100M\t=\t200\t150\t"
        "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\t"
        "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\t"
        "NM:i:0\tYS:Z:AAAAA";
    int rc = meth_rewrite_record(in, &out, /*set_as_failed*/ 0, /*no_chim*/ 0);
    CHECK(rc == 1, "rewrite_record rc");
    /* Should strip 'f', drop YS, append YD:Z:f. */
    CHECK(strstr(out.s, "\tchr1\t") != NULL, "f prefix stripped");
    CHECK(strstr(out.s, "YS:Z") == NULL, "YS stripped");
    CHECK(strstr(out.s, "YD:Z:f") != NULL, "YD:Z:f appended");
    free(out.s);
}

static void test_rewrite_record_chimera(void)
{
    kstring_t out = {0, 0, NULL};
    /* 100bp read with longest M = 30 (<44%). Must get 0x200, lose 0x2, mapq capped at 1. */
    const char *in =
        "r2\t99\trchr1\t100\t60\t30M70S\t=\t200\t0\t"
        "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN\t"
        "*\tNM:i:0";
    int rc = meth_rewrite_record(in, &out, 0, 0);
    CHECK(rc == 1, "chimera rewrite_record rc");
    CHECK(strstr(out.s, "YD:Z:r") != NULL, "YD:Z:r (from rchr1)");
    /* Assert flag field (2nd col) has 0x200 set and 0x2 cleared. */
    char *flag_start = strchr(out.s, '\t') + 1;
    int flag = atoi(flag_start);
    CHECK((flag & 0x200) != 0, "0x200 set");
    CHECK((flag & 0x2) == 0,   "0x2 cleared");
    free(out.s);
}

static void test_group_propagation(void)
{
    /* Two records with same QNAME. First is chimeric (will fail), second is clean.
     * Expected: both emerge with 0x200 set and 0x2 cleared. */
    std::string input =
        "r1\t99\tfchr1\t100\t60\t30M70S\t=\t200\t0\t"
        "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN\t*\n"
        "r1\t147\tfchr1\t200\t60\t100M\t=\t100\t-150\t"
        "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\t*\n";
    kstring_t out = {0, 0, NULL};
    int rc = meth_process_stream_from_string(input.c_str(), &out, 0, 0);
    CHECK(rc == 0, "stream process rc");

    /* Both record lines should have 0x200 and not 0x2.
     * Skip header lines (start with '@') like the injected @PG. */
    const char *p = out.s;
    int records_seen = 0;
    while (*p && records_seen < 2) {
        if (*p == '@') {
            const char *nl = strchr(p, '\n');
            p = nl ? nl + 1 : p + strlen(p);
            continue;
        }
        const char *tab = strchr(p, '\t');
        CHECK(tab != NULL, "found first tab");
        if (!tab) break;
        int flag = atoi(tab + 1);
        CHECK((flag & 0x200) != 0, "group-propagated 0x200");
        CHECK((flag & 0x2) == 0,   "group-propagated not-proper");
        const char *nl = strchr(p, '\n');
        CHECK(nl != NULL, "newline");
        if (!nl) break;
        p = nl + 1;
        records_seen += 1;
    }
    CHECK(records_seen == 2, "two records seen");
    free(out.s);
}

static void test_pg_injection(void)
{
    std::string input =
        "@HD\tVN:1.6\n"
        "@SQ\tSN:fchr1\tLN:1000\n"
        "@PG\tID:bwa-mem2\tPN:bwa-mem2\tVN:2.2.1\n"
        "r1\t4\t*\t0\t0\t*\t*\t0\t0\tAAAA\t!!!!\n";
    kstring_t out = {0, 0, NULL};
    meth_process_stream_from_string(input.c_str(), &out, 0, 0);
    CHECK(strstr(out.s, "\n@PG\tID:bwa-mem2-meth\t") != NULL, "bwa-mem2-meth @PG injected");
    /* Must be between last header line and first record. */
    const char *pg = strstr(out.s, "@PG\tID:bwa-mem2-meth");
    const char *r1 = strstr(out.s, "\nr1\t");
    CHECK(pg != NULL && r1 != NULL && pg < r1, "PG before records");
    free(out.s);
}

int main(void)
{
    test_longest_m();
    test_rewrite_sq();
    test_rewrite_passthrough();
    test_rewrite_record_forward();
    test_rewrite_record_chimera();
    test_group_propagation();
    test_pg_injection();
    if (failures > 0) { fprintf(stderr, "%d test(s) failed\n", failures); return 1; }
    fprintf(stderr, "OK\n");
    return 0;
}
