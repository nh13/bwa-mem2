/* SPDX-License-Identifier: MIT */
#include "meth_postproc.h"
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

int main(void)
{
    test_longest_m();
    if (failures > 0) { fprintf(stderr, "%d test(s) failed\n", failures); return 1; }
    fprintf(stderr, "OK\n");
    return 0;
}
