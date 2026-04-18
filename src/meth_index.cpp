/* SPDX-License-Identifier: MIT */
#include "meth_index.h"
#include "FMI_search.h"
#include "bntseq.h"
#include "bwa.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>

extern "C" int bwa_idx_build(const char *fa, const char *prefix);

int meth_index_main(int argc, char *argv[])
{
    const char *prefix = NULL;
    int opt;
    while ((opt = getopt(argc, argv, "p:")) >= 0) {
        if (opt == 'p') prefix = optarg;
        else return 1;
    }
    if (optind + 1 > argc) {
        fprintf(stderr, "Usage: bwa-mem2 meth-index [-p prefix] <in.fasta>\n");
        fprintf(stderr, "\n"
                "Builds the standard bwa-mem2 index for <in.fasta> AND an\n"
                "additional BS-aware FMI over a C->T-projected version of the\n"
                "reference. Produces these files alongside the input FASTA:\n"
                "  <prefix>.pac          (original alphabet, shared)\n"
                "  <prefix>.ann\n"
                "  <prefix>.amb\n"
                "  <prefix>.0123         (normal 4-letter binary ref)\n"
                "  <prefix>.bwt.2bit.64  (normal FM-index)\n"
                "  <prefix>.meth.0123        (C->T projected binary ref)\n"
                "  <prefix>.meth.bwt.2bit.64 (BS-aware FM-index)\n");
        return 1;
    }

    const char *fa = argv[optind];
    if (prefix == NULL) prefix = fa;

    /* Phase 1: emit .pac / .ann / .amb / .0123 / .bwt.2bit.64 exactly as
     * `bwa-mem2 index` does. Reuses bwa_idx_build so the normal index is
     * bit-for-bit identical to what the non-meth path produces. */
    if (bwa_idx_build(fa, prefix) != 0) {
        fprintf(stderr, "ERROR: meth-index: bwa_idx_build failed\n");
        return 2;
    }

    /* Phase 2: build the BS-aware FMI over a C->T projection of the .pac.
     * Emits <prefix>.meth.0123 and <prefix>.meth.bwt.2bit.64. */
    FMI_search *fmi = new FMI_search(prefix);
    int rc = fmi->build_index_bs(".meth");
    delete fmi;
    if (rc != 0) {
        fprintf(stderr, "ERROR: meth-index: build_index_bs failed (rc=%d)\n", rc);
        return 3;
    }

    fprintf(stderr, "[meth-index] OK: wrote %s.meth.0123 and %s.meth.bwt.2bit.64\n",
            prefix, prefix);
    return 0;
}
