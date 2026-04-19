/* SPDX-License-Identifier: MIT
 *
 * bwa-mem2 meth — thin wrapper over `main_mem` that auto-injects `--meth`
 * and `--meth-index` so the user gets native BS-aware alignment + BAM
 * post-processing without having to remember two long options.
 *
 * Phase B (this file): pure argv rewrite. The FMI prefix split (loading
 * the BS-aware `.meth.*` FMI alongside the original-alphabet `.pac`) is
 * done inside `main_mem` via `opt->meth_dual_index`. No seed/extension
 * changes yet — those land in Phase C.
 */

#include "meth_index.h"
#include "fastmap.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>

int meth_align_main(int argc, char *argv[])
{
    /* Build a new argv with `--meth-index` prepended as the first flag.
     * (`--meth-index` implies `--meth`, so one injection suffices.) We
     * allocate a (argc + 2) slot argv so the final NULL sentinel is kept.
     */
    int new_argc = argc + 1;
    char **new_argv = (char **)calloc(new_argc + 1, sizeof(char *));
    if (new_argv == NULL) {
        fprintf(stderr, "ERROR: meth: out of memory building argv\n");
        return 1;
    }

    new_argv[0] = argv[0];
    new_argv[1] = (char *)"--meth-index";
    for (int i = 1; i < argc; ++i) new_argv[i + 1] = argv[i];
    new_argv[new_argc] = NULL;

    int rc = main_mem(new_argc, new_argv);
    free(new_argv);
    return rc;
}
