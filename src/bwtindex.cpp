/* The MIT License

   Copyright (c) 2008 Genome Research Ltd (GRL).

   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:

   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.

   Modified Copyright (C) 2019 Intel Corporation, Heng Li.
   Contacts: Vasimuddin Md <vasimuddin.md@intel.com>; Sanchit Misra <sanchit.misra@intel.com>;
   Heng Li <hli@jimmy.harvard.edu> 
*/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <time.h>
#include <getopt.h>
#include <limits.h>
#include <zlib.h>
#include "bntseq.h"
#include "bwa.h"
#include "bwt.h"
#include "utils.h"
#include "FMI_search.h"
#include "kseq.h"

KSEQ_DECLARE(gzFile)

/* --- bwameth-style c2t doubled-reference builder ---
 *
 * Emits `<fa>.bwameth.c2t` with two contigs per input chromosome:
 *   >r<name>    original seq with G→A applied (reverse-strand BS target)
 *   >f<name>    original seq with C→T applied (forward-strand BS target)
 * Then runs the standard bwa_idx_build on that file. Byte-identical to
 * the FASTA emitted by bwameth.py's `index-mem2` subcommand, so the
 * resulting FMI works with any existing bwameth.py tooling. */
static int meth_index_c2t_build(const char *fa)
{
    char out_fa[PATH_MAX];
    int n = snprintf(out_fa, sizeof(out_fa), "%s.bwameth.c2t", fa);
    if (n <= 0 || (size_t)n >= sizeof(out_fa)) {
        fprintf(stderr, "ERROR: reference path too long\n");
        return 1;
    }

    gzFile in = gzopen(fa, "r");
    if (in == NULL) {
        fprintf(stderr, "ERROR: cannot open %s\n", fa);
        return 2;
    }
    FILE *out = fopen(out_fa, "w");
    if (out == NULL) {
        fprintf(stderr, "ERROR: cannot open %s for writing\n", out_fa);
        gzclose(in);
        return 3;
    }

    fprintf(stderr, "[bwa_index:--meth] writing %s ...\n", out_fa);
    kseq_t *seq = kseq_init(in);
    int64_t total_bases = 0, n_seqs = 0;
    while (kseq_read(seq) >= 0) {
        /* Upper-case in place — matches bwameth.py's fasta_iter which upper-
         * cases each sequence before applying the substitutions. Without
         * this, soft-masked reference regions would keep lowercase bases
         * and our FASTA would not be byte-identical to bwameth.py's. */
        for (size_t i = 0; i < seq->seq.l; ++i) {
            char c = seq->seq.s[i];
            if (c >= 'a' && c <= 'z') seq->seq.s[i] = (char)(c - 'a' + 'A');
        }

        /* >r<name>  with G→A */
        fprintf(out, ">r%s\n", seq->name.s);
        for (size_t i = 0; i < seq->seq.l; ++i) {
            char c = seq->seq.s[i];
            if (c == 'G') c = 'A';
            fputc(c, out);
            if (((i + 1) % 100) == 0) fputc('\n', out);
        }
        if (seq->seq.l % 100 != 0) fputc('\n', out);

        /* >f<name>  with C→T */
        fprintf(out, ">f%s\n", seq->name.s);
        for (size_t i = 0; i < seq->seq.l; ++i) {
            char c = seq->seq.s[i];
            if (c == 'C') c = 'T';
            fputc(c, out);
            if (((i + 1) % 100) == 0) fputc('\n', out);
        }
        if (seq->seq.l % 100 != 0) fputc('\n', out);

        total_bases += (int64_t)seq->seq.l;
        ++n_seqs;
    }
    kseq_destroy(seq);
    gzclose(in);
    if (fclose(out) != 0) {
        fprintf(stderr, "ERROR: failed to close %s\n", out_fa);
        return 4;
    }
    fprintf(stderr, "[bwa_index:--meth] emitted %lld seqs, %lld bases (each as r+f = %lld bp of c2t text)\n",
            (long long)n_seqs, (long long)total_bases, (long long)(2 * total_bases));

    /* Build the FMI over the doubled c2t FASTA, with prefix == the c2t file
     * (so downstream `bwa-mem2 mem --meth` can find it at `<fa>.bwameth.c2t.*`). */
    if (bwa_idx_build(out_fa, out_fa) != 0) {
        fprintf(stderr, "ERROR: bwa_idx_build failed on %s\n", out_fa);
        return 5;
    }
    return 0;
}

int bwa_index(int argc, char *argv[]) // the "index" command
{
	int c;
	char *prefix = 0;
	int meth = 0;
	static struct option long_opts[] = {
		{"meth", no_argument, 0, 1000},
		{0, 0, 0, 0}
	};
	while ((c = getopt_long(argc, argv, "p:", long_opts, NULL)) >= 0) {
		if (c == 'p') prefix = optarg;
		else if (c == 1000) meth = 1;
		else return 1;
	}

	if (optind + 1 > argc) {
		fprintf(stderr, "Usage: bwa-mem2 index [-p prefix] [--meth] <in.fasta>\n");
		fprintf(stderr, "\n"
		        "  -p STR    output prefix (default: <in.fasta>)\n"
		        "  --meth    build a bwameth-style doubled c2t reference + FMI.\n"
		        "            Writes <in.fasta>.bwameth.c2t and the FMI alongside it.\n"
		        "            Use with `bwa-mem2 mem --meth <in.fasta> R1.fq [R2.fq]`.\n");
		return 1;
	}
	if (meth) {
		if (prefix != 0) {
			fprintf(stderr, "ERROR: --meth does not accept -p (prefix is <in.fasta>.bwameth.c2t)\n");
			return 1;
		}
		return meth_index_c2t_build(argv[optind]);
	}
	if (prefix == 0) prefix = argv[optind];
	bwa_idx_build(argv[optind], prefix);
	return 0;
}

int bwa_idx_build(const char *fa, const char *prefix)
{
	extern void bwa_pac_rev_core(const char *fn, const char *fn_rev);

	clock_t t;
	int64_t l_pac;

	{ // nucleotide indexing
		gzFile fp = xzopen(fa, "r");
		t = clock();
		fprintf(stderr, "[bwa_index] Pack FASTA... ");
		l_pac = bns_fasta2bntseq(fp, prefix, 1);
		fprintf(stderr, "%.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
		err_gzclose(fp);
        FMI_search *fmi = new FMI_search(prefix);
        fmi->build_index();
        delete fmi;
	}
	return 0;
}
