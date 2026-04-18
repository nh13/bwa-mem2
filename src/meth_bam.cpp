/* SPDX-License-Identifier: MIT */

/* htslib headers must come before any bwa-mem2 header that pulls in
 * bwa-mem2's kstring.h (they share the KSTRING_H include guard). */
#include "htslib/sam.h"
#include "htslib/kstring.h"

#include "meth_bam.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>

/* Private writer struct — opaque to callers via meth_bam.h. */
struct meth_bam_writer_s {
    htsFile          *fp;
    sam_hdr_t        *hdr;
    meth_chrom_map_t *cmap; /* non-owning */
};

/* Declared in src/bwa.h (without extern "C"); match that linkage here. */
extern char bwa_rg_id[256];

/* Globals (declared extern in meth_bam.h) */
meth_bam_writer_t *g_meth_bam_writer = NULL;
/* g_meth_cmap is defined in bwamem.cpp (so the worker hook can access it
 * without a link dependency on meth_bam.cpp; both files see the header). */

/* --- Allocation wrappers ------------------------------------------- */

extern "C" {
struct bam1_t *meth_bam_alloc(void) { return bam_init1(); }
void meth_bam_free(struct bam1_t *b) { if (b) bam_destroy1(b); }
}

/* ------------------------------------------------------------------- */
/* Chrom map                                                            */
/* ------------------------------------------------------------------- */

static char *meth_strdup_local(const char *s)
{
    size_t n = strlen(s) + 1;
    char *p = (char *)malloc(n);
    if (p) memcpy(p, s, n);
    return p;
}

meth_chrom_map_t *meth_chrom_map_build_from_bns(const bntseq_t *bns)
{
    if (bns == NULL) return NULL;
    meth_chrom_map_t *m = (meth_chrom_map_t *)calloc(1, sizeof(*m));
    if (m == NULL) return NULL;
    m->n_internal = bns->n_seqs;
    if (m->n_internal > 0) {
        m->out_tid      = (int *)    calloc((size_t)m->n_internal, sizeof(int));
        m->direction    = (char *)   calloc((size_t)m->n_internal, sizeof(char));
        m->output_names = (char **)  calloc((size_t)m->n_internal, sizeof(char *));
        m->output_lens  = (int64_t *)calloc((size_t)m->n_internal, sizeof(int64_t));
        if (!m->out_tid || !m->direction || !m->output_names || !m->output_lens) {
            meth_chrom_map_free(m);
            return NULL;
        }
    }

    for (int i = 0; i < m->n_internal; ++i) {
        const char *name = bns->anns[i].name;
        char prefix = name[0];
        const char *stripped = name;
        char dir = 0;
        if (prefix == 'f' || prefix == 'r') {
            stripped = name + 1;
            dir = prefix;
        }
        m->direction[i] = dir;

        int existing = -1;
        for (int j = 0; j < m->n_output; ++j) {
            if (strcmp(m->output_names[j], stripped) == 0) { existing = j; break; }
        }
        if (existing >= 0) {
            m->out_tid[i] = existing;
        } else {
            int idx = m->n_output++;
            m->output_names[idx] = meth_strdup_local(stripped);
            if (m->output_names[idx] == NULL) { meth_chrom_map_free(m); return NULL; }
            m->output_lens[idx] = bns->anns[i].len;
            m->out_tid[i] = idx;
        }
    }
    return m;
}

void meth_chrom_map_free(meth_chrom_map_t *m)
{
    if (m == NULL) return;
    free(m->out_tid);
    free(m->direction);
    if (m->output_names != NULL) {
        for (int i = 0; i < m->n_output; ++i) free(m->output_names[i]);
        free(m->output_names);
    }
    free(m->output_lens);
    free(m);
}

/* ------------------------------------------------------------------- */
/* BAM writer                                                           */
/* ------------------------------------------------------------------- */

/* Replace embedded tabs in a string with spaces (BAM CL tag fields can't
 * contain literal tabs — safer to just sanitize). Caller's buffer. */
static void sanitize_cl(char *s)
{
    for (; *s; ++s) if (*s == '\t') *s = ' ';
}

meth_bam_writer_t *meth_bam_writer_open(const char *path_or_dash,
                                        meth_chrom_map_t *cmap,
                                        const char *bwa_pg,
                                        const char *meth_pg_cl)
{
    if (path_or_dash == NULL || cmap == NULL) return NULL;
    meth_bam_writer_t *w = (meth_bam_writer_t *)calloc(1, sizeof(*w));
    if (w == NULL) return NULL;
    w->cmap = cmap;

    /* "wb0" = uncompressed BAM (BGZF level 0). samtools/htslib accept this. */
    w->fp = hts_open(path_or_dash, "wb0");
    if (w->fp == NULL) { free(w); return NULL; }

    w->hdr = sam_hdr_init();
    if (w->hdr == NULL) { hts_close(w->fp); free(w); return NULL; }

    /* @HD */
    if (sam_hdr_add_line(w->hdr, "HD", "VN", "1.6", "SO", "unsorted", NULL) < 0) goto fail;

    /* @SQ from consolidated chrom map */
    for (int i = 0; i < cmap->n_output; ++i) {
        char len_buf[32];
        snprintf(len_buf, sizeof(len_buf), "%lld", (long long)cmap->output_lens[i]);
        if (sam_hdr_add_line(w->hdr, "SQ",
                             "SN", cmap->output_names[i],
                             "LN", len_buf,
                             NULL) < 0) goto fail;
    }

    /* Original bwa-mem2 @PG */
    if (bwa_pg != NULL && bwa_pg[0] != '\0') {
        if (sam_hdr_add_lines(w->hdr, bwa_pg, 0) < 0) goto fail;
    }

    /* bwa-mem2-meth @PG */
    {
        char pg_buf[4096];
        char cl_copy[2048] = "bwa-mem2 mem --meth";
        if (meth_pg_cl != NULL && meth_pg_cl[0] != '\0') {
            strncpy(cl_copy, meth_pg_cl, sizeof(cl_copy) - 1);
            cl_copy[sizeof(cl_copy) - 1] = '\0';
            sanitize_cl(cl_copy);
        }
        snprintf(pg_buf, sizeof(pg_buf),
                 "@PG\tID:bwa-mem2-meth\tPN:bwa-mem2-meth\tVN:2.2.1-meth\tCL:%s\n",
                 cl_copy);
        if (sam_hdr_add_lines(w->hdr, pg_buf, 0) < 0) goto fail;
    }

    if (sam_hdr_write(w->fp, w->hdr) < 0) goto fail;
    return w;

fail:
    sam_hdr_destroy(w->hdr);
    hts_close(w->fp);
    free(w);
    return NULL;
}

int meth_bam_writer_write(meth_bam_writer_t *w, bam1_t *b)
{
    if (w == NULL || b == NULL) return -1;
    return sam_write1(w->fp, w->hdr, b);
}

int meth_bam_writer_close(meth_bam_writer_t *w)
{
    if (w == NULL) return 0;
    int rc = 0;
    if (w->hdr) { sam_hdr_destroy(w->hdr); w->hdr = NULL; }
    if (w->fp) {
        if (hts_close(w->fp) < 0) rc = -1;
        w->fp = NULL;
    }
    free(w);
    return rc;
}

/* ------------------------------------------------------------------- */
/* mem_aln_t -> bam1_t                                                  */
/* ------------------------------------------------------------------- */

/* bwa-mem2 stores CIGAR ops as 0=M, 1=I, 2=D, 3=S, 4=H (see add_cigar's
 * "MIDSH" lookup in bwamem.cpp). BAM's canonical encoding is 0=M, 1=I,
 * 2=D, 3=N, 4=S, 5=H. Remap on emit. */
static const uint32_t BAM_OP_FROM_MEM[5] = { 0, 1, 2, 4, 5 };

/* Length consumed on the reference by a CIGAR op (BAM-standard ops). */
static inline int bam_op_consumes_ref(int op)
{
    return op == 0 || op == 2 || op == 3 || op == 7 || op == 8;
}

/* Length consumed on the query by a CIGAR op (BAM-standard ops). */
/* unused helper — keep for future use */

static int64_t cigar_ref_len_bam(const uint32_t *cigar, int n)
{
    int64_t l = 0;
    for (int i = 0; i < n; ++i) {
        int op = cigar[i] & 0xf;
        int len = cigar[i] >> 4;
        if (bam_op_consumes_ref(op)) l += len;
    }
    return l;
}

/* Longest run of M/=/X (BAM-standard ops) for the chimera QC heuristic. */
static int cigar_longest_m_bam(const uint32_t *cigar, int n)
{
    int longest = 0;
    for (int i = 0; i < n; ++i) {
        int op = cigar[i] & 0xf;
        int len = cigar[i] >> 4;
        if ((op == 0 || op == 7 || op == 8) && len > longest) longest = len;
    }
    return longest;
}

int meth_mem_aln_to_bam(bam1_t *b,
                        const mem_opt_t *opt, const bntseq_t *bns,
                        const bseq1_t *s, int n_alns,
                        const mem_aln_t *list, int which,
                        const mem_aln_t *m_,
                        const meth_chrom_map_t *cmap)
{
    (void)n_alns; (void)bns;
    if (b == NULL || opt == NULL || s == NULL || list == NULL || cmap == NULL) return -1;

    /* Copies so we can mutate safely (parallels mem_aln2sam pattern). */
    mem_aln_t p = list[which];
    mem_aln_t m;
    const mem_aln_t *mp = NULL;
    if (m_ != NULL) { m = *m_; mp = &m; }

    /* Flag munging — identical to mem_aln2sam:1599-1608 */
    p.flag |= mp ? 0x1 : 0;
    p.flag |= p.rid < 0 ? 0x4 : 0;
    p.flag |= mp && mp->rid < 0 ? 0x8 : 0;
    if (p.rid < 0 && mp && mp->rid >= 0) {
        p.rid = mp->rid; p.pos = mp->pos; p.is_rev = mp->is_rev; p.n_cigar = 0;
    }
    if (mp && mp->rid < 0 && p.rid >= 0) {
        m.rid = p.rid; m.pos = p.pos; m.is_rev = p.is_rev; m.n_cigar = 0;
    }
    p.flag |= p.is_rev ? 0x10 : 0;
    p.flag |= mp && mp->is_rev ? 0x20 : 0;

    /* Supp/alt bit folded into high byte of flag (parallels mem_aln2sam:1614) */
    uint16_t flag16 = (uint16_t)((p.flag & 0xffff) | (p.flag & 0x10000 ? 0x100 : 0));

    /* Resolve output tids and direction */
    int32_t tid = -1, mtid = -1;
    char direction = 0;
    if (p.rid >= 0 && p.rid < cmap->n_internal) {
        tid       = cmap->out_tid[p.rid];
        direction = cmap->direction[p.rid];
    }
    if (mp && mp->rid >= 0 && mp->rid < cmap->n_internal) {
        mtid = cmap->out_tid[mp->rid];
    }

    /* Remap primary CIGAR: bwa-mem2 ops -> BAM ops, + soft->hard for supp */
    uint32_t *bam_cigar = NULL;
    size_t    bam_n_cigar = 0;
    if (p.n_cigar > 0) {
        bam_cigar = (uint32_t *)malloc((size_t)p.n_cigar * sizeof(uint32_t));
        if (bam_cigar == NULL) return -1;
        for (int i = 0; i < p.n_cigar; ++i) {
            int op  = p.cigar[i] & 0xf;
            int len = p.cigar[i] >> 4;
            if (!(opt->flag & MEM_F_SOFTCLIP) && !p.is_alt && (op == 3 || op == 4))
                op = which ? 4 : 3;              /* hard clip for supp */
            uint32_t bam_op = (op >= 0 && op < 5) ? BAM_OP_FROM_MEM[op] : 0;
            bam_cigar[i] = ((uint32_t)len << 4) | bam_op;
        }
        bam_n_cigar = (size_t)p.n_cigar;
    }

    /* TLEN — parallels mem_aln2sam:1641-1646 */
    hts_pos_t tlen = 0;
    if (mp && mp->rid >= 0 && p.rid == mp->rid && p.n_cigar > 0 && m.n_cigar > 0) {
        /* Build a BAM-remapped mate CIGAR temp for ref-len computation */
        uint32_t *mbc = (uint32_t *)malloc((size_t)m.n_cigar * sizeof(uint32_t));
        if (mbc == NULL) { free(bam_cigar); return -1; }
        for (int i = 0; i < m.n_cigar; ++i) {
            int op = m.cigar[i] & 0xf;
            int len = m.cigar[i] >> 4;
            uint32_t bop = (op >= 0 && op < 5) ? BAM_OP_FROM_MEM[op] : 0;
            mbc[i] = ((uint32_t)len << 4) | bop;
        }
        int64_t p_rlen = cigar_ref_len_bam(bam_cigar, (int)bam_n_cigar);
        int64_t m_rlen = cigar_ref_len_bam(mbc, m.n_cigar);
        free(mbc);
        int64_t p0 = p.pos + (p.is_rev ? p_rlen - 1 : 0);
        int64_t p1 = m.pos + (m.is_rev ? m_rlen - 1 : 0);
        tlen = -(p0 - p1 + (p0 > p1 ? 1 : p0 < p1 ? -1 : 0));
    }

    /* Compute SEQ/QUAL range with supp soft-clip trim */
    int qb = 0, qe = s->l_seq;
    if (p.n_cigar && which && !(opt->flag & MEM_F_SOFTCLIP) && !p.is_alt) {
        if (!p.is_rev) {
            int c0 = p.cigar[0] & 0xf;
            int cN = p.cigar[p.n_cigar-1] & 0xf;
            if (c0 == 3 || c0 == 4) qb += p.cigar[0] >> 4;
            if (cN == 3 || cN == 4) qe -= p.cigar[p.n_cigar-1] >> 4;
        } else {
            int c0 = p.cigar[0] & 0xf;
            int cN = p.cigar[p.n_cigar-1] & 0xf;
            if (c0 == 3 || c0 == 4) qe -= p.cigar[0] >> 4;
            if (cN == 3 || cN == 4) qb += p.cigar[p.n_cigar-1] >> 4;
        }
    }

    int emit_seq = !(p.flag & 0x100);
    size_t l_emit = 0;
    char *seq_text = NULL;
    char *qual_bin = NULL;
    if (emit_seq && qe > qb) {
        l_emit = (size_t)(qe - qb);
        seq_text = (char *)malloc(l_emit + 1);
        if (seq_text == NULL) { free(bam_cigar); return -1; }
        if (!p.is_rev) {
            for (size_t i = 0; i < l_emit; ++i) seq_text[i] = "ACGTN"[(int)s->seq[qb + (int)i]];
        } else {
            for (size_t i = 0; i < l_emit; ++i) seq_text[i] = "TGCAN"[(int)s->seq[qe - 1 - (int)i]];
        }
        seq_text[l_emit] = '\0';
        if (s->qual) {
            qual_bin = (char *)malloc(l_emit);
            if (qual_bin == NULL) { free(seq_text); free(bam_cigar); return -1; }
            if (!p.is_rev) {
                for (size_t i = 0; i < l_emit; ++i) qual_bin[i] = (char)((unsigned char)s->qual[qb + (int)i] - 33);
            } else {
                for (size_t i = 0; i < l_emit; ++i) qual_bin[i] = (char)((unsigned char)s->qual[qe - 1 - (int)i] - 33);
            }
        }
    }

    /* Chimera QC + set-as-failed (applied here so chimera propagation upstream
     * only needs to scan flags across a group). */
    uint8_t mapq = p.mapq;
    int mapped = !(flag16 & 0x4) && direction != 0;
    if (mapped) {
        if (opt->meth_set_as_failed != 0 && opt->meth_set_as_failed == direction) {
            flag16 |= 0x200;
        }
        if (!opt->meth_no_chim && bam_cigar && bam_n_cigar > 0 && s->l_seq > 0) {
            int lm = cigar_longest_m_bam(bam_cigar, (int)bam_n_cigar);
            if (100 * lm < 44 * s->l_seq) {
                flag16 |= 0x200;
                flag16 &= ~0x2;
                if (mapq > 1) mapq = 1;
            }
        }
    }

    /* Build the bam1_t. bam_set1 handles 4-bit packing, name storage, etc. */
    int ret = bam_set1(b,
                       strlen(s->name), s->name,
                       flag16,
                       tid,
                       (hts_pos_t)p.pos,
                       mapq,
                       bam_n_cigar, bam_cigar,
                       mtid,
                       mp ? (hts_pos_t)mp->pos : -1,
                       tlen,
                       l_emit, seq_text, qual_bin,
                       /* l_aux */ 0);

    free(bam_cigar);
    free(seq_text);
    free(qual_bin);
    if (ret < 0) return -1;

    /* Aux tags — roughly match mem_aln2sam emission order */
    if (p.n_cigar > 0) {
        int32_t nm = (int32_t)p.NM;
        bam_aux_append(b, "NM", 'i', sizeof(nm), (const uint8_t *)&nm);
        const char *md = (const char *)(p.cigar + p.n_cigar);
        bam_aux_append(b, "MD", 'Z', (int)strlen(md) + 1, (const uint8_t *)md);
    }
    if (mp && mp->n_cigar > 0) {
        char mc_buf[4096]; int mc_len = 0;
        for (int i = 0; i < mp->n_cigar && mc_len < (int)sizeof(mc_buf) - 16; ++i) {
            int op = mp->cigar[i] & 0xf;
            int len = mp->cigar[i] >> 4;
            if (!(opt->flag & MEM_F_SOFTCLIP) && !mp->is_alt && (op == 3 || op == 4))
                op = which ? 4 : 3;
            mc_len += snprintf(mc_buf + mc_len, sizeof(mc_buf) - mc_len,
                               "%d%c", len, "MIDSH"[op]);
        }
        if (mc_len > 0)
            bam_aux_append(b, "MC", 'Z', mc_len + 1, (const uint8_t *)mc_buf);
    }
    if (p.score >= 0) {
        int32_t as = (int32_t)p.score;
        bam_aux_append(b, "AS", 'i', sizeof(as), (const uint8_t *)&as);
    }
    if (p.sub >= 0) {
        int32_t xs = (int32_t)p.sub;
        bam_aux_append(b, "XS", 'i', sizeof(xs), (const uint8_t *)&xs);
    }
    if (bwa_rg_id[0]) {
        bam_aux_append(b, "RG", 'Z', (int)strlen(bwa_rg_id) + 1, (const uint8_t *)bwa_rg_id);
    }
    if (p.alt_sc > 0) {
        float pa_f = (float)((double)p.score / (double)p.alt_sc);
        bam_aux_append(b, "pa", 'f', sizeof(pa_f), (const uint8_t *)&pa_f);
    }
    if (p.XA != NULL) {
        bam_aux_append(b, "XA", 'Z', (int)strlen(p.XA) + 1, (const uint8_t *)p.XA);
    }
    /* YD:Z — meth strand hypothesis */
    if (mapped) {
        char yd[2] = { direction, '\0' };
        bam_aux_append(b, "YD", 'Z', 2, (const uint8_t *)yd);
    }

    return 0;
}

void meth_bam_group_propagate_qcfail(bam1_t **group, int n)
{
    if (group == NULL || n <= 0) return;
    int any_fail = 0;
    for (int i = 0; i < n; ++i) {
        if (group[i] != NULL && (group[i]->core.flag & 0x200)) { any_fail = 1; break; }
    }
    if (!any_fail) return;
    for (int i = 0; i < n; ++i) {
        if (group[i] == NULL) continue;
        group[i]->core.flag |= 0x200;
        group[i]->core.flag &= (uint16_t)~0x2;
    }
}
