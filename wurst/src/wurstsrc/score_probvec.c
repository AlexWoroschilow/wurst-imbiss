/*
 * 14 June 2005
 * This will take a pair of class membership probability vectors
 * and fill out a score matrix.
 * An optimising note. I just checked the assembler output and
 * with -O2, the compiler is recognising the inner loops and
 * seems to be using the streaming instructions. This is despite
 * the indirect indexing used in the compact forms.
 *
 * $Id: score_probvec.c,v 1.9 2014/10/21 15:08:28 torda Exp $
 */

#include <stdio.h>
#include <stdlib.h>

#include "e_malloc.h"
#include "mprintf.h"
#include "prob_vec.h"
#include "prob_vec_i.h"
#include "score_mat.h"
#include "score_mat_i.h"
#include "score_probvec.h"

/* ---------------- enum_types  -------------------------------
 */
enum compact_state {
    BOTH_EXPANDED,
    BOTH_COMPACT,
    ONE_COMPACT,
    TWO_COMPACT
};

/* ---------------- sanity_checks -----------------------------
 * Perform any checks on the arguments.
 * Of the errors we check for, only one seems to be a problem during
 * searches, so we put a cap on it.
 */
static int
sanity_checks (const struct score_mat *score_mat,
            const struct prob_vec *p_v1, const struct prob_vec *p_v2)
{
    const size_t frag_len1 = p_v1->frag_len;
    const size_t frag_len2 = p_v2->frag_len;
    const size_t n1 = p_v1->prot_len;      /* The length of the */
    const size_t n2 = p_v2->prot_len;      /* original protein */
    const unsigned short MAX_MISMATCH_ERR = 10;

    const char *this_sub = "score_pvec";       /* This is a lie, not a typo */
    const char *frag_mismatch =
        "Fragment length in prob vector 1 and 2 different. %d != %d\n";
    const char *class_mismatch =
        "Prob vectors have different numbers of classes. %d != %d\n";
    const char *max_err_print = "Max num error messages. Not printing them any more.\n";
    extern const char *mismatch;
    if (p_v1->norm_type == PVEC_CRAP || p_v1->norm_type == PVEC_CRAP) {
        err_printf (this_sub, "old code warning.  I must die.\n");
        exit (EXIT_FAILURE);
    }
    if (( score_mat->n_rows != n1 + 2) || ( score_mat->n_cols != n2 + 2)) {
        static unsigned short n_printed = 0;
        if (n_printed < MAX_MISMATCH_ERR) {
            err_printf (this_sub, mismatch,
                        score_mat->n_rows-2, score_mat->n_cols-2, n1, n2);
            n_printed++;
            if (n_printed == MAX_MISMATCH_ERR)
                err_printf (this_sub, max_err_print, MAX_MISMATCH_ERR);
        }
        return EXIT_FAILURE;
    }
    if ( p_v1->n_class != p_v2->n_class) {
        err_printf (this_sub, class_mismatch, p_v1->n_class, p_v2->n_class);
        return EXIT_FAILURE;
    }
    if (frag_len1 != frag_len2) {
        err_printf (this_sub, frag_mismatch, (int) frag_len1, frag_len2);
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}

/* ---------------- one_cmpct_fill_matrix ---------------------
 * Of the two probability vectors, the first argument is in
 * compact form.
 */
static void
one_cmpct_fill_matrix (struct score_mat *score_mat,
                       const struct prob_vec *p_v1,
                       const struct prob_vec *p_v2)
{
    float **mship2 = p_v2->mship;
    float **scores = score_mat->mat;
    float *valtmp;
    unsigned *ndxtmp;
    const unsigned short *pi    = p_v1->cmpct_n;
    const unsigned short *p_ndx = p_v1->cmpct_ndx;
    const float *p_prob         = p_v1->cmpct_prob;
    const size_t frag_len = p_v1->frag_len;
    size_t i, j, k;
    unsigned short a, b;
    
    valtmp = E_MALLOC (p_v1->n_class * sizeof (valtmp[0]));
    ndxtmp = E_MALLOC (p_v1->n_class * sizeof (ndxtmp[0]));
    for (i = 0; i < p_v1->n_pvec; i++) {     /* Loop over sites in prot1 */
        const unsigned short n_used = pi[i];
        for (a = 0; a < n_used; a++) {
            ndxtmp[a] = *p_ndx++;             /* Store the class index and */
            valtmp[a] = *p_prob++;            /* corresponding probability */
        }
        for (j = 0; j < p_v2->n_pvec; j++) {      /* Over sites in prot2 */
            float score = 0.0;
            for (b = 0; b < n_used; b++)
                score += valtmp[b] * mship2[j][ndxtmp[b]];
            for (k = 0; k < frag_len; k++)
                scores [i + k + 1][j + k + 1] += score;
        }
    }
    free (ndxtmp);
    free (valtmp);
}

/* ---------------- two_cmpct_fill_matrix ---------------------
 * Of the two probability vectors, the second is in compact
 * form. What I do not like is that the inner loop over scores[][]
 * is now indexed the wrong way [j][i].
 */
static void
two_cmpct_fill_matrix (struct score_mat *score_mat,
                       const struct prob_vec *p_v1,
                       const struct prob_vec *p_v2)
{
    float **mship1 = p_v1->mship;
    float **scores = score_mat->mat;
    float *valtmp;
    unsigned *ndxtmp;
    const unsigned short *pi    = p_v2->cmpct_n;
    const unsigned short *p_ndx = p_v2->cmpct_ndx;
    const float *p_prob         = p_v2->cmpct_prob;
    const size_t frag_len = p_v2->frag_len;
    size_t i, j, k;
    unsigned short a, b;

    valtmp = E_MALLOC (p_v2->n_class * sizeof (valtmp[0]));
    ndxtmp = E_MALLOC (p_v2->n_class * sizeof (ndxtmp[0]));
    for (i = 0; i < p_v2->n_pvec; i++) {     /* Loop over sites in prot2 */
        const unsigned short n_used = pi[i];
        for (a = 0; a < n_used; a++) {
            ndxtmp[a] = *p_ndx++;             /* Store the class index and */
            valtmp[a] = *p_prob++;            /* corresponding probability */
        }
        for (j = 0; j < p_v1->n_pvec; j++) {      /* Over sites in prot1 */
            float score = 0.0;
            for (b = 0; b < n_used; b++)
                score += valtmp[b] * mship1[j][ndxtmp[b]];
            for (k = 0; k < frag_len; k++)
                scores [j + k + 1][i + k + 1] += score;
        }
    }
    free (ndxtmp);
    free (valtmp);
}

/* ---------------- exp_fill_matrix ---------------------------
 * For the case of an both vector lists being expanded, fill
 * the matrix.
 */
static void
exp_fill_matrix (struct score_mat *score_mat,
                 const struct prob_vec *p_v1, const struct prob_vec *p_v2)
{
    float **mship1 = p_v1->mship;
    float **mship2 = p_v2->mship;
    float **scores = score_mat->mat; 

    size_t i;
    const size_t frag_len  = p_v1->frag_len;
    const size_t n_class = p_v1->n_class;

    for (i = 0; i < p_v1->n_pvec; i++) {
        size_t j;
        for (j = 0; j < p_v2->n_pvec; j++) {
            size_t k, m;
            float score = 0;
            for (m = 0; m < n_class; m++)        /* Dot product of class */
                score += mship1[i][m] * mship2[j][m];   /* probabilities */
            for (k = 0; k < frag_len; k++)          /* Put elements into */
                scores [i + k + 1] [j + k + 1] += score; /* score matrix */
        }
    }
}

/* ---------------- cmpct_fill_matrix -------------------------
 * Both of our vectors are in compact form. The easiest way to
 * treat this is to expand the second vector and call the function
 * which already handles this case.
 */
static void
cmpct_fill_matrix(struct score_mat *score_mat,
            const struct prob_vec *p_v1, struct prob_vec *p_v2)
{
    prob_vec_expand (p_v2);
    one_cmpct_fill_matrix (score_mat, p_v1, p_v2);
}

/* ---------------- score_pvec  -------------------------------
 * This is the core matrix scoring routine used by many functions.
 * It is complicated because the probability vectors could be passed
 * in as expanded or compacted.
 */
int
score_pvec (struct score_mat *score_mat,
            struct prob_vec *p_v1, struct prob_vec *p_v2)
{
    float **scores;
    size_t i;
    enum compact_state c_state;
    const size_t n_rows = p_v1->prot_len + 2;
    const size_t n_cols = p_v2->prot_len + 2;
    const char *this_sub = "score_pvec";
    extern const char *prog_bug;
    if (sanity_checks (score_mat, p_v1, p_v2) == EXIT_FAILURE)
        return EXIT_FAILURE;
    if (p_v1->mship && p_v2->mship)
        c_state = BOTH_EXPANDED;
    else if ( p_v1->mship == NULL && p_v2->mship == NULL)
        c_state = BOTH_COMPACT;
    else if ( p_v1->mship)
        c_state = TWO_COMPACT;
    else
        c_state = ONE_COMPACT;

    scores = score_mat->mat;
    for (i = 0; i < n_cols ; i++)
        scores[0][i] = 0;                   /* First row */
    for (i = 0; i < n_cols ; i++)
        scores[n_rows - 1][i] = 0;          /* Last row */
    for (i = 0; i < n_rows ; i++)
        scores [i][0] = 0;                  /* First column */
    for (i = 0; i < n_rows ; i++)
        scores [i][n_cols - 1 ] = 0;        /* Last column */

    prob_vec_unit_vec (p_v1);      /* We have to make sure both vectors */
    prob_vec_unit_vec (p_v2);      /* are properly normalised */

    switch ( c_state) {
    case BOTH_EXPANDED: exp_fill_matrix (score_mat, p_v1, p_v2);        break;
    case BOTH_COMPACT:  cmpct_fill_matrix (score_mat, p_v1, p_v2);      break;
    case ONE_COMPACT:   one_cmpct_fill_matrix (score_mat, p_v1, p_v2);  break;
    case TWO_COMPACT:   two_cmpct_fill_matrix (score_mat, p_v1, p_v2);  break;
    default: err_printf (this_sub, prog_bug, __FILE__, __LINE__);
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}



