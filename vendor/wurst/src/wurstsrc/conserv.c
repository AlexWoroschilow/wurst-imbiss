/*
 * conserv.c
 *
 *  Created on: Jan 4, 2012
 *      Author: ibondarenko
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "e_malloc.h"
#include "mprintf.h"

#include "prob_vec.h"
#include "pair_set.h"
#include "seq.h"
#include "conserv.h"
#include "matrix.h"
#include "multialign.h"
#include "prob_vec_i.h"
#include "scratch.h"
#include "read_seq_i.h"

/* ------------------- initveclist -----------------------------
 */
struct prob_vecs *
initveclist(size_t num) {
    struct prob_vecs * res = memset(E_MALLOC ( sizeof(struct prob_vecs)), 0,
            sizeof(struct prob_vecs));
    size_t siz = sizeof(struct prob_vec *) * num;
    res->pv = memset(E_MALLOC ( siz), 0, siz);
    res->maxnum = num;
    res->curnum = 0;

    return res;
}

/* ------------------- cleanveclist ----------------------------
 */
void
cleanveclist(struct prob_vecs *ps) {

    if (ps && (ps->curnum > 0)) {
        while (ps->curnum > 0) {
            ps->curnum--;
            if (ps->pv[ps->curnum]) {
                prob_vec_destroy(ps->pv[ps->curnum]);
                ps->pv[ps->curnum] = 0;
            }
        }
    }
}

/* ------------------- deleteveclist ---------------------------
 */
void
prob_vecs_destroy(struct prob_vecs *ps) {

    if (ps) {
        cleanveclist(ps);
        free(ps->pv);
        free(ps);
//        ps = 0;
    }
}

/* ------------------- initseqlist -----------------------------
 */
struct sequences_list *
initseqlist(size_t num) {

    struct sequences_list * res = memset(
            E_MALLOC ( sizeof(struct sequences_list)), 0,
            sizeof(struct sequences_list));
    size_t siz = sizeof(struct seq *) * num;
    res->ps = memset(E_MALLOC ( siz), 0, siz);
    res->maxnum = num;
    res->curnum = 0;

    return res;

}

/* ------------------- cleanseqlist  ---------------------------
 */
void
cleanseqlist(struct sequences_list *ps) {

    if (ps) {
        while (ps->curnum > 0) {
            ps->curnum--;
            if (ps->ps[ps->curnum]) {
                seq_destroy(ps->ps[ps->curnum]);
                ps->ps[ps->curnum] = 0;
            }
        }
    }
}


/* ------------------- deleteseqlist ---------------------------
 */
void
sequences_list_destroy(struct sequences_list *ps) {

    if (ps) {
        cleanseqlist(ps);
        free(ps->ps);
        free(ps);
//        ps = 0;
    }
}

/* ------------------- addvec ----------------------------------
 */
void
addvec(struct prob_vecs* pvs, struct prob_vec *vec) {

    if (pvs && (pvs->curnum < pvs->maxnum)) {
        pvs->pv[pvs->curnum] = prob_vec_copy(vec);
        pvs->curnum++;
    }
    mprintf("current num = %u \n", pvs->curnum);
}

/* ------------------- addseq ----------------------------------
 */
void
addseq(struct sequences_list* pvs, struct seq *seq) {

    if (pvs && (pvs->curnum < pvs->maxnum)) {
        pvs->ps[pvs->curnum] = seq_copy(seq);
        pvs->curnum++;
    }

}

/* ------------------- conserv_str -----------------------------
 */
char *
conserv_str(float * consv, size_t len) {

    size_t i, count = 0;
    scr_reset();
    for (i = 0; i < len; i++) {
        if (count > 0) {
            scr_printf(", ");
        }
        scr_printf("%1.3f", consv[i]);
        count++;
    }
    return (scr_printf("%c", '\0'));

}

/* ------------------- printconserv  ---------------------------
 */
char *
printconserv(float *consv, const size_t len, struct pair_set *pair_set,
        const size_t startpos, const size_t maxlen, const size_t strNum1)
{
    size_t i, count;
    extern const char *null_point;
    const char *this_sub = "printconserv";
    if (!pair_set) {
        err_printf(this_sub, null_point);
        return NULL;
    }

    if (startpos >= len) {
        err_printf(this_sub, "startpos is to big.");
        return NULL;
    }

    count = 0;
    scr_reset();
    for (i = 0; i < len; i++) {
        if (startpos <= pair_set->indices[i][strNum1]) {
            if (pair_set->indices[i][strNum1] != GAP_INDEX) {
                if (maxlen > count) {
                    if (count > 0) {
                        scr_printf(", ");
                    }
                    scr_printf("%1.3f", consv[i]);
                    count++;
                } else {
                    break;
                }
            }
        }
    }

    return (scr_printf("%c", '\0'));
}

/* ------------------- get_conserve  ---------------------------
 */
float *
get_conserv(struct prob_vec* pv) {

    char *this_sub = "get_conserv";
    float *res;
    size_t j, k, siz;
    if (prob_vec_expand(pv) == EXIT_FAILURE) {
        err_printf(this_sub, "fail on vec \n");
        return NULL;
    }

    siz = sizeof(float) * pv->n_pvec + 1;
    res = memset(E_MALLOC (siz), 0, siz);

    for (j = 0; j < pv->n_pvec; j++) { /* align len */
        for (k = 0; k < pv->n_class - 1; k++) {
            res[j] += (pv->mship[j][k] * pv->mship[j][k]);
        }
        res[j] = sqrt(res[j]);
    }

    return res;
}

/* ------------------- get_seq_conserv -------------------------
 */
float *
get_seq_conserv(struct pair_set *pair_set, struct sequences_list* ps) {
    float **mtx;
    float *res;
    size_t i, m, siz, t;
    char inx;
    const size_t letter_NUM = 27;
    const float Residue_NUM = 21;

    extern const char *null_point;
    const char *this_sub = "get_seq_sonserv";
    if (!pair_set) {
        err_printf(this_sub, null_point);
        return NULL;
    }

    for (i = 0; i < ps->curnum; i++) {
        seq_thomas2std(ps->ps[i]);
    }

    mtx = f_matrix(letter_NUM, pair_set->n);
    memset(mtx[0], 0, letter_NUM * pair_set->n * sizeof(mtx[0][0]));
    for (i = 0; i < pair_set->n; i++) {
        for (m = 0; m < pair_set->m; m++) {
            if (pair_set->indices[i][m] != GAP_INDEX) {
                t = pair_set->indices[i][m];
                inx = ps->ps[m]->seq[t];
                mtx[inx - '`'][i]++;
            } else {
                mtx[0][i]++;
            }
        }
    }

    /* normalisierung */
    for (i = 0; i < pair_set->n; i++) {
        for (m = 0; m < letter_NUM; m++) {
            mtx[m][i] /= pair_set->m;
        }
    }

    /* conservierung */
    siz = sizeof(res[0]) * pair_set->n + 1;
    res = memset(E_MALLOC (siz), 0, siz);

    for (i = 0; i < pair_set->n; i++) {
        for (m = 0; m < letter_NUM; m++) {
            if (mtx[m][i]) {
                res[i] += mtx[m][i] * log(mtx[m][i]);
            }
        }
        res[i] = 1 + (res[i] / log(Residue_NUM));
    }

    kill_f_matrix(mtx);
    return (res);

}

/* ------------------- getconservvec ---------------------------
 */
float *
getconservvec(struct pair_set *pair_set, struct prob_vecs* pvs) {

    struct prob_vec *pvec;
    size_t i, j, jj, m, k;
    long tt;
    float *res, norm;
    size_t *counter;

    extern const char *null_point;
    const char *this_sub = "getconservvec";
    if (!pair_set) {
        err_printf(this_sub, null_point);
        return NULL;
    }

    for (i = 0; i < pair_set->m - 1; i++) {
        if (pvs->pv[i]->frag_len != pvs->pv[i + 1]->frag_len) {
            err_printf(this_sub, "can not add probability vectors with ");
            err_printf(this_sub, "different fragment lengths!\n");
            return NULL;
        }
        if (pvs->pv[i]->n_class != pvs->pv[i + 1]->n_class) {
            err_printf(this_sub, "can not add probability vectors with ");
            err_printf(this_sub, "different number of classes!\n");
            return NULL;
        }

        prob_vec_unit_vec(pvs->pv[i]);
        if (prob_vec_expand(pvs->pv[i]) == EXIT_FAILURE) {
            err_printf(this_sub, "fail on vec %d \n", i);
            return NULL;
        }
    }

    prob_vec_unit_vec(pvs->pv[pair_set->m - 1]);
    if (prob_vec_expand(pvs->pv[pair_set->m - 1]) == EXIT_FAILURE) {
        err_printf(this_sub, "fail on vec %d \n", pair_set->m - 1);
        return NULL;
    }

    pvec = new_pvec(pvs->pv[0]->frag_len, pair_set->n + pvs->pv[0]->frag_len
            - 1, pair_set->n, pvs->pv[0]->n_class);
    if (!pvec) {
        err_printf(this_sub, "new_pvec returned NULL. \n");
        return NULL;
    }

    counter = memset(E_MALLOC (pair_set->n *sizeof(size_t)), 0, pair_set->n
            * sizeof(size_t));
    /*------------------------------------- TODO ----------------------------------------------
     * check that mapping between residues and pvecs is preserved during the averaging.
     * if the averaging starts to slide the alignment even slightly, all sorts of nasty things
     * happen. -especially if chains have different lengths or the alignment contains overhangs
     * at the ends.\
    * ----------------------------------------------------------------------------------------*/

    memset(pvec->mship[0], 0, pvec->n_pvec * pvec->n_class
            * sizeof(pvec->mship[0][0]));
    for (i = 1; i < pvec->n_pvec; i++)
        pvec->mship[i] = pvec->mship[i - 1] + pvec->n_class;

    for (i = 0; i < pvec->n_pvec - pvec->frag_len + 1; i++) { /* align len  - fragment length */
        for (jj = 0; jj < pvec->frag_len; jj++) { /* throw a fragment */
            for (m = 0; m < pair_set->m; m++) { /* alle seqs im algnm */
                if (pair_set->indices[i + jj][m] != GAP_INDEX) {
                    tt
                            = (pair_set->indices[i][m] == GAP_INDEX) ? (pair_set->indices[i
                                    + jj][m] - jj)
                                    : pair_set->indices[i][m];
                    if (tt >= 0 && pvs->pv[m]->n_pvec > tt) {
                        for (k = 0; k < pvs->pv[0]->n_class; k++) {
                            pvec->mship[i + jj][k] += pvs->pv[m]->mship[tt][k];
                        }
                        counter[i + jj]++;
                    }
                }

            }

        }
    }

    for (i = 0; i < pvec->n_pvec; i++) { /* align len  - fragment length */
        for (k = 0; k < pvs->pv[0]->n_class; k++) {
            if (counter[i]) {
                pvec->mship[i][k] /= counter[i];
            }
        }
    }
    res = get_conserv(pvec);

    norm = pvec->frag_len * pair_set->m;
    for (i = 0; i < pvec->n_pvec; i++) { /* normalisierung for GAPs*/
        res[i] *= counter[i] / norm;
    }

    pvec->cmpct_n = NULL;
    pvec->cmpct_prob = NULL;
    pvec->cmpct_ndx = NULL;
    pvec->norm_type = PVEC_TRUE_PROB;
    prob_vec_destroy(pvec);
    free(counter);

    return (res);
}

