#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <signal.h>

#include "align_i.h"
#include "e_malloc.h"
#include "matrix.h"
#include "mprintf.h"
#include "multialign.h"
#include "prob_vec.h"
#include "prob_vec_i.h"
#include "score_mat.h"
#include "score_mat_i.h"
#include "score_probvec.h"
#include "pair_set.h"
#include "pair_set_i.h"

/* --------------------------- pvec_avg ---------------------------------
 * compute average pvecs. this shouldn't be here and will be moved to pvec.c
 * later. I'm only keeping it here for convinience while debugging the bastard.
 * Some more of Akira's code will probably be imported later.
 */

struct prob_vec *
pvec_avg(struct prob_vec *p_vec1, struct prob_vec *p_vec2,
        struct pair_set *p_set, float weight) {
    struct prob_vec *pvec;
    int **pairs = p_set->indices;
    size_t frag_len = p_vec1->frag_len;
    size_t n_class = p_vec1->n_class;
    size_t align_len = p_set->n;
    size_t i, j, k;
    /*float weight = cur_step / (cur_step + 1);*/
    const char *this_sub = "prob_vec_add2";

    if (frag_len != p_vec2->frag_len) {
        err_printf(this_sub, "can not add probability vectors with ");
        err_printf(this_sub, "different fragment lengths!\n");
    }
    if (n_class != p_vec2->n_class) {
        err_printf(this_sub, "can not add probability vectors with ");
        err_printf(this_sub, "different number of classes!\n");
    }

    prob_vec_unit_vec(p_vec1);
    prob_vec_unit_vec(p_vec2);

    pvec = new_pvec(frag_len, align_len + frag_len - 1, align_len, n_class);
    memset(pvec->mship[0], 0, align_len * n_class * sizeof(float));
    j = 0;

    if (prob_vec_expand(p_vec1) == EXIT_FAILURE) {
        err_printf(this_sub, "fail on vec 1\n");
        return NULL;
    }
    if (prob_vec_expand(p_vec2) == EXIT_FAILURE) {
        err_printf(this_sub, "fail on vec 2\n");
        return NULL;
    }

    /*------------------------------------- TODO ----------------------------------------------
     * check that mapping between residues and pvecs is preserved during the averaging.
     * if the averaging starts to slide the alignment even slightly, all sorts of nasty things
     * happen. -especially if chains have different lengths or the alignment contains overhangs
     * at the ends.\
 * ----------------------------------------------------------------------------------------*/

    for (i = 0; i < align_len; i++) {
        if (pairs[i][0] != GAP_INDEX && pairs[i][0] < p_vec1->n_pvec) { /* no gap in struct 1 */
            if (pairs[i][1] == GAP_INDEX || pairs[i][1] >= p_vec1->n_pvec) { /* gap in struct2     */
                for (k = 0; k < pvec->n_class; k++) {
                    pvec->mship[j][k] = p_vec1->mship[pairs[i][0]][k];
                }
            }
        } else if (pairs[i][0] == GAP_INDEX || pairs[i][0] >= p_vec1->n_pvec) { /* gap in struct 1 */
            if (pairs[i][1] != GAP_INDEX && pairs[i][1] < p_vec2->n_pvec) {
                for (k = 0; k < pvec->n_class; k++) {
                    pvec->mship[j][k] = p_vec2->mship[pairs[i][1]][k];
                }
            }
        } else {
            for (k = 0; k < pvec->n_class; k++) {
                pvec->mship[j][k] = (p_vec1->mship[pairs[i][0]][k]) * weight
                        + (p_vec2->mship[pairs[i][1]][k]) * (1 - weight);
            }
        }
        j++;
    }

    pvec->cmpct_n = NULL;
    pvec->cmpct_prob = NULL;
    pvec->cmpct_ndx = NULL;

    pvec->norm_type = PVEC_TRUE_PROB;
    prob_vec_unit_vec(pvec);

    return (pvec);
}

/* ---------------- merge_align ------------------------
 * combine two alignments into a multiple alignment
 */

struct pair_set *
merge_align(struct pair_set *align1, struct pair_set *align2,
        struct pair_set *alignment) {
    struct pair_set *result;
    size_t j, col, idx1, idx2;
    j = col = idx1 = idx2 = 0;

    result = E_MALLOC (sizeof (struct pair_set));
    result->m = align1->m + align2->m;
    result->n = alignment->n;
    result->indices = i_matrix(result->n, result->m);
    /* now we copy the values from the old pairsets to the new one and insert gaps where the consensus
     * alignment tells us to*/

    for (col = 0; col < alignment->n; col++) {
        if (alignment->indices[col][0] != GAP_INDEX && idx1 < align1->n) {
            for (j = 0; j < align1->m; j++) {
                result->indices[col][j] = align1->indices[idx1][j];
            }
            idx1++;
        } else {
            for (j = 0; j < align1->m; j++) {
                result->indices[col][j] = GAP_INDEX;
            }
        }
        if (alignment->indices[col][1] != GAP_INDEX && idx2 < align2->n) {
            for (j = 0; j < align2->m; j++) {
                result->indices[col][j + align1->m] = align2->indices[idx2][j];
            }
            idx2++;
        } else {
            for (j = 0; j < align2->m; j++) {
                result->indices[col][j + align1->m] = GAP_INDEX;
            }
        }
    }
    return (result);

}

/* ---------------- merge_align ------------------------
 * combine two alignments into a multiple alignment
 */

struct pair_set *
merge_localigns(struct pair_set *align1, struct pair_set *align2, int c1,
        int c2) {
    struct pair_set *result;
    size_t j, col, idx1, idx2;
    j = col = idx1 = idx2 = 0;

    result = E_MALLOC (sizeof (struct pair_set));
    result->m = align1->m + align2->m;

    /* now we copy the values from the old pairsets to the new one and insert gaps where the consensus
     * alignment tells us to*/

    while ((idx1 < align1->n) && (idx2 < align2->n)) {
        if (align1->indices[idx1][c1] == GAP_INDEX && align2->indices[idx2][c2]
                == GAP_INDEX) {
            idx1++;
            idx2++;
        } else if (align1->indices[idx1][c1] == GAP_INDEX) {
            idx1++;
        } else if (align2->indices[idx2][c2] == GAP_INDEX) {
            idx2++;
        } else if (align1->indices[idx1][c1] < align2->indices[idx2][c2]) {
            idx1++;
        } else if (align1->indices[idx1][c1] > align2->indices[idx2][c2]) {
            idx2++;
        } else {
            idx1++;
            idx2++;
        }
        col++;
    }
    while (idx1 < align1->n) {
        idx1++;
        col++;
    }
    while (idx2 < align2->n) {
        idx2++;
        col++;
    }
    result->n = col;
    result->indices = i_matrix(result->n, result->m);
    col = idx1 = idx2 = 0;
    while ((idx1 < align1->n) && (idx2 < align2->n)) {
        if (align1->indices[idx1][c1] == GAP_INDEX && align2->indices[idx2][c2]
                == GAP_INDEX) {
            for (j = 0; j < align1->m; j++) {
                result->indices[col][j] = align1->indices[idx1][j];
            }
            for (j = 0; j < align2->m; j++) {
                result->indices[col][j + align1->m] = align2->indices[idx2][j];
            }
            idx1++;
            idx2++;
        } else if (align1->indices[idx1][c1] == GAP_INDEX) {
            for (j = 0; j < align1->m; j++) {
                result->indices[col][j] = align1->indices[idx1][j];
            }
            for (j = 0; j < align2->m; j++) {
                result->indices[col][j + align1->m] = GAP_INDEX;
            }
            idx1++;
        } else if (align2->indices[idx2][c2] == GAP_INDEX) {
            for (j = 0; j < align1->m; j++) {
                result->indices[col][j] = GAP_INDEX;
            }
            for (j = 0; j < align2->m; j++) {
                result->indices[col][j + align1->m] = align2->indices[idx2][j];
            }
            idx2++;
        } else if (align1->indices[idx1][c1] < align2->indices[idx2][c2]) {
            for (j = 0; j < align1->m; j++) {
                result->indices[col][j] = align1->indices[idx1][j];
            }
            for (j = 0; j < align2->m; j++) {
                result->indices[col][j + align1->m] = GAP_INDEX;
            }
            idx1++;
        } else if (align1->indices[idx1][c1] > align2->indices[idx2][c2]) {
            for (j = 0; j < align1->m; j++) {
                result->indices[col][j] = GAP_INDEX;
            }
            for (j = 0; j < align2->m; j++) {
                result->indices[col][j + align1->m] = align2->indices[idx2][j];
            }
            idx2++;
        } else {
            for (j = 0; j < align1->m; j++) {
                result->indices[col][j] = align1->indices[idx1][j];
            }
            for (j = 0; j < align2->m; j++) {
                result->indices[col][j + align1->m] = align2->indices[idx2][j];
            }
            idx1++;
            idx2++;
        }
        col++;
    }
    while (idx1 < align1->n) {
        for (j = 0; j < align1->m; j++) {
            result->indices[col][j] = align1->indices[idx1][j];
        }
        for (j = 0; j < align2->m; j++) {
            result->indices[col][j + align1->m] = GAP_INDEX;
        }
        idx1++;
        col++;
    }
    while (idx2 < align2->n) {
        for (j = 0; j < align1->m; j++) {
            result->indices[col][j] = GAP_INDEX;
        }
        for (j = 0; j < align2->m; j++) {
            result->indices[col][j + align1->m] = align2->indices[idx2][j];
        }
        idx2++;
        col++;
    }

    return (result);

}

/* ---------------- merge_alignments ------------------------
 * combine two alignments into a multiple alignment
 */

struct pair_set *
merge_alignments(struct pair_set *align1, struct pair_set *align2,
        struct pair_set *alignment, int c1, int c2) {
    struct pair_set *result;
    size_t i, j, idx, idx1, idx2;
    int count;
    count = 0;
    /* now we copy the values from the old pairsets to the new one and insert gaps where the consensus
     * alignment tells us to*/
    i = j = idx = idx1 = idx2 = 0;
    /* checking size of merged alignment */
    while ((idx < alignment->n) && (idx1 < align1->n) && (idx2 < align2->n)) {
        if ((alignment->indices[idx][0] != GAP_INDEX)
                && (alignment->indices[idx][1] != GAP_INDEX)) {
            if ((align1->indices[idx1][c1] != GAP_INDEX)
                    && (align2->indices[idx2][c2] != GAP_INDEX)) {
                idx1++;
                idx2++;
                idx++;
            } else if ((align1->indices[idx1][c1] == GAP_INDEX)
                    && (align2->indices[idx2][c2] == GAP_INDEX)) {
                idx1++;
                idx2++;
            } else if (align1->indices[idx1][c1] == GAP_INDEX) {
                idx1++;
            } else {
                idx2++;
            }
        } else if (alignment->indices[idx][0] == GAP_INDEX) {
            if ((align1->indices[idx1][c1] != GAP_INDEX)
                    && (align2->indices[idx2][c2] != GAP_INDEX)) {
                idx2++;
                idx++;
            } else if ((align1->indices[idx1][c1] == GAP_INDEX)
                    && (align2->indices[idx2][c2] != GAP_INDEX)) {
                idx1++;
                idx2++;
                idx++;

            } else if (align2->indices[idx2][c2] == GAP_INDEX) {
                idx2++;
            } else { /* ??? */
                idx1++;
            }
        } else if (alignment->indices[idx][1] == GAP_INDEX) {
            if ((align1->indices[idx1][c1] != GAP_INDEX)
                    && (align2->indices[idx2][c2] != GAP_INDEX)) {
                idx1++;
                idx++;
            } else if ((align1->indices[idx1][c1] != GAP_INDEX)
                    && (align2->indices[idx2][c2] == GAP_INDEX)) {
                idx1++;
                idx2++;
                idx++;
            } else if (align1->indices[idx1][c1] == GAP_INDEX) {
                idx1++;
            } else { /* ??? */
                idx2++;
            }
        } else {
            mprintf("missed a spot...\n");
        }
        count++;
    }
    result = E_MALLOC (sizeof (struct pair_set));
    result->m = (align1->m) + (align2->m);
    result->n = count;
    result->indices = i_matrix(result->n, result->m);
    i = j = idx = idx1 = idx2 = 0;
    for (i = 0; i < result->n; i++) {
        for (j = 0; j < result->m; j++) {
            result->indices[j][i] = GAP_INDEX;
        }
    }
    i = j = idx = idx1 = idx2 = 0;
    /* run through all columns of the alignment */
    for (i = 0; i < result->n; i++) {
        if ((alignment->indices[idx][0] != GAP_INDEX)
                && (alignment->indices[idx][1] != GAP_INDEX)) {
            if ((align1->indices[idx1][c1] != GAP_INDEX)
                    && (align2->indices[idx2][c2] != GAP_INDEX)) {
                for (j = 0; j < align1->m; j++) {
                    result->indices[i][j] = align1->indices[idx1][j];
                }
                idx1++;
                for (j = 0; j < align2->m; j++) {
                    result->indices[i][j + align1->m]
                            = align2->indices[idx2][j];
                }
                idx2++;
                idx++;
            } else if ((align1->indices[idx1][c1] == GAP_INDEX)
                    && (align2->indices[idx2][c2] == GAP_INDEX)) {
                for (j = 0; j < align1->m; j++) {
                    result->indices[i][j] = align1->indices[idx1][j];
                }
                idx1++;
                for (j = 0; j < align2->m; j++) {
                    result->indices[i][j + align1->m]
                            = align2->indices[idx2][j];
                }
                idx2++;
            } else if (align1->indices[idx1][c1] == GAP_INDEX) {
                for (j = 0; j < align1->m; j++) {
                    result->indices[i][j] = align1->indices[idx1][j];
                }
                for (j = align1->m; j < result->m; j++) {
                    result->indices[i][j] = GAP_INDEX;
                }
                idx1++;
            } else {
                for (j = 0; j < align1->m; j++) {
                    result->indices[i][j] = GAP_INDEX;
                }
                for (j = align1->m; j < result->m; j++) {
                    result->indices[i][j] = align2->indices[idx2][j];
                }
                idx2++;
            }
        } else if (alignment->indices[idx][0] == GAP_INDEX) {
            if ((align1->indices[idx1][c1] != GAP_INDEX)
                    && (align2->indices[idx2][c2] != GAP_INDEX)) {
                for (j = 0; j < align1->m; j++) {
                    result->indices[i][j] = GAP_INDEX;
                }
                for (j = align1->m; j < result->m; j++) {
                    result->indices[i][j] = align2->indices[idx2][j];
                }
                idx2++;
                idx++;
            } else if ((align1->indices[idx1][c1] == GAP_INDEX)
                    && (align2->indices[idx2][c2] != GAP_INDEX)) {
                for (j = 0; j < align1->m; j++) {
                    result->indices[i][j] = align1->indices[idx1][j];
                }
                idx1++;
                for (j = 0; j < align2->m; j++) {
                    result->indices[i][j + align1->m]
                            = align2->indices[idx2][j];
                }
                idx2++;
                idx++;

            } else if (align2->indices[idx2][c2] == GAP_INDEX) {
                for (j = 0; j < align1->m; j++) {
                    result->indices[i][j] = GAP_INDEX;
                }
                for (j = align1->m; j < result->m; j++) {
                    result->indices[i][j] = align2->indices[idx2][j];
                }
                idx2++;
            } else { /* ??? */
                for (j = 0; j < align1->m; j++) {
                    result->indices[i][j] = GAP_INDEX;
                }
                idx1++;
                for (j = 0; j < align2->m; j++) {
                    result->indices[i][j + align1->m]
                            = align2->indices[idx2][j];
                }
            }
        } else if (alignment->indices[idx][1] == GAP_INDEX) {
            if ((align1->indices[idx1][c1] != GAP_INDEX)
                    && (align2->indices[idx2][c2] != GAP_INDEX)) {
                for (j = 0; j < align1->m; j++) {
                    result->indices[i][j] = align1->indices[idx1][j];
                }
                for (j = align1->m; j < result->m; j++) {
                    result->indices[i][j] = GAP_INDEX;
                }
                idx1++;
                idx++;
            } else if ((align1->indices[idx1][c1] != GAP_INDEX)
                    && (align2->indices[idx2][c2] == GAP_INDEX)) {
                for (j = 0; j < align1->m; j++) {
                    result->indices[i][j] = align1->indices[idx1][j];
                }
                idx1++;
                for (j = 0; j < align2->m; j++) {
                    result->indices[i][j + align1->m]
                            = align2->indices[idx2][j];
                }
                idx2++;
                idx++;
            } else if (align1->indices[idx1][c1] == GAP_INDEX) {
                for (j = 0; j < align1->m; j++) {
                    result->indices[i][j] = align1->indices[idx1][j];
                }
                for (j = align1->m; j < result->m; j++) {
                    result->indices[i][j] = GAP_INDEX;
                }
                idx1++;
            } else { /* ??? */

                for (j = 0; j < align1->m; j++) {
                    result->indices[i][j] = align1->indices[idx1][j];
                }
                for (j = 0; j < align2->m; j++) {
                    result->indices[i][j + align1->m] = GAP_INDEX;
                }
                idx2++;
            }
        } else {
            mprintf("missed a spot...\n");
        }
    }
    return (result);
}

/* -------------------- remove_gapcols -----------------------
 * remove columns which contain nothing but gaps            */
struct pair_set *
remove_gapcols(struct pair_set *align) {
    size_t i, j;
    int allgaps, gapcols, target;
    struct pair_set *result;
    gapcols = 0;
    target = 0;
    for (i = 0; i < align->n; i++) {
        allgaps = 1;
        for (j = 0; j < align->m; j++) {
            if (align->indices[i][j] != GAP_INDEX) {
                allgaps = 0;
                break;
            }
        }
        if (allgaps) {
            gapcols++;
        }
    }
    result = E_MALLOC (sizeof (struct pair_set));
    result->m = align->m;
    result->n = align->n - gapcols;
    result->indices = i_matrix(result->n, result->m);

    for (i = 0; i < align->n; i++) {
        allgaps = 1;
        for (j = 0; j < result->m; j++) {
            if (align->indices[i][j] != GAP_INDEX) {
                allgaps = 0;
                break;
            }
        }
        if (!allgaps) {
            for (j = 0; j < result->m; j++) {
                result->indices[target][j] = align->indices[i][j];
            }
            target++;
        }
    }
    return (result);
}

/* --------------------- remove_seq ---------------------------
 * remove duplicate sequences/structures from an alignment.  */
struct pair_set *
remove_seq(struct pair_set *align, int idx) {
    struct pair_set *result;
    struct pair_set *newres;
    int src = 0;
    int i;
    int j;
    int ind = idx;
    if (ind < 0) {
        ind += align->m;
    }
    result = E_MALLOC (sizeof (struct pair_set));
    result->m = align->m - 1;
    result->n = align->n;
    result->indices = i_matrix(result->n, result->m);
    for (i = 0; i < (int) result->m; i++) {
        if (i == (int) ind) {
            src++;
        }
        for (j = 0; j < (int) result->n; j++) {
            result->indices[j][i] = align->indices[j][src];
        }
        src++;
    }
    newres = remove_gapcols(result);
    pair_set_destroy(result);
    return (newres);

    /*	return(align);*/
}

/* --------------------- remove_gaps ---------------------------
 * remove gaps from one line of an alignment.  */
struct pair_set *
remove_gaps(struct pair_set *align, int idx) {
    struct pair_set *result;
    int target = 0;
    int i;
    int j;
    int gapcount = 0;
    if (idx < 0) {
        idx += align->m;
    }
    for (i = 0; i < (int) align->n; i++) {
        if (align->indices[i][idx] == GAP_INDEX) {
            gapcount++;
        }
    }
    result = E_MALLOC (sizeof (struct pair_set));
    result->m = align->m;
    result->n = align->n - gapcount;
    result->indices = i_matrix(result->n, result->m);
    for (i = 0; i < (int) align->n; i++) {
        if (align->indices[i][idx] != GAP_INDEX) {
            for (j = 0; j < (int) result->m; j++) {
                result->indices[target][j] = align->indices[i][j];
            }
            target++;
        }

    }
    /* pair_set_destroy(align);*/
    return (result);

}

/* --------------------- split_multal --------------------------
 * split out a pairwise ailgnment from a pairset               */
struct pair_set *
split_multal(struct pair_set *align, int a, int b) {
    struct pair_set *result;
    size_t i;
    result = E_MALLOC ( sizeof( align ) );
    result->m = 2;
    result->n = align->n;
    result->indices = i_matrix(result->n, result->m);
    for (i = 0; i < result->n; i++) {
        result->indices[i][0] = align->indices[i][a];
        result->indices[i][1] = align->indices[i][b];
    }

    return (result);
}
