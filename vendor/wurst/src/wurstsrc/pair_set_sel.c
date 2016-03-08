/*
 * pair_set_sel.c
 *
 *  Created on: Sep 8, 2011
 *      Author: ibondarenko
 */

#include <ctype.h>
#include <limits.h> /* necessary for CHAR_BIT */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "bit_set.h"
#include "pair_set.h"
#include "coord.h"
#include "pair_set_sel.h"
#include "e_malloc.h"
#include "string.h"
#include "seq.h"
#include "read_seq_i.h"
#include "matrix.h"
#include "mprintf.h"

/* ---------------- pair_set_sel_geti -------------------------
 */
char *
pair_set_sel_geti(struct pair_set *align)
{
    size_t i;
    int a, b;
    char *res;
    if (align == NULL) {
        err_printf ("pair_set_sel_geti", " align is NULL.\n");
        return NULL;
    }
    res = bit_set_calloc(align->n);
    memset(res, 0, bit_set_getsize(align->n)); 
    for (i = 0; i < align->n; i++) {
        a = align->indices[i][0];
        b = align->indices[i][1];
        if ( a != GAP_INDEX && b != GAP_INDEX) {
            bit_set_bitset(res, i);
        } 
    }
    return res;
}

/*
 * @todo
 * returns a pair_set, that contains only the subset of residues
**/
/* ---------------- selected_pair_set_get ---------------------
 */
struct pair_set *
selected_pair_set_get(struct pair_set *source, char* pss)
{
    size_t i, j;
    unsigned int counter = 0;
    struct pair_set * res = E_MALLOC (sizeof (*source));
    res->m = source->m;
    res->score = source->score;
    res->smpl_score = source->smpl_score;
    res->indices = i_matrix(source->n, source->m);
    for (i = 0; i < source->n; i++) {
        if (bit_set_ifbitset(pss, i)) {
            for (j = 0; j < source->m; j++) {
                res->indices[counter][j] = source->indices[i][j];
            }
            counter++;
         }
    }
    res->n = counter;
    return res;
}

/* ---------------- pair_set_print_prepare --------------------
 */
void
pair_set_print_prepare(struct seq *s1, struct seq *s2,
                       struct pair_set *align, char* pss)
{
    size_t i;
    int a, b;
    int **pairs;
    static const char *this_sub = "pair_set_print_prepare";
    if ((s1 == NULL) || (s2 == NULL)) {
      err_printf (this_sub, "Given a NULL sequence fix me !!\n");
      return;
    }
    if (align == NULL) {
      err_printf (this_sub, "Given a NULL alignment\n");
      return;
    }
    pairs = align->indices;
    if (s1->format == THOMAS) seq_thomas2std (s1);
    if (s2->format == THOMAS) seq_thomas2std (s2);
    for (i = 0; i < align->n; i++) {
        a = pairs[i][0];
        b = pairs[i][1];
        if (bit_set_ifbitset(pss, i)) {
            s1->seq[a] = toupper(s1->seq[a]);
            s2->seq[b] = toupper(s2->seq[b]);
        }
        else {
            s1->seq[a] = tolower(s1->seq[a]);
            s2->seq[b] = tolower(s2->seq[b]);
        }
    }
}

/* ---------------- pair_set_sel_sel_print --------------------
 */
char *
pair_set_sel_print(char* pss, size_t bs_len)
{
    size_t i;
    char* res;
    res = NULL;

    res = (char*) E_MALLOC (bs_len+1);
    memset(res, 0, bs_len+1);
    for (i=0; i<bs_len; i++){
        if (bit_set_ifbitset(pss, i))
           res[i] = '1';
        else
            res[i] = '0';
    }
    return res;
}


/* ---------------- pair_set_del_maxdistance ------------------
 * @todo
 * takes an alignment and two structures and
 * deletes the flag of a pair with max distance from the source.
 **/
char *
pair_set_sel_delmaxdistance(struct coord *c1, struct coord *c2,
                            struct pair_set *align, char* source)
{

    float maxd, td;
    size_t i, maxi;
    struct RPoint dr;
    int **pairs = align->indices;

    if (!source) return source;

    maxd = 0.0;
    maxi = align->n;
    for (i=0; i<align->n; i++){
        if (bit_set_ifbitset(source, i)) {
            dr.x = c1->rp_ca[pairs[i][0]].x - c2->rp_ca[pairs[i][1]].x;
            dr.y = c1->rp_ca[pairs[i][0]].y - c2->rp_ca[pairs[i][1]].y;
            dr.z = c1->rp_ca[pairs[i][0]].z - c2->rp_ca[pairs[i][1]].z;
            td = sqrt(dr.x * dr.x + dr.y * dr.y + dr.z * dr.z);
            if (td > maxd) {
                maxi = i;
                maxd = td;
            }
        }
    }

    if (maxi < align->n) {
        bit_set_clean(source, maxi);
    }
    return source;
}

/* ---------------- printchar  --------------------------------
 */
void
printchar (int t) {
   mprintf("%d", CHAR_BIT);
}

/* ---------------- pair_set_sel_len  -------------------------
 */
size_t
pair_set_sel_len(struct pair_set *align){
    return align->n;
}
