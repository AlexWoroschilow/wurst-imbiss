/* scoranlys.c
 *  
 * 23rd January 2004
 * Routines to do with detailed analysis of threading models
 *  This groups together functions accessing data also used in wurst/src/align.c
 *  and wurst/src/rescore.c
 */

#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


#include "amino_a.h"
#include "e_malloc.h"
#include "matrix.h"
#include "mprintf.h"
#include "pair_set.h"
#include "scor_set.h"
#include "score_mat.h"
#include "score_mat_i.h"
#include "sub_mat.h"
#include "misc.h"
#include "fio.h"
#include "seq.h"
#include "read_seq_i.h"
#include "coord.h"
#include "coord_i.h"
#include "rmacros.h"
#include "mcigraph_io_i.h"
#include "scoranlys_i.h"


#if !defined (lint) && !defined (DONT_SEE_RCS)
    static const char *rcsid =
    "$Id: scoranlys.c,v 1.1 2007/09/28 16:57:16 mmundry Exp $";
#endif /* !defined (lint) && !defined (DONT_SEE_RCS) */




/* --------------- matrix calculation routines --------------- */

/* ---------------- Constants    ------------------------------
 */
#ifndef BUFSIZ
    enum {BUFSIZ = 1024};
#endif
#ifndef ZERO_AA
    enum {ZERO_AA = 20};
#endif
#include "cp_cc_allat+0.h"


/* ---------------- STANDART_BLOCK   --------------------------
 * Once I turned this into a function. It ran at the same speed
 * on a sun. I guess it made smaller code.
 */

#define STANDART_BLOCK(ri, rj, r0, n) \
            dr.x = ri->x - rj->x; \
            dr.y = ri->y - rj->y; \
            dr.z = ri->z - rj->z; \
            d = VECTOR_SQR_LENGTH(dr); \
            if (d < CUTOFF_SQR) { \
                d = sqrt(d); \
                d = (d - r0) * WIDTH_FACTOR; \
                d = 1.0 - tanh(d); \
                Xf[n] += d; \
            }


/* A simple edge-weighting scheme to see if we can find 
 * clusters of residues with high likelihood.
 */


#ifdef COMPILE_BUGS

/* ---------------- rs_pairmatrix -----------------------------
   returns a matrix containing the pairwise interactions,
   if any as measured by the forcefield.
 */
static float
ContactEmatrix( const struct coord *c, const float *P)
{
    char            *aa;
    int             i, j, aai, aaj, aai2, n, nr_aa;
    float           *cnt, *Xf;           /* scratch arrays, see note above */
    float           *Cemat; /* a symmetric score matrix containing the contact values */
    float           d, cnti, E;
    struct RPoint   dr;
    struct RPoint   *ri__N, *ri_CA, *ri_CB, *ri__C, *ri__O;
    struct RPoint   *rj__N, *rj_CA, *rj_CB, *rj__C, *rj__O;
    struct RPoint   *r__N, *r_CA, *r_CB, *r__C, *r__O;
    size_t          tmpcnt, tmpXf;
    struct seq      *s;
    const char      *this_sub = "ContactEmatrix";
    
    nr_aa = c->size;
    s = c->seq;
    seq_std2thomas (s);  /* Force Thomas style names for amino acids */
    aa = s->seq;
    r__N = c->rp_n;
    r_CA = c->rp_ca;
    r_CB = c->rp_cb;
    r__C = c->rp_c;
    r__O = c->rp_o;
    
    cnt = (float *) E_MALLOC ( tmpcnt = (nr_aa * sizeof (float)));
    memset (cnt, 0, tmpcnt);
    Xf =  (float *) E_MALLOC ( tmpXf = (NR_PARAM * sizeof (float)));
    memset (Xf, 0, tmpXf);


    for (i = 0; i < nr_aa-2; i++) {
        j = i + 2;
        aai = aa[i];
        aaj = aa[j];
        if ((aai != ZERO_AA) && (aaj != ZERO_AA)) {
            ri__N = &r__N[i];
            ri_CA = &r_CA[i];
            ri_CB = &r_CB[i];
            ri__C = &r__C[i];
            ri__O = &r__O[i];
            rj__N = &r__N[j];
            rj_CA = &r_CA[j];
            rj_CB = &r_CB[j];
            rj__C = &r__C[j];
            rj__O = &r__O[j];

            STANDART_BLOCK(ri__N, rj__N, R01__N__N,     START1__N__N)
            STANDART_BLOCK(ri__N, rj_CA, R01__N_CA,     START1__N_CA)
            STANDART_BLOCK(ri__N, rj_CB, R01__N_CB, aaj+START1__N_CB)
            STANDART_BLOCK(ri__N, rj__C, R01__N__C,     START1__N__C)
            STANDART_BLOCK(ri__N, rj__O, R01__N__O,     START1__N__O)

            STANDART_BLOCK(ri_CA, rj__N, R01_CA__N,     START1_CA__N)
            STANDART_BLOCK(ri_CA, rj_CA, R01_CA_CA,     START1_CA_CA)
            STANDART_BLOCK(ri_CA, rj_CB, R01_CA_CB, aaj+START1_CA_CB)
            STANDART_BLOCK(ri_CA, rj__C, R01_CA__C,     START1_CA__C)
            STANDART_BLOCK(ri_CA, rj__O, R01_CA__O,     START1_CA__O)

            STANDART_BLOCK(ri_CB, rj__N, R01_CB__N, aai+START1_CB__N)
            STANDART_BLOCK(ri_CB, rj_CA, R01_CB_CA, aai+START1_CB_CA)
            if (aai < aaj)
                n = aai + aaj*(aaj+1)/2;
            else
            n = aaj + aai*(aai+1)/2;
            STANDART_BLOCK(ri_CB, rj_CB, R01_CB_CB,   n+START1_CB_CB)
            STANDART_BLOCK(ri_CB, rj__C, R01_CB__C, aai+START1_CB__C)
            STANDART_BLOCK(ri_CB, rj__O, R01_CB__O, aai+START1_CB__O)

            STANDART_BLOCK(ri__C, rj__N, R01__C__N,     START1__C__N)
            STANDART_BLOCK(ri__C, rj_CA, R01__C_CA,     START1__C_CA)
            STANDART_BLOCK(ri__C, rj_CB, R01__C_CB, aaj+START1__C_CB)
            STANDART_BLOCK(ri__C, rj__C, R01__C__C,     START1__C__C)
            STANDART_BLOCK(ri__C, rj__O, R01__C__O,     START1__C__O)

            STANDART_BLOCK(ri__O, rj__N, R01__O__N,     START1__O__N)
            STANDART_BLOCK(ri__O, rj_CA, R01__O_CA,     START1__O_CA)
            STANDART_BLOCK(ri__O, rj_CB, R01__O_CB, aaj+START1__O_CB)
            STANDART_BLOCK(ri__O, rj__C, R01__O__C,     START1__O__C)
            STANDART_BLOCK(ri__O, rj__O, R01__O__O,     START1__O__O)
        }
    }

    for (i = 0; i < nr_aa-3; i++) {
        j = i + 3;
        aai = aa[i];
        aaj = aa[j];
        if ((aai != ZERO_AA) && (aaj != ZERO_AA)) {
            ri__N = &r__N[i];
            ri_CA = &r_CA[i];
            ri_CB = &r_CB[i];
            ri__C = &r__C[i];
            ri__O = &r__O[i];
            rj__N = &r__N[j];
            rj_CA = &r_CA[j];
            rj_CB = &r_CB[j];
            rj__C = &r__C[j];
            rj__O = &r__O[j];

            STANDART_BLOCK(ri__N, rj__N, R02__N__N,     START2__N__N)
            STANDART_BLOCK(ri__N, rj_CA, R02__N_CA,     START2__N_CA)
            STANDART_BLOCK(ri__N, rj_CB, R02__N_CB, aaj+START2__N_CB)
            STANDART_BLOCK(ri__N, rj__C, R02__N__C,     START2__N__C)
            STANDART_BLOCK(ri__N, rj__O, R02__N__O,     START2__N__O)

            STANDART_BLOCK(ri_CA, rj__N, R02_CA__N,     START2_CA__N)
            STANDART_BLOCK(ri_CA, rj_CA, R02_CA_CA,     START2_CA_CA)
            STANDART_BLOCK(ri_CA, rj_CB, R02_CA_CB, aaj+START2_CA_CB)
            STANDART_BLOCK(ri_CA, rj__C, R02_CA__C,     START2_CA__C)
            STANDART_BLOCK(ri_CA, rj__O, R02_CA__O,     START2_CA__O)

            STANDART_BLOCK(ri_CB, rj__N, R02_CB__N, aai+START2_CB__N)
            STANDART_BLOCK(ri_CB, rj_CA, R02_CB_CA, aai+START2_CB_CA)
            if (aai < aaj)
                n = aai + aaj*(aaj+1)/2;
            else
            n = aaj + aai*(aai+1)/2;
            STANDART_BLOCK(ri_CB, rj_CB, R02_CB_CB,   n+START2_CB_CB)
            STANDART_BLOCK(ri_CB, rj__C, R02_CB__C, aai+START2_CB__C)
            STANDART_BLOCK(ri_CB, rj__O, R02_CB__O, aai+START2_CB__O)

            STANDART_BLOCK(ri__C, rj__N, R02__C__N,     START2__C__N)
            STANDART_BLOCK(ri__C, rj_CA, R02__C_CA,     START2__C_CA)
            STANDART_BLOCK(ri__C, rj_CB, R02__C_CB, aaj+START2__C_CB)
            STANDART_BLOCK(ri__C, rj__C, R02__C__C,     START2__C__C)
            STANDART_BLOCK(ri__C, rj__O, R02__C__O,     START2__C__O)

            STANDART_BLOCK(ri__O, rj__N, R02__O__N,     START2__O__N)
            STANDART_BLOCK(ri__O, rj_CA, R02__O_CA,     START2__O_CA)
            STANDART_BLOCK(ri__O, rj_CB, R02__O_CB, aaj+START2__O_CB)
            STANDART_BLOCK(ri__O, rj__C, R02__O__C,     START2__O__C)
            STANDART_BLOCK(ri__O, rj__O, R02__O__O,     START2__O__O)
        }
    }

    for (i = 0; i < nr_aa-4; i++) {
        aai = aa[i];
        if (aai != ZERO_AA) {
            aai2 = aai * (aai + 1) / 2;
            ri__N = &r__N[i];
            ri_CA = &r_CA[i];
            ri_CB = &r_CB[i];
            ri__C = &r__C[i];
            ri__O = &r__O[i];
            cnti = 0.0;
            for (j = i+4; j < nr_aa; j++) {
                aaj = aa[j];
                if (aaj != ZERO_AA) {
                    rj_CA = &r_CA[j];
                    dr.x = ri_CA->x - rj_CA->x;      /* CA - CA interaction */
                    dr.y = ri_CA->y - rj_CA->y;
                    dr.z = ri_CA->z - rj_CA->z;
                    d = VECTOR_SQR_LENGTH(dr);
                    if (d < CA_CA_CUTOFF_SQR) {
                        if (d < CUTOFF_SQR) {
                            d = sqrt(d);
                            d = (d - R03_CA_CA) * WIDTH_FACTOR;
                            d = 1.0 - tanh(d);
                            cnti += d;
                            cnt[j] += d;
                            Xf[START3_CA_CA] += d;
                        }

                        rj__N = &r__N[j];
                        rj_CB = &r_CB[j];
                        rj__C = &r__C[j];
                        rj__O = &r__O[j];
                        STANDART_BLOCK(ri__N, rj__N, R03__N__N,     START3__N__N)
                        STANDART_BLOCK(ri__N, rj_CA, R03__N_CA,     START3__N_CA)
                        STANDART_BLOCK(ri__N, rj_CB, R03__N_CB, aaj+START3__N_CB)
                        STANDART_BLOCK(ri__N, rj__C, R03__N__C,     START3__N__C)
                        STANDART_BLOCK(ri__N, rj__O, R03__N__O,     START3__N__O)

                       STANDART_BLOCK(ri_CA, rj__N, R03_CA__N,     START3_CA__N)
/*                        CA--CA see above  */
                       STANDART_BLOCK(ri_CA, rj_CB, R03_CA_CB, aaj+START3_CA_CB)
                       STANDART_BLOCK(ri_CA, rj__C, R03_CA__C,     START3_CA__C)
                       STANDART_BLOCK(ri_CA, rj__O, R03_CA__O,     START3_CA__O)

                       STANDART_BLOCK(ri_CB, rj__N, R03_CB__N, aai+START3_CB__N)
                       STANDART_BLOCK(ri_CB, rj_CA, R03_CB_CA, aai+START3_CB_CA)
                       if (aai < aaj)
                           n = aai + aaj*(aaj+1)/2;
                       else
                           n = aaj + aai2;
                       STANDART_BLOCK(ri_CB, rj_CB, R03_CB_CB,   n+START3_CB_CB)
                       STANDART_BLOCK(ri_CB, rj__C, R03_CB__C, aai+START3_CB__C)
                       STANDART_BLOCK(ri_CB, rj__O, R03_CB__O, aai+START3_CB__O)

                       STANDART_BLOCK(ri__C, rj__N, R03__C__N,     START3__C__N)
                       STANDART_BLOCK(ri__C, rj_CA, R03__C_CA,     START3__C_CA)
                       STANDART_BLOCK(ri__C, rj_CB, R03__C_CB, aaj+START3__C_CB)
                       STANDART_BLOCK(ri__C, rj__C, R03__C__C,     START3__C__C)
                       STANDART_BLOCK(ri__C, rj__O, R03__C__O,     START3__C__O)

                       STANDART_BLOCK(ri__O, rj__N, R03__O__N,     START3__O__N)
                       STANDART_BLOCK(ri__O, rj_CA, R03__O_CA,     START3__O_CA)
                       STANDART_BLOCK(ri__O, rj_CB, R03__O_CB, aaj+START3__O_CB)
                       STANDART_BLOCK(ri__O, rj__C, R03__O__C,     START3__O__C)
                       STANDART_BLOCK(ri__O, rj__O, R03__O__O,     START3__O__O)

                    } /* if (CA_CA_distance < CA_CA_CUTOFF_SQR) */
                }     /* if (aaj != ZERO_AA) */
            }     /* for all j aa */
            cnt[i] += cnti;
        }      /* if (aai != ZERO_AA) */
    }         /* for all i aa */

    for (i = 0; i < nr_aa; i++) {
        d = (cnt[i] - NR_NEIGHBOUR);
        d = 1.0 - tanh(d);
        n = START4_NEIGHBOUR + aa[i];
        if (n >= NR_PARAM) {
            err_printf (this_sub, "n: %d >= NR_PARAM: %d. Disaster.\n", n, NR_PARAM);
            err_printf (this_sub, "Stopping at %s: %d\n", __FILE__, __LINE__);
            exit (EXIT_FAILURE);
        }
        Xf[n] += d;
    }

    E = 0.0;
    for (i = 0; i < NR_PARAM; i++)
        E += Xf[i] * P[i];

    free (Xf);
    free (cnt);

    return -E;
}

#endif

/* ------------scor_set_rs --------------------------------------
   returns the vector of local scores from all
   enumerating over pairwise interactions,
   if any, as measured by the forcefield.
 */

/* STANDART_BLOCK
   (sic) assumes that we are in an i,j interaction
   loop and building a score set awght 
   Xf is still maintained because of the neighbour counting term
*/
#define STANDARDT_BLOCK(ri, rj, r0, n) \
            dr.x = ri->x - rj->x; \
            dr.y = ri->y - rj->y; \
            dr.z = ri->z - rj->z; \
            d = VECTOR_SQR_LENGTH(dr); \
            if (d < CUTOFF_SQR) { \
                d = sqrt(d); \
                d = (d - r0) * WIDTH_FACTOR; \
                d = 1.0 - tanh(d); \
                awght->scores[i]+=-P[n]* d; \
                awght->scores[j]+=-P[n]* d; \
                Xf[n] += d; \
            }

struct scor_set *
scor_set_rs( struct coord *c, const float *P)
{
    char            *aa;
    int             i, j, aai, aaj, aai2, n, nr_aa;
    float           *cnt, *Xf;           /* scratch arrays, see note above */
    float           d, cnti;
    struct scor_set *awght;
    struct RPoint   dr;
    struct RPoint   *ri__N, *ri_CA, *ri_CB, *ri__C, *ri__O;
    struct RPoint   *rj__N, *rj_CA, *rj_CB, *rj__C, *rj__O;
    struct RPoint   *r__N, *r_CA, *r_CB, *r__C, *r__O;
    size_t          tmpcnt, tmpXf;
    struct seq      *s;
    const char      *this_sub = "scor_set_rs";
    coord_a_2_nm(c);
    nr_aa = c->size;
    s = c->seq;
    awght = E_MALLOC(sizeof(*awght));
    awght->n = nr_aa;
    awght->scores = E_MALLOC(sizeof(*awght->scores)*nr_aa);
    for (i=0; i<nr_aa; i++)
        awght->scores[i]=0.0;
    seq_std2thomas (s);  /* Force Thomas style names for amino acids */
    aa = s->seq;
    r__N = c->rp_n;
    r_CA = c->rp_ca;
    r_CB = c->rp_cb;
    r__C = c->rp_c;
    r__O = c->rp_o;
    
    cnt = (float *) E_MALLOC ( tmpcnt = (nr_aa * sizeof (float)));
    memset (cnt, 0, tmpcnt);
    Xf =  (float *) E_MALLOC ( tmpXf = (NR_PARAM * sizeof (float)));
    memset (Xf, 0, tmpXf);


    for (i = 0; i < nr_aa-2; i++) {
        j = i + 2;
        aai = aa[i];
        aaj = aa[j];
        if ((aai != ZERO_AA) && (aaj != ZERO_AA)) {
            ri__N = &r__N[i];
            ri_CA = &r_CA[i];
            ri_CB = &r_CB[i];
            ri__C = &r__C[i];
            ri__O = &r__O[i];
            rj__N = &r__N[j];
            rj_CA = &r_CA[j];
            rj_CB = &r_CB[j];
            rj__C = &r__C[j];
            rj__O = &r__O[j];

            STANDARDT_BLOCK(ri__N, rj__N, R01__N__N,     START1__N__N)
            STANDARDT_BLOCK(ri__N, rj_CA, R01__N_CA,     START1__N_CA)
            STANDARDT_BLOCK(ri__N, rj_CB, R01__N_CB, aaj+START1__N_CB)
            STANDARDT_BLOCK(ri__N, rj__C, R01__N__C,     START1__N__C)
            STANDARDT_BLOCK(ri__N, rj__O, R01__N__O,     START1__N__O)

            STANDARDT_BLOCK(ri_CA, rj__N, R01_CA__N,     START1_CA__N)
            STANDARDT_BLOCK(ri_CA, rj_CA, R01_CA_CA,     START1_CA_CA)
            STANDARDT_BLOCK(ri_CA, rj_CB, R01_CA_CB, aaj+START1_CA_CB)
            STANDARDT_BLOCK(ri_CA, rj__C, R01_CA__C,     START1_CA__C)
            STANDARDT_BLOCK(ri_CA, rj__O, R01_CA__O,     START1_CA__O)

            STANDARDT_BLOCK(ri_CB, rj__N, R01_CB__N, aai+START1_CB__N)
            STANDARDT_BLOCK(ri_CB, rj_CA, R01_CB_CA, aai+START1_CB_CA)
            if (aai < aaj)
                n = aai + aaj*(aaj+1)/2;
            else
            n = aaj + aai*(aai+1)/2;
            STANDARDT_BLOCK(ri_CB, rj_CB, R01_CB_CB,   n+START1_CB_CB)
            STANDARDT_BLOCK(ri_CB, rj__C, R01_CB__C, aai+START1_CB__C)
            STANDARDT_BLOCK(ri_CB, rj__O, R01_CB__O, aai+START1_CB__O)

            STANDARDT_BLOCK(ri__C, rj__N, R01__C__N,     START1__C__N)
            STANDARDT_BLOCK(ri__C, rj_CA, R01__C_CA,     START1__C_CA)
            STANDARDT_BLOCK(ri__C, rj_CB, R01__C_CB, aaj+START1__C_CB)
            STANDARDT_BLOCK(ri__C, rj__C, R01__C__C,     START1__C__C)
            STANDARDT_BLOCK(ri__C, rj__O, R01__C__O,     START1__C__O)

            STANDARDT_BLOCK(ri__O, rj__N, R01__O__N,     START1__O__N)
            STANDARDT_BLOCK(ri__O, rj_CA, R01__O_CA,     START1__O_CA)
            STANDARDT_BLOCK(ri__O, rj_CB, R01__O_CB, aaj+START1__O_CB)
            STANDARDT_BLOCK(ri__O, rj__C, R01__O__C,     START1__O__C)
            STANDARDT_BLOCK(ri__O, rj__O, R01__O__O,     START1__O__O)
        }
    }

    for (i = 0; i < nr_aa-3; i++) {
        j = i + 3;
        aai = aa[i];
        aaj = aa[j];
        if ((aai != ZERO_AA) && (aaj != ZERO_AA)) {
            ri__N = &r__N[i];
            ri_CA = &r_CA[i];
            ri_CB = &r_CB[i];
            ri__C = &r__C[i];
            ri__O = &r__O[i];
            rj__N = &r__N[j];
            rj_CA = &r_CA[j];
            rj_CB = &r_CB[j];
            rj__C = &r__C[j];
            rj__O = &r__O[j];

            STANDARDT_BLOCK(ri__N, rj__N, R02__N__N,     START2__N__N)
            STANDARDT_BLOCK(ri__N, rj_CA, R02__N_CA,     START2__N_CA)
            STANDARDT_BLOCK(ri__N, rj_CB, R02__N_CB, aaj+START2__N_CB)
            STANDARDT_BLOCK(ri__N, rj__C, R02__N__C,     START2__N__C)
            STANDARDT_BLOCK(ri__N, rj__O, R02__N__O,     START2__N__O)

            STANDARDT_BLOCK(ri_CA, rj__N, R02_CA__N,     START2_CA__N)
            STANDARDT_BLOCK(ri_CA, rj_CA, R02_CA_CA,     START2_CA_CA)
            STANDARDT_BLOCK(ri_CA, rj_CB, R02_CA_CB, aaj+START2_CA_CB)
            STANDARDT_BLOCK(ri_CA, rj__C, R02_CA__C,     START2_CA__C)
            STANDARDT_BLOCK(ri_CA, rj__O, R02_CA__O,     START2_CA__O)

            STANDARDT_BLOCK(ri_CB, rj__N, R02_CB__N, aai+START2_CB__N)
            STANDARDT_BLOCK(ri_CB, rj_CA, R02_CB_CA, aai+START2_CB_CA)
            if (aai < aaj)
                n = aai + aaj*(aaj+1)/2;
            else
            n = aaj + aai*(aai+1)/2;
            STANDARDT_BLOCK(ri_CB, rj_CB, R02_CB_CB,   n+START2_CB_CB)
            STANDARDT_BLOCK(ri_CB, rj__C, R02_CB__C, aai+START2_CB__C)
            STANDARDT_BLOCK(ri_CB, rj__O, R02_CB__O, aai+START2_CB__O)

            STANDARDT_BLOCK(ri__C, rj__N, R02__C__N,     START2__C__N)
            STANDARDT_BLOCK(ri__C, rj_CA, R02__C_CA,     START2__C_CA)
            STANDARDT_BLOCK(ri__C, rj_CB, R02__C_CB, aaj+START2__C_CB)
            STANDARDT_BLOCK(ri__C, rj__C, R02__C__C,     START2__C__C)
            STANDARDT_BLOCK(ri__C, rj__O, R02__C__O,     START2__C__O)

            STANDARDT_BLOCK(ri__O, rj__N, R02__O__N,     START2__O__N)
            STANDARDT_BLOCK(ri__O, rj_CA, R02__O_CA,     START2__O_CA)
            STANDARDT_BLOCK(ri__O, rj_CB, R02__O_CB, aaj+START2__O_CB)
            STANDARDT_BLOCK(ri__O, rj__C, R02__O__C,     START2__O__C)
            STANDARDT_BLOCK(ri__O, rj__O, R02__O__O,     START2__O__O)
        }
    }

    for (i = 0; i < nr_aa-4; i++) {
        aai = aa[i];
        if (aai != ZERO_AA) {
            aai2 = aai * (aai + 1) / 2;
            ri__N = &r__N[i];
            ri_CA = &r_CA[i];
            ri_CB = &r_CB[i];
            ri__C = &r__C[i];
            ri__O = &r__O[i];
            cnti = 0.0;
            for (j = i+4; j < nr_aa; j++) {
                aaj = aa[j];
                if (aaj != ZERO_AA) {
                    rj_CA = &r_CA[j];
                    dr.x = ri_CA->x - rj_CA->x;      /* CA - CA interaction */
                    dr.y = ri_CA->y - rj_CA->y;
                    dr.z = ri_CA->z - rj_CA->z;
                    d = VECTOR_SQR_LENGTH(dr);
                    if (d < CA_CA_CUTOFF_SQR) {
                        if (d < CUTOFF_SQR) {
                            d = sqrt(d);
                            d = (d - R03_CA_CA) * WIDTH_FACTOR;
                            d = 1.0 - tanh(d);
                            cnti += d;
                            cnt[j] += d;
                            Xf[START3_CA_CA] += d;
                            awght->scores[i]+=-P[START3_CA_CA]*d;
                            awght->scores[j]+=-P[START3_CA_CA]*d;
                            
                        }

                        rj__N = &r__N[j];
                        rj_CB = &r_CB[j];
                        rj__C = &r__C[j];
                        rj__O = &r__O[j];
                        STANDARDT_BLOCK(ri__N, rj__N, R03__N__N,     START3__N__N)
                        STANDARDT_BLOCK(ri__N, rj_CA, R03__N_CA,     START3__N_CA)
                        STANDARDT_BLOCK(ri__N, rj_CB, R03__N_CB, aaj+START3__N_CB)
                        STANDARDT_BLOCK(ri__N, rj__C, R03__N__C,     START3__N__C)
                        STANDARDT_BLOCK(ri__N, rj__O, R03__N__O,     START3__N__O)

                       STANDARDT_BLOCK(ri_CA, rj__N, R03_CA__N,     START3_CA__N)
/*                        CA--CA see above  */
                       STANDARDT_BLOCK(ri_CA, rj_CB, R03_CA_CB, aaj+START3_CA_CB)
                       STANDARDT_BLOCK(ri_CA, rj__C, R03_CA__C,     START3_CA__C)
                       STANDARDT_BLOCK(ri_CA, rj__O, R03_CA__O,     START3_CA__O)

                       STANDARDT_BLOCK(ri_CB, rj__N, R03_CB__N, aai+START3_CB__N)
                       STANDARDT_BLOCK(ri_CB, rj_CA, R03_CB_CA, aai+START3_CB_CA)
                       if (aai < aaj)
                           n = aai + aaj*(aaj+1)/2;
                       else
                           n = aaj + aai2;
                       STANDARDT_BLOCK(ri_CB, rj_CB, R03_CB_CB,   n+START3_CB_CB)
                       STANDARDT_BLOCK(ri_CB, rj__C, R03_CB__C, aai+START3_CB__C)
                       STANDARDT_BLOCK(ri_CB, rj__O, R03_CB__O, aai+START3_CB__O)

                       STANDARDT_BLOCK(ri__C, rj__N, R03__C__N,     START3__C__N)
                       STANDARDT_BLOCK(ri__C, rj_CA, R03__C_CA,     START3__C_CA)
                       STANDARDT_BLOCK(ri__C, rj_CB, R03__C_CB, aaj+START3__C_CB)
                       STANDARDT_BLOCK(ri__C, rj__C, R03__C__C,     START3__C__C)
                       STANDARDT_BLOCK(ri__C, rj__O, R03__C__O,     START3__C__O)

                       STANDARDT_BLOCK(ri__O, rj__N, R03__O__N,     START3__O__N)
                       STANDARDT_BLOCK(ri__O, rj_CA, R03__O_CA,     START3__O_CA)
                       STANDARDT_BLOCK(ri__O, rj_CB, R03__O_CB, aaj+START3__O_CB)
                       STANDARDT_BLOCK(ri__O, rj__C, R03__O__C,     START3__O__C)
                       STANDARDT_BLOCK(ri__O, rj__O, R03__O__O,     START3__O__O)

                    } /* if (CA_CA_distance < CA_CA_CUTOFF_SQR) */
                }     /* if (aaj != ZERO_AA) */
            }     /* for all j aa */
            cnt[i] += cnti;
        }      /* if (aai != ZERO_AA) */
    }         /* for all i aa */

    for (i = 0; i < nr_aa; i++) {
        d = (cnt[i] - NR_NEIGHBOUR);
        d = 1.0 - tanh(d);
        n = START4_NEIGHBOUR + aa[i];
        if (n >= NR_PARAM) {
            err_printf (this_sub, "n: %d >= NR_PARAM: %d. Disaster.\n", n, NR_PARAM);
            err_printf (this_sub, "Stopping at %s: %d\n", __FILE__, __LINE__);
            exit (EXIT_FAILURE);
        }
        awght->scores[i] += -Xf[n]*P[n];
    }
    free (cnt);
    free(Xf);
    coord_nm_2_a(c);

    return awght;
}

/* STANDARXT_BLOCK
   (e.sic) assumes that we are in an i,j interaction
   loop and building a score matrix awght 
   And what do we do with the neighbour count ?
*/
#define STANDARXT_BLOCK(ri, rj, r0, n) \
            dr.x = ri->x - rj->x; \
            dr.y = ri->y - rj->y; \
            dr.z = ri->z - rj->z; \
            d = VECTOR_SQR_LENGTH(dr); \
            if (d < CUTOFF_SQR) { \
                d = sqrt(d); \
                d = (d - r0) * WIDTH_FACTOR; \
                d = 1.0 - tanh(d); \
                awght[i][j]-=P[n]* d; \
                awght[j][i]-=P[n]* d; \
                Xf[n] += d; \
            }



static float **
scor_contact_rs( struct coord *c, const float *P)
{
    char            *aa;
    int             i, j, aai, aaj, aai2, n, nr_aa;
    float           *cnt, *Xf;           /* scratch arrays, see note above */
    float           **awght; /* The weight matrix we are calculating */
    float           d, cnti;
    struct RPoint   dr;
    struct RPoint   *ri__N, *ri_CA, *ri_CB, *ri__C, *ri__O;
    struct RPoint   *rj__N, *rj_CA, *rj_CB, *rj__C, *rj__O;
    struct RPoint   *r__N, *r_CA, *r_CB, *r__C, *r__O;
    size_t          tmpcnt, tmpXf;
    struct seq      *s;
    const char      *this_sub = "scor_contact_rs";
    coord_a_2_nm(c);
    nr_aa = c->size;
    s = c->seq;
    awght = f_matrix(nr_aa, nr_aa);
    for (i=0; i<nr_aa; i++)
        for (j=0; j<nr_aa; j++)
            awght[i][j]=0.0;

    seq_std2thomas (s);  /* Force Thomas style names for amino acids */
    aa = s->seq;
    r__N = c->rp_n;
    r_CA = c->rp_ca;
    r_CB = c->rp_cb;
    r__C = c->rp_c;
    r__O = c->rp_o;
    
    cnt = (float *) E_MALLOC ( tmpcnt = (nr_aa * sizeof (float)));
    memset (cnt, 0, tmpcnt);
    Xf =  (float *) E_MALLOC ( tmpXf = (NR_PARAM * sizeof (float)));
    memset (Xf, 0, tmpXf);


    for (i = 0; i < nr_aa-2; i++) {
        j = i + 2;
        aai = aa[i];
        aaj = aa[j];
        if ((aai != ZERO_AA) && (aaj != ZERO_AA)) {
            ri__N = &r__N[i];
            ri_CA = &r_CA[i];
            ri_CB = &r_CB[i];
            ri__C = &r__C[i];
            ri__O = &r__O[i];
            rj__N = &r__N[j];
            rj_CA = &r_CA[j];
            rj_CB = &r_CB[j];
            rj__C = &r__C[j];
            rj__O = &r__O[j];

            STANDARXT_BLOCK(ri__N, rj__N, R01__N__N,     START1__N__N)
            STANDARXT_BLOCK(ri__N, rj_CA, R01__N_CA,     START1__N_CA)
            STANDARXT_BLOCK(ri__N, rj_CB, R01__N_CB, aaj+START1__N_CB)
            STANDARXT_BLOCK(ri__N, rj__C, R01__N__C,     START1__N__C)
            STANDARXT_BLOCK(ri__N, rj__O, R01__N__O,     START1__N__O)

            STANDARXT_BLOCK(ri_CA, rj__N, R01_CA__N,     START1_CA__N)
            STANDARXT_BLOCK(ri_CA, rj_CA, R01_CA_CA,     START1_CA_CA)
            STANDARXT_BLOCK(ri_CA, rj_CB, R01_CA_CB, aaj+START1_CA_CB)
            STANDARXT_BLOCK(ri_CA, rj__C, R01_CA__C,     START1_CA__C)
            STANDARXT_BLOCK(ri_CA, rj__O, R01_CA__O,     START1_CA__O)

            STANDARXT_BLOCK(ri_CB, rj__N, R01_CB__N, aai+START1_CB__N)
            STANDARXT_BLOCK(ri_CB, rj_CA, R01_CB_CA, aai+START1_CB_CA)
            if (aai < aaj)
                n = aai + aaj*(aaj+1)/2;
            else
            n = aaj + aai*(aai+1)/2;
            STANDARXT_BLOCK(ri_CB, rj_CB, R01_CB_CB,   n+START1_CB_CB)
            STANDARXT_BLOCK(ri_CB, rj__C, R01_CB__C, aai+START1_CB__C)
            STANDARXT_BLOCK(ri_CB, rj__O, R01_CB__O, aai+START1_CB__O)

            STANDARXT_BLOCK(ri__C, rj__N, R01__C__N,     START1__C__N)
            STANDARXT_BLOCK(ri__C, rj_CA, R01__C_CA,     START1__C_CA)
            STANDARXT_BLOCK(ri__C, rj_CB, R01__C_CB, aaj+START1__C_CB)
            STANDARXT_BLOCK(ri__C, rj__C, R01__C__C,     START1__C__C)
            STANDARXT_BLOCK(ri__C, rj__O, R01__C__O,     START1__C__O)

            STANDARXT_BLOCK(ri__O, rj__N, R01__O__N,     START1__O__N)
            STANDARXT_BLOCK(ri__O, rj_CA, R01__O_CA,     START1__O_CA)
            STANDARXT_BLOCK(ri__O, rj_CB, R01__O_CB, aaj+START1__O_CB)
            STANDARXT_BLOCK(ri__O, rj__C, R01__O__C,     START1__O__C)
            STANDARXT_BLOCK(ri__O, rj__O, R01__O__O,     START1__O__O)
        }
    }

    for (i = 0; i < nr_aa-3; i++) {
        j = i + 3;
        aai = aa[i];
        aaj = aa[j];
        if ((aai != ZERO_AA) && (aaj != ZERO_AA)) {
            ri__N = &r__N[i];
            ri_CA = &r_CA[i];
            ri_CB = &r_CB[i];
            ri__C = &r__C[i];
            ri__O = &r__O[i];
            rj__N = &r__N[j];
            rj_CA = &r_CA[j];
            rj_CB = &r_CB[j];
            rj__C = &r__C[j];
            rj__O = &r__O[j];

            STANDARXT_BLOCK(ri__N, rj__N, R02__N__N,     START2__N__N)
            STANDARXT_BLOCK(ri__N, rj_CA, R02__N_CA,     START2__N_CA)
            STANDARXT_BLOCK(ri__N, rj_CB, R02__N_CB, aaj+START2__N_CB)
            STANDARXT_BLOCK(ri__N, rj__C, R02__N__C,     START2__N__C)
            STANDARXT_BLOCK(ri__N, rj__O, R02__N__O,     START2__N__O)

            STANDARXT_BLOCK(ri_CA, rj__N, R02_CA__N,     START2_CA__N)
            STANDARXT_BLOCK(ri_CA, rj_CA, R02_CA_CA,     START2_CA_CA)
            STANDARXT_BLOCK(ri_CA, rj_CB, R02_CA_CB, aaj+START2_CA_CB)
            STANDARXT_BLOCK(ri_CA, rj__C, R02_CA__C,     START2_CA__C)
            STANDARXT_BLOCK(ri_CA, rj__O, R02_CA__O,     START2_CA__O)

            STANDARXT_BLOCK(ri_CB, rj__N, R02_CB__N, aai+START2_CB__N)
            STANDARXT_BLOCK(ri_CB, rj_CA, R02_CB_CA, aai+START2_CB_CA)
            if (aai < aaj)
                n = aai + aaj*(aaj+1)/2;
            else
            n = aaj + aai*(aai+1)/2;
            STANDARXT_BLOCK(ri_CB, rj_CB, R02_CB_CB,   n+START2_CB_CB)
            STANDARXT_BLOCK(ri_CB, rj__C, R02_CB__C, aai+START2_CB__C)
            STANDARXT_BLOCK(ri_CB, rj__O, R02_CB__O, aai+START2_CB__O)

            STANDARXT_BLOCK(ri__C, rj__N, R02__C__N,     START2__C__N)
            STANDARXT_BLOCK(ri__C, rj_CA, R02__C_CA,     START2__C_CA)
            STANDARXT_BLOCK(ri__C, rj_CB, R02__C_CB, aaj+START2__C_CB)
            STANDARXT_BLOCK(ri__C, rj__C, R02__C__C,     START2__C__C)
            STANDARXT_BLOCK(ri__C, rj__O, R02__C__O,     START2__C__O)

            STANDARXT_BLOCK(ri__O, rj__N, R02__O__N,     START2__O__N)
            STANDARXT_BLOCK(ri__O, rj_CA, R02__O_CA,     START2__O_CA)
            STANDARXT_BLOCK(ri__O, rj_CB, R02__O_CB, aaj+START2__O_CB)
            STANDARXT_BLOCK(ri__O, rj__C, R02__O__C,     START2__O__C)
            STANDARXT_BLOCK(ri__O, rj__O, R02__O__O,     START2__O__O)
        }
    }

    for (i = 0; i < nr_aa-4; i++) {
        aai = aa[i];
        if (aai != ZERO_AA) {
            aai2 = aai * (aai + 1) / 2;
            ri__N = &r__N[i];
            ri_CA = &r_CA[i];
            ri_CB = &r_CB[i];
            ri__C = &r__C[i];
            ri__O = &r__O[i];
            cnti = 0.0;
            for (j = i+4; j < nr_aa; j++) {
                aaj = aa[j];
                if (aaj != ZERO_AA) {
                    rj_CA = &r_CA[j];
                    dr.x = ri_CA->x - rj_CA->x;      /* CA - CA interaction */
                    dr.y = ri_CA->y - rj_CA->y;
                    dr.z = ri_CA->z - rj_CA->z;
                    d = VECTOR_SQR_LENGTH(dr);
                    if (d < CA_CA_CUTOFF_SQR) {
                        if (d < CUTOFF_SQR) {
                            d = sqrt(d);
                            d = (d - R03_CA_CA) * WIDTH_FACTOR;
                            d = 1.0 - tanh(d);
                            cnti += d;
                            cnt[j] += d;
                            Xf[START3_CA_CA] += d;
                            awght[i][j]-=P[n]*d;
                            awght[j][i]-=P[n]*d;
                        }

                        rj__N = &r__N[j];
                        rj_CB = &r_CB[j];
                        rj__C = &r__C[j];
                        rj__O = &r__O[j];
                        STANDARXT_BLOCK(ri__N, rj__N, R03__N__N,     START3__N__N)
                        STANDARXT_BLOCK(ri__N, rj_CA, R03__N_CA,     START3__N_CA)
                        STANDARXT_BLOCK(ri__N, rj_CB, R03__N_CB, aaj+START3__N_CB)
                        STANDARXT_BLOCK(ri__N, rj__C, R03__N__C,     START3__N__C)
                        STANDARXT_BLOCK(ri__N, rj__O, R03__N__O,     START3__N__O)

                       STANDARXT_BLOCK(ri_CA, rj__N, R03_CA__N,     START3_CA__N)
/*                        CA--CA see above  */
                       STANDARXT_BLOCK(ri_CA, rj_CB, R03_CA_CB, aaj+START3_CA_CB)
                       STANDARXT_BLOCK(ri_CA, rj__C, R03_CA__C,     START3_CA__C)
                       STANDARXT_BLOCK(ri_CA, rj__O, R03_CA__O,     START3_CA__O)

                       STANDARXT_BLOCK(ri_CB, rj__N, R03_CB__N, aai+START3_CB__N)
                       STANDARXT_BLOCK(ri_CB, rj_CA, R03_CB_CA, aai+START3_CB_CA)
                       if (aai < aaj)
                           n = aai + aaj*(aaj+1)/2;
                       else
                           n = aaj + aai2;
                       STANDARXT_BLOCK(ri_CB, rj_CB, R03_CB_CB,   n+START3_CB_CB)
                       STANDARXT_BLOCK(ri_CB, rj__C, R03_CB__C, aai+START3_CB__C)
                       STANDARXT_BLOCK(ri_CB, rj__O, R03_CB__O, aai+START3_CB__O)

                       STANDARXT_BLOCK(ri__C, rj__N, R03__C__N,     START3__C__N)
                       STANDARXT_BLOCK(ri__C, rj_CA, R03__C_CA,     START3__C_CA)
                       STANDARXT_BLOCK(ri__C, rj_CB, R03__C_CB, aaj+START3__C_CB)
                       STANDARXT_BLOCK(ri__C, rj__C, R03__C__C,     START3__C__C)
                       STANDARXT_BLOCK(ri__C, rj__O, R03__C__O,     START3__C__O)

                       STANDARXT_BLOCK(ri__O, rj__N, R03__O__N,     START3__O__N)
                       STANDARXT_BLOCK(ri__O, rj_CA, R03__O_CA,     START3__O_CA)
                       STANDARXT_BLOCK(ri__O, rj_CB, R03__O_CB, aaj+START3__O_CB)
                       STANDARXT_BLOCK(ri__O, rj__C, R03__O__C,     START3__O__C)
                       STANDARXT_BLOCK(ri__O, rj__O, R03__O__O,     START3__O__O)

                    } /* if (CA_CA_distance < CA_CA_CUTOFF_SQR) */
                }     /* if (aaj != ZERO_AA) */
            }     /* for all j aa */
            cnt[i] += cnti;
        }      /* if (aai != ZERO_AA) */
    }         /* for all i aa */

    for (i = 0; i < nr_aa; i++) {
        d = (cnt[i] - NR_NEIGHBOUR);
        d = 1.0 - tanh(d);
        n = START4_NEIGHBOUR + aa[i];
        if (n >= NR_PARAM) {
            err_printf (this_sub, "n: %d >= NR_PARAM: %d. Disaster.\n", n, NR_PARAM);
            err_printf (this_sub, "Stopping at %s: %d\n", __FILE__, __LINE__);
            exit (EXIT_FAILURE);
        }
        awght[i][i] -= d*P[n];

/*         awght->scores[i] = -d*P[n]; No contrib from solvation at the moment - pure interaction only */
    }
    free (cnt);
    free(Xf);

    coord_nm_2_a(c);

    return awght;
}


/* ---------mci_contact_rs ------------------------------------
 * writes matrix of tanh contact intensities for a model
 */

int mci_contact_rs( const char *fname, struct coord *coord, const float *P) {
    float **cm;
    const char *this_sub="mci_contact_rs";
    cm = scor_contact_rs(coord, P);
    if (write_mci_mat(fname, coord->size, cm)==EXIT_FAILURE) {
        kill_f_matrix(cm);
        err_printf(this_sub, "failed to write the contact energy map to %s", fname);
        return EXIT_FAILURE;
    }
    kill_f_matrix(cm);
    return EXIT_SUCCESS;
}

/* ---------mci_sc_n_rs ---------------------------------------
 * writes matrix of tanh contact intensities summed with
 * sequence structure fitness for both ends of the interaction.
 */

int 
mci_sc_n_rs( const char *fname, struct scor_set *scrs, struct coord *coord, const float *P) {
    FILE *tfile;
    float **cm;
    const char *this_sub = "mci_sc_n_rs";
    size_t i,j;

    if ((tfile = mfopen(fname, "w", this_sub))==NULL)
        return EXIT_FAILURE;

    cm = scor_contact_rs(coord, P);

    /* non-statistician's joint probability calculation */
    for (i=0;i<coord->size; i++) {
        for (j=i+1; j<coord->size; j++) {
            cm[i][j] += cm[i][i]+cm[j][j];
            cm[j][i] = cm[i][j];
        }
        cm[i][i]=scrs->scores[i];
    }
    if (write_mci_mat(fname, coord->size, cm)==EXIT_FAILURE) {
        kill_f_matrix(cm);
        err_printf(this_sub, "failed to write the contact energy map to %s", fname);
        return EXIT_FAILURE;
    }
    kill_f_matrix(cm);
    return EXIT_SUCCESS;
}
