/*
 * Fitting code, originally written by Wilfred van Gunsteren
 */
#ifndef _XOPEN_SOURCE
#   define _XOPEN_SOURCE 500    /* Necessary to get maths constants */
#endif
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "coord_i.h"
#include "coord.h"
#include "pair_set.h"
#include "lsqf.h"
#include "mprintf.h"
#include "seq.h"
#include "read_seq_i.h"

/* Temp debugging lines Please delete me later */
#include <unistd.h>
#include <sys/types.h>
#include <signal.h>
#if !defined (lint) && !defined (DONT_SEE_RCS)
    static const char *rcsid =
    "$Id: lsqf.c,v 1.8 2014/07/04 11:58:08 torda Exp $";
#endif /* !defined (lint) && !defined (DONT_SEE_RCS) */


/* ---------------- constants ---------------------------------
 * We use the square root of two so often, we can help the compiler
 * by saying it is just one value we re-use.
 */
static const double m_sqrt2 = M_SQRT2;

/* ---------------- eigen -------------------------------------
 * This calculates eigenvalues and corresponding eigenvectors.
 * It came from Thomas Huber, who took it from the GROMOS
 * source code...
CCCCC W.F. VAN GUNSTEREN, CAMBRIDGE, JUNE 1979 CCCCCCCCCCCCCCCCCCCCCCCC
                                                                      C
     SUBROUTINE EIGEN (A,R,N,MV)                                      C
                                                                      C
         EIGEN COMPUTES EIGENVALUES AND EIGENVECTORS OF THE REAL      C
     SYMMETRIC N*N MATRIX A, USING THE DIAGONALIZATION METHOD         C
     DESCRIBED IN "MATHEMATICAL METHODS FOR DIGITAL COMPUTERS", EDS.  C
     A.RALSTON AND H.S.WILF, WILEY, NEW YORK, 1962, CHAPTER 7.        C
     IT HAS BEEN COPIED FROM THE IBM SCIENTIFIC SUBROUTINE PACKAGE.   C
                                                                      C
     A(1..N*(N+1)/2) = MATRIX TO BE DIAGONALIZED, STORED IN SYMMETRIC C
                       STORAGE MODE, VIZ. THE I,J-TH ELEMENT (I.GE.J) C
                       IS STORED AT THE LOCATION K=I*(I-1)/2+J IN A;  C
                       THE EIGENVALUES ARE DELIVERED IN DESCENDING    C
                       ORDER ON THE DIAGONAL, VIZ. AT THE LOCATIONS   C
                       K=I*(I+1)/2                                    C
     R(1..N,1..N) = DELIVERED WITH THE CORRESPONDING EIGENVECTORS     C
                    STORED COLUMNWISE                                 C
     N = ORDER OF MATRICES A AND R                                    C
     MV = 0 : EIGENVALUES AND EIGENVECTORS ARE COMPUTED               C
        = 1 : ONLY EIGENVALUES ARE COMPUTED                           C
                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
*/

static void
eigen(float *a, float *r__, int *n, int *mv)
{
    int i__1, i__2, i__3;
    float d__1;

    float cosx, sinx, cosx2, sinx2;
    int i__, j, k, l, m;
    float x, y, anorm, sincs, anrmx;
    int ia, ij, il, im, ll, lm, iq, mm, jq, lq, mq, ind, ilq, imq, ilr,
        imr;
    float thr;
    static const float range = (float) 1e-12;

    /* Parameter adjustments */
    --r__;
    --a;

    if (*mv != 1) {
        iq = -(*n);
        i__1 = *n;
        for (j = 1; j <= i__1; ++j) {
            iq += *n;
            i__2 = *n;
            for (i__ = 1; i__ <= i__2; ++i__) {
                ij = iq + i__;
                r__[ij] = 0.;
                if (i__ - j == 0)
                    r__[ij] = 1.;
            }
        }
    }

    /****COMPUTE INITIAL AND FINAL NORMS (ANORM AND ANRMX) */
    anorm = 0.0;
    i__2 = *n;
    for (i__ = 1; i__ <= i__2; ++i__) {
        i__1 = *n;
        for (j = i__; j <= i__1; ++j) {
            if (i__ - j != 0) {
                ia = i__ + (j * j - j) / 2;
                anorm += a[ia] * a[ia];
            }
        }
    }
    if (anorm > 0.0) {
        anorm = sqrt(anorm) * m_sqrt2;
        anrmx = anorm * range / (float) (*n);

        /*****INITIALIZE INDICATORS AND COMPUTE THRESHOLD, THR */
        thr = anorm;

        do {                       /* ----- while (thr - anrmx > 0.0) ------ */
            thr /= (float) (*n);
            ind = 1 /*EXIT_SUCCESS */ ;
            while (ind) {
                ind = 0 /*EXIT_FAILURE */ ;
                l = 1;
                for (l = 1; l <= *n - 1; l++) {
                    for (m = l + 1; m <= *n; m++) {
                        /*****COMPUT SIN AND COS */
                        mq = (m * m - m) / 2;
                        lq = (l * l - l) / 2;
                        lm = l + mq;
                        d__1 = a[lm];
                        if (fabs(d__1) - thr >= 0.0) {
                            ind = 1 /*EXIT_SUCCESS */ ;
                            ll = l + lq;
                            mm = m + mq;
                            x = (a[ll] - a[mm]) * (float) 0.5;
                            y = -a[lm] / sqrt(a[lm] * a[lm] + x * x);
                            if (x < 0.0)
                                y = -y;
                            sinx =
                                y / sqrt((sqrt(1.0 - y * y) + 1.0) * 2.0);
                            sinx2 = sinx * sinx;
                            cosx = sqrt(1.0 - sinx2);
                            cosx2 = cosx * cosx;
                            sincs = sinx * cosx;

                            /* *****ROTATE L AND M COLUMNS */
                            ilq = *n * (l - 1);
                            imq = *n * (m - 1);
                            i__1 = *n;
                            for (i__ = 1; i__ <= i__1; ++i__) {
                                iq = (i__ * i__ - i__) / 2;
                                if (i__ - l != 0) {
                                    i__2 = i__ - m;
                                    if (i__2 != 0) {
                                        if (i__2 < 0)
                                            im = i__ + mq;
                                        else
                                            im = m + iq;
                                        if (i__ - l >= 0)
                                            il = l + iq;
                                        else
                                            il = i__ + lq;
                                        x = a[il] * cosx - a[im] * sinx;
                                        a[im] =
                                            a[il] * sinx + a[im] * cosx;
                                        a[il] = x;
                                    }     /* ---- (i__2 != 0) ---- */
                                }
                                /* ------ if (i__ - l != 0) ---- */
                                if (*mv != 1) {
                                    ilr = ilq + i__;
                                    imr = imq + i__;
                                    x = r__[ilr] * cosx - r__[imr] * sinx;
                                    r__[imr] =
                                        r__[ilr] * sinx + r__[imr] * cosx;
                                               r__[ilr] = x;
                                }
                            }
                            x = a[lm] * 2. * sincs;
                            y = a[ll] * cosx2 + a[mm] * sinx2 - x;
                            x = a[ll] * sinx2 + a[mm] * cosx2 + x;
                            a[lm] =
                                (a[ll] - a[mm]) * sincs + a[lm] * (cosx2 -
                                                                   sinx2);
                            a[ll] = y;
                            a[mm] = x;
                        }
                        /* *****TESTS FOR COMPLETION */
                        /* *****TEST FOR M = LAST COLUMN */
                    }    /* ---- for (m = l+1; m <= *n; m++)  ---- */
                    /* *****TEST FOR L = SECOND FROM LAST COLUMN */
                }        /* ----- for (l = 1; l < *n-1; l++) ------ */
            }               /* --- while (ind) --- */
        } while (thr > anrmx); /* *****COMPARE THRESHOLD WITH FINAL NORM */
    }
    /* ---- if (anorm > 0) ------ */
    /* *****SORT EIGENVALUES AND EIGENVECTORS */
    iq = -(*n);
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
        iq += *n;
        ll = i__ + (i__ * i__ - i__) / 2;
        jq = *n * (i__ - 2);
        i__2 = *n;
        for (j = i__; j <= i__2; ++j) {
            jq += *n;
            mm = j + (j * j - j) / 2;
            if (a[ll] - a[mm] < 0.0) {
                x = a[ll];
                a[ll] = a[mm];
                a[mm] = x;
                if (*mv != 1) {
                    i__3 = *n;
                    for (k = 1; k <= i__3; ++k) {
                        ilr = iq + k;
                        imr = jq + k;
                        x = r__[ilr];
                        r__[ilr] = r__[imr];
                        r__[imr] = x;
                    }
                }
            }
        }
    }
}

/* ---------------- apply_rot  ----------------------------------
 * apply rotation matrix rmat to structure of size nr_atoms.
 */
static void
apply_rot(float rmat[3][3], struct coord *structure)
{
    size_t n;
    struct RPoint *r1 = structure->rp_c;
    struct RPoint dr;
    for (n = 0; n < structure->size; n++) {
        dr.x =
            rmat[0][0] * r1[n].x + rmat[0][1] * r1[n].y +
            rmat[0][2] * r1[n].z;
        dr.y =
            rmat[1][0] * r1[n].x + rmat[1][1] * r1[n].y +
            rmat[1][2] * r1[n].z;
        dr.z =
            rmat[2][0] * r1[n].x + rmat[2][1] * r1[n].y +
            rmat[2][2] * r1[n].z;
        r1[n].x = dr.x;
        r1[n].y = dr.y;
        r1[n].z = dr.z;
    }
    r1 = structure->rp_ca;
    for (n = 0; n < structure->size; n++) {
        dr.x =
            rmat[0][0] * r1[n].x + rmat[0][1] * r1[n].y +
            rmat[0][2] * r1[n].z;
        dr.y =
            rmat[1][0] * r1[n].x + rmat[1][1] * r1[n].y +
            rmat[1][2] * r1[n].z;
        dr.z =
            rmat[2][0] * r1[n].x + rmat[2][1] * r1[n].y +
            rmat[2][2] * r1[n].z;
        r1[n].x = dr.x;
        r1[n].y = dr.y;
        r1[n].z = dr.z;
    }
    r1 = structure->rp_cb;
    for (n = 0; n < structure->size; n++) {
        dr.x =
            rmat[0][0] * r1[n].x + rmat[0][1] * r1[n].y +
            rmat[0][2] * r1[n].z;
        dr.y =
            rmat[1][0] * r1[n].x + rmat[1][1] * r1[n].y +
            rmat[1][2] * r1[n].z;
        dr.z =
            rmat[2][0] * r1[n].x + rmat[2][1] * r1[n].y +
            rmat[2][2] * r1[n].z;
        r1[n].x = dr.x;
        r1[n].y = dr.y;
        r1[n].z = dr.z;
    }
    r1 = structure->rp_n;
    for (n = 0; n < structure->size; n++) {
        dr.x =
            rmat[0][0] * r1[n].x + rmat[0][1] * r1[n].y +
            rmat[0][2] * r1[n].z;
        dr.y =
            rmat[1][0] * r1[n].x + rmat[1][1] * r1[n].y +
            rmat[1][2] * r1[n].z;
        dr.z =
            rmat[2][0] * r1[n].x + rmat[2][1] * r1[n].y +
            rmat[2][2] * r1[n].z;
        r1[n].x = dr.x;
        r1[n].y = dr.y;
        r1[n].z = dr.z;
    }
    r1 = structure->rp_o;
    for (n = 0; n < structure->size; n++) {
        dr.x =
            rmat[0][0] * r1[n].x + rmat[0][1] * r1[n].y +
            rmat[0][2] * r1[n].z;
        dr.y =
            rmat[1][0] * r1[n].x + rmat[1][1] * r1[n].y +
            rmat[1][2] * r1[n].z;
        dr.z =
            rmat[2][0] * r1[n].x + rmat[2][1] * r1[n].y +
            rmat[2][2] * r1[n].z;
        r1[n].x = dr.x;
        r1[n].y = dr.y;
        r1[n].z = dr.z;
    }
    return;
}

/* ---------------- lsq_fit  ----------------------------------
 * Least square fit routine to determine the rotation matrix to
 * superimpose coordinates r1 onto coordinates r2.
 * omega is a symmetric matrix in symmetric storage mode:
 * The element [i][j] is stored at position [i*(i+1)/2+j] with i>j
 * Thanks to Wilfred and Thomas Huber.
 */
static int
lsq_fit(const int nr_atoms, const struct RPoint *r1, const struct RPoint *r2,
        float R[3][3])
{

    int i, j, ii, jj, n;
    float U[3][3], det_U, sign_detU, sigma;
    float H[3][3], K[3][3];     /*, R[3][3]; */
    float omega[21], eve_omega[36], eva_omega[6];
    static const float TINY  = (float) 1.e-10;
    static const float SMALL = (float) 1.e-5;
    static const char *this_sub = "lsq_fit";

/* ----- CALCULATE THE MATRIX U AND ITS DETERMINANT ----- */
    for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++)
            U[i][j] = 0.0;

    for (n = 0; n < nr_atoms; n++) {
        U[0][0] += r1[n].x * r2[n].x;
        U[0][1] += r1[n].x * r2[n].y;
        U[0][2] += r1[n].x * r2[n].z;
        U[1][0] += r1[n].y * r2[n].x;
        U[1][1] += r1[n].y * r2[n].y;
        U[1][2] += r1[n].y * r2[n].z;
        U[2][0] += r1[n].z * r2[n].x;
        U[2][1] += r1[n].z * r2[n].y;
        U[2][2] += r1[n].z * r2[n].z;
    }

    det_U = U[0][0] * U[1][1] * U[2][2] + U[0][2] * U[1][0] * U[2][1] +
        U[0][1] * U[1][2] * U[2][0] - U[2][0] * U[1][1] * U[0][2] -
        U[2][2] * U[1][0] * U[0][1] - U[2][1] * U[1][2] * U[0][0];

    if (fabs(det_U) < TINY) {
        err_printf(this_sub, "determinant of U equal to zero\n");
        return EXIT_FAILURE;
    }

    sign_detU = det_U / fabs(det_U);    /* sign !!! */

/* ----- CONSTRUCT OMEGA, DIAGONALIZE IT AND DETERMINE H AND K --- */

    for (i = 0; i < 6; i++)
        for (j = i; j < 6; j++)
            omega[(j * (j + 1) / 2) + i] = 0.0;
    for (j = 3; j < 6; j++) {
        jj = j * (j + 1) / 2;
        for (i = 0; i < 3; i++) {
            ii = jj + i;
            omega[ii] = U[i][j - 3];
        }
    }

    i = 6;                      /* dimension of omega matrix */
    j = 0;                      /* both, eigenvalues and eigenvectors are calculated */
    eigen(omega, eve_omega, &i, &j);

    for (i = 0; i < 6; i++)
        eva_omega[i] = omega[i * (i + 1) / 2 + i];

    if (det_U < 0.0) {
        if (fabs(eva_omega[1] - eva_omega[2]) < SMALL) {
            err_printf(this_sub,
                       "determinant of U < 0 && degenerated eigenvalues\n");
            return EXIT_FAILURE;
        }
    }

    for (i = 0; i < 3; i++){
        for (j = 0; j < 3; j++) {
            H[i][j] = m_sqrt2 * eve_omega[j * 6 + i];
            K[i][j] = m_sqrt2 * eve_omega[j * 6 + i + 3];
        }
    }
    sigma = (H[1][0] * H[2][1] - H[2][0] * H[1][1]) * H[0][2] +
        (H[2][0] * H[0][1] - H[0][0] * H[2][1]) * H[1][2] +
        (H[0][0] * H[1][1] - H[1][0] * H[0][1]) * H[2][2];

    if (sigma <= 0.0) {
        for (i = 0; i < 3; i++) {
            H[i][2] = -H[i][2];
            K[i][2] = -K[i][2];
        }
    }

/* --------- DETERMINE R AND ROTATE X ----------- */
    for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++)
            R[j][i] = K[j][0] * H[i][0] + K[j][1] * H[i][1] +
                sign_detU * K[j][2] * H[i][2];


#ifdef DEBUG
    mfprintf(stdout, "Rotation matrix:\n");
    for (j = 0; j < 3; j++) {
        for (i = 0; i < 3; i++)
            mfprintf(stdout, "%f ", R[i][j]);
        mfprintf(stdout, "\n");
    }
#endif

    return EXIT_SUCCESS;
}

/* ---------------- calc_RMSD ------------------------------------
 * Calculates the RMSD between nr_atoms Points in arrays r1 and r2
 */
static float
calc_RMSD(const int nr_atoms, const struct RPoint *r1, const struct RPoint *r2)
{
    float rmsd = 0.0;
    float dr_sqrlength = 0.0;
    struct RPoint dr;
    int n;

    for (n = 0; n < nr_atoms; n++) {
        dr.x = r1[n].x - r2[n].x;
        dr.y = r1[n].y - r2[n].y;
        dr.z = r1[n].z - r2[n].z;
        dr_sqrlength = dr.x * dr.x + dr.y * dr.y + dr.z * dr.z;
        rmsd += dr_sqrlength;
    }
    rmsd /= nr_atoms;
    rmsd = sqrt(rmsd);

    return (rmsd);
}

/* ---------------- CM_Translate ------------------------------
 * Translates the coordinates of two molecules to superimpose such
 * that the center of mass is in the origin.
 */
static struct RPoint
calc_CM(const size_t nr_atoms, const struct RPoint *r1)
{
    unsigned int i;
    float total_mass1 = nr_atoms;
    struct RPoint cm1;

    cm1.x = 0.0;
    cm1.y = 0.0;
    cm1.z = 0.0;

    for (i = 0; i < nr_atoms; i++) {
        cm1.x += r1[i].x;
        cm1.y += r1[i].y;
        cm1.z += r1[i].z;
    }
    cm1.x /= total_mass1;
    cm1.y /= total_mass1;
    cm1.z /= total_mass1;

    return cm1;
}

/* ---------------- apply_trans -------------------------------
 * Apply the transformation specified by the vector trans to structure
 * of size nr_atoms
 */
static void
apply_trans(const struct RPoint *trans, struct RPoint *structure,
            size_t nr_atoms)
{
    struct RPoint *r1;
    size_t i;
    const struct RPoint cm1 = *trans;
    r1 = structure;
    for (i = 0; i < nr_atoms; i++) {
        r1[i].x -= cm1.x;
        r1[i].y -= cm1.y;
        r1[i].z -= cm1.z;
    }
    return;
}

/* ---------------- coord_trans -------------------------------
 * apply the transformation specified by the vector trans to structure
 * of length size. Uses apply_trans().
 */
static void
coord_trans(const struct RPoint *trans, struct coord *structure,
            const size_t size)
{
    apply_trans(trans, structure->rp_ca, size);
    apply_trans(trans, structure->rp_cb, size);
    apply_trans(trans, structure->rp_n, size);
    apply_trans(trans, structure->rp_c, size);
    apply_trans(trans, structure->rp_o, size);
}

/* ---------------- coord_centre -------------------------------
 * transform all atom coordinates so that the centre of mass is
 * in the origin
 */
void
coord_centre(struct coord *c)
{
    struct RPoint t = calc_CM(coord_size(c), c->rp_ca);
    coord_trans(&t, c, coord_size(c));
}


/*-------------------- copy_coord_elem   ----------------
 * We have an atom / site from one coord structure, to
 * be copied into a second.
 */
static void
copy_coord_elem(struct coord *dst, const struct coord *src, const int dst_ndx,
                const int src_ndx)
{
    /*struct RPoint spoint;
    struct RPoint dpoint;
    spoint = src->rp_ca[src_ndx];
    dpoint = dst->rp_ca[dst_ndx];
    dpoint.x = spoint.x;
    dpoint.y = spoint.y;
    dpoint.z = spoint.z;
    dst->rp_ca[dst_ndx] = dpoint;*/
    dst->rp_ca[dst_ndx].x = src->rp_ca[src_ndx].x;
    dst->rp_ca[dst_ndx].y = src->rp_ca[src_ndx].y;
    dst->rp_ca[dst_ndx].z = src->rp_ca[src_ndx].z;
    dst->rp_cb[dst_ndx].x = src->rp_cb[src_ndx].x;
    dst->rp_cb[dst_ndx].y = src->rp_cb[src_ndx].y;
    dst->rp_cb[dst_ndx].z = src->rp_cb[src_ndx].z;
    dst->rp_n[dst_ndx].x = src->rp_n[src_ndx].x;
    dst->rp_n[dst_ndx].y = src->rp_n[src_ndx].y;
    dst->rp_n[dst_ndx].z = src->rp_n[src_ndx].z;
    dst->rp_c[dst_ndx].x = src->rp_c[src_ndx].x;
    dst->rp_c[dst_ndx].y = src->rp_c[src_ndx].y;
    dst->rp_c[dst_ndx].z = src->rp_c[src_ndx].z;
    dst->rp_o[dst_ndx].x = src->rp_o[src_ndx].x;
    dst->rp_o[dst_ndx].y = src->rp_o[src_ndx].y;
    dst->rp_o[dst_ndx].z = src->rp_o[src_ndx].z;
    dst->orig[dst_ndx] = src->orig[src_ndx];
    dst->icode[dst_ndx] = src->icode[src_ndx];
    if (src->psi)
        dst->psi[dst_ndx] = src->psi[src_ndx];
    if (src->phi)
        dst->phi[dst_ndx] = src->phi[src_ndx];
    if (src->sec_typ)
        dst->sec_typ[dst_ndx] = src->sec_typ[src_ndx];
    dst->seq->seq[dst_ndx] = src->seq->seq[src_ndx];
}

/*-------------------- coord_rmsd ------------------------------------------
 * Takes a pairset and two structures and moves coord1 on top of coord2 so
 * that the RMSD is minimized. The superimposed structures are returned as
 * c1_new and c2_new. If sub_flag is set, the returned structures contain only
 * the subset of residues specified by the pairset.
 */
int
coord_rmsd(const struct pair_set *const pairset, const struct coord *source,
           struct coord *target, const int sub_flag, float *rmsd,
           struct coord **c1_new, struct coord **c2_new)
{
    float rmat[3][3];
    size_t tmp_idx = 0;
    size_t size, i;
    int r, a, b;
    struct RPoint translation1, translation2;
    int **pairs = pairset->indices;
    struct coord *temp_struct_1 =
        coord_template(source, coord_size(source));
    struct coord *temp_struct_2 =
        coord_template(target, coord_size(target));
    const char *this_sub = "coord_rmsd";
    /* Create temporary structures to hold the aligned parts of the two structures */
    temp_struct_1->seq = seq_copy(source->seq);
    temp_struct_2->seq = seq_copy(target->seq);
    coord_nm_2_a(temp_struct_1);
    coord_nm_2_a(temp_struct_2);
    if(sub_flag > 4){
        for (i = 0; i < pairset->n; i++) {
            a = pairs[i][0];
            b = pairs[i][1];
            if (a != GAP_INDEX && b != GAP_INDEX) {
                copy_coord_elem(temp_struct_1, source, tmp_idx, pairs[i][1]);
                copy_coord_elem(temp_struct_2, target, tmp_idx, pairs[i][0]);
                tmp_idx++;
            }
        }
    } else {
        for (i = 0; i < pairset->n; i++) {
            a = pairs[i][0];
            b = pairs[i][1];
            if (a != GAP_INDEX && b != GAP_INDEX) {
                copy_coord_elem(temp_struct_1, source, tmp_idx, pairs[i][0]);
                copy_coord_elem(temp_struct_2, target, tmp_idx, pairs[i][1]);
                tmp_idx++;
            }
        }
    }
    size = tmp_idx;
    coord_trim(temp_struct_1, size);
    coord_trim(temp_struct_2, size);
    /* Determine the center of mass of the two structures
     * and move the CoMs to the coordinate origin.
     */
    translation1 = calc_CM(size, temp_struct_1->rp_ca);
    translation2 = calc_CM(size, temp_struct_2->rp_ca);
    coord_trans(&translation1, temp_struct_1, size);
    coord_trans(&translation2, temp_struct_2, size);

    /* Determine the rotation matrix and apply the
     * rotation to structure 2. Then calculate the RMSD
     */

    r = lsq_fit(tmp_idx, temp_struct_1->rp_ca, temp_struct_2->rp_ca,
                rmat);
    if (r == EXIT_FAILURE) {
        err_printf(this_sub, "lsq_fit fail\n");
        *rmsd = -1;
        goto escape;
    }
    apply_rot(rmat, temp_struct_1);
    *rmsd = calc_RMSD(size, temp_struct_1->rp_ca, temp_struct_2->rp_ca);

    /* Move the structures' CoMs back to the original CoM of structure 2.*/
    translation2.x *= -1;
    translation2.y *= -1;
    translation2.z *= -1;
    /* If only the aligned subset of the structures is needed,
     * translate and return the temporary structures.
     */
    if (sub_flag % 2) {
        coord_trans(&translation2, temp_struct_1, size);
        coord_trans(&translation2, temp_struct_2, size);
        *c1_new = temp_struct_1;
        *c2_new = temp_struct_2;
    }
    /* Otherwise create a copy of the original structures, apply translation
     * and rotation to the copies and return those.
     */
    else {
        coord_destroy(temp_struct_1);
        coord_destroy(temp_struct_2);
        *c1_new = coord_template(source, coord_size(source));
        *c2_new = coord_template(target, coord_size(target));
        (*c1_new)->seq = seq_copy(source->seq);
        (*c2_new)->seq = seq_copy(target->seq);
        for (i = 0; i < coord_size(source); i++)
            copy_coord_elem(*c1_new, source, i, i);
        for (i = 0; i < coord_size(target); i++)
            copy_coord_elem(*c2_new, target, i, i);
        coord_trans(&translation1, *c1_new, (*c1_new)->size);
        apply_rot(rmat, *c1_new);
        coord_trans(&translation2, *c1_new, (*c1_new)->size);
    }

    return (EXIT_SUCCESS);
  escape:
    return (EXIT_FAILURE);
}


/* ---------------- get_rmsd      -----------------------------
 * We should change the return type to unsigned int and check
 * if count can be a pointer to an unsigned in.
 */
int
get_rmsd(const struct pair_set *pairset, const struct coord *r1,
    const struct coord *r2, float *rmsd_ptr, int *count){
    float dr_sqrlength = 0.0;
    float rmsd = 0.0;
    struct RPoint dr;
    size_t n;

    unsigned pairs = 0;
    for (n = 0; n < pairset->n; n++) {
        if(pairset->indices[n][0] != GAP_INDEX && pairset->indices[n][1] != GAP_INDEX){
            dr.x = r1->rp_ca[pairset->indices[n][0]].x - r2->rp_ca[pairset->indices[n][1]].x;
            dr.y = r1->rp_ca[pairset->indices[n][0]].y - r2->rp_ca[pairset->indices[n][1]].y;
            dr.z = r1->rp_ca[pairset->indices[n][0]].z - r2->rp_ca[pairset->indices[n][1]].z;
            dr_sqrlength = dr.x * dr.x + dr.y * dr.y + dr.z * dr.z;
            rmsd += sqrt(dr_sqrlength);
            pairs++;
        }
    }

    *rmsd_ptr = rmsd;
    *count = pairs;
    return pairs;
}

/* ---------------- constants ---------------------------------
 * There are some magic numbers used in the calculation of the
 * TM score, Zhang, Y, Skolnick, J. Proteins, 57, 702-710 (2004).
 * Define them once since they are used twice below.
 */
static const float one24   = (float) 1.24;
static const float fifteen = (float) 15.;
static const float one8    = (float) 1.8;

/* ---------------- tm_score ----------------------------------
 */
double
tm_score ( const struct pair_set *pairset, const struct coord *r1,
           const struct coord *r2 )
{
    double tmscore = 0.0;
    double dO = 0.0;
    double di = 0.0;
    size_t ln = 0;
    size_t lt = 0;
    struct RPoint dr;
    unsigned int n = 0;
    double dr_sqrlength = 0.0;
    if(!pairset || !r1 || !r2)
        return(0.0);

    ln = ( r1->seq->length > r2->seq->length ) ? r1->seq->length : r2->seq->length;

    for (n = 0; n < pairset->n; n++)
        if(pairset->indices[n][0] != GAP_INDEX && pairset->indices[n][1] != GAP_INDEX)
            lt++;

    dO = (one24 * cbrt(ln - fifteen)) - one8;  /* Constants from Skolnick paper */

    for (n = 0; n < pairset->n; n++) {
        if(pairset->indices[n][0] != GAP_INDEX && pairset->indices[n][1] != GAP_INDEX){
            dr.x = r1->rp_ca[pairset->indices[n][0]].x - r2->rp_ca[pairset->indices[n][1]].x;
            dr.y = r1->rp_ca[pairset->indices[n][0]].y - r2->rp_ca[pairset->indices[n][1]].y;
            dr.z = r1->rp_ca[pairset->indices[n][0]].z - r2->rp_ca[pairset->indices[n][1]].z;
            dr_sqrlength = dr.x * dr.x + dr.y * dr.y + dr.z * dr.z;
            di = sqrt(dr_sqrlength);
            tmscore += 1 / ( 1 + (di / dO) * (di / dO));
        }
    }

    tmscore = tmscore / ln;
    return(tmscore);
}

/* ---------------- tm_score_s    -----------------------------
 */
double
tm_score_s ( const struct pair_set *pairset, const struct coord *r1,
             const struct coord *r2 )
{
    double tmscore = 0.0;
    double dO = 0.0;
    double di = 0.0;
    size_t ln = 0;
    size_t lt = 0;
    struct RPoint dr;
    unsigned int n = 0;
    double dr_sqrlength = 0.0;
    if(!pairset || !r1 || !r2)
        return(0.0);

    ln = ( r1->seq->length > r2->seq->length ) ? r2->seq->length : r1->seq->length;

    for (n = 0; n < pairset->n; n++)
        if(pairset->indices[n][0] != GAP_INDEX && pairset->indices[n][1] != GAP_INDEX)
            lt++;

    dO = (one24 * cbrt(ln - fifteen)) - one8;  /* Constants from Skolnick paper */
    for (n = 0; n < pairset->n; n++) {
        if(pairset->indices[n][0] != GAP_INDEX && pairset->indices[n][1] != GAP_INDEX){
            dr.x = r1->rp_ca[pairset->indices[n][0]].x - r2->rp_ca[pairset->indices[n][1]].x;
            dr.y = r1->rp_ca[pairset->indices[n][0]].y - r2->rp_ca[pairset->indices[n][1]].y;
            dr.z = r1->rp_ca[pairset->indices[n][0]].z - r2->rp_ca[pairset->indices[n][1]].z;
            dr_sqrlength = dr.x * dr.x + dr.y * dr.y + dr.z * dr.z;
            di = sqrt(dr_sqrlength);
            tmscore += 1 / ( 1 + (di / dO) * (di / dO));
        }
    }
    
    tmscore = tmscore / ln;
    return(tmscore);
}
