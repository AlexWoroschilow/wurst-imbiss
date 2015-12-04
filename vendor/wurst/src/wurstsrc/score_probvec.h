/*
 * 14 June 2005
 * rcsid = $Id: score_probvec.h,v 1.3 2013/10/16 14:24:09 mosisch Exp $
 */

#ifndef SCORE_PROBVEC_H
#define SCORE_PROBVEC_H

struct score_mat;
struct prob_vec;
int
score_pvec (struct score_mat *score_mat,
            struct prob_vec *p_v1, struct prob_vec *p_v2);
#endif /* SCORE_PROBVEC_H */
