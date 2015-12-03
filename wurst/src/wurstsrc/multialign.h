/*
 * $Id: multialign.h,v 1.3 2012/04/18 12:02:04 ibondarenko Exp $
 */

#ifndef MULTRIALIGN_H
#define MULTRIALIGN_H

#include "pair_set.h"
#include "score_mat_i.h"

struct prob_vec *
pvec_avg(struct prob_vec *p_vec1,
             struct prob_vec *p_vec2,
             struct pair_set *p_set,
             float weight);

struct pair_set *
merge_alignments(struct pair_set *align1, struct pair_set *align2,
                 struct pair_set *alignment, int c1, int c2);

struct pair_set *
merge_localigns(struct pair_set *align1, struct pair_set *align2, int c1, int c2);

struct pair_set *
merge_align(struct pair_set *align1, struct pair_set *align2,
                 struct pair_set *alignment);

struct pair_set *
remove_seq(struct pair_set *align, int idx);

struct pair_set *
remove_gaps(struct pair_set *align, int idx);

struct pair_set *
split_multal(struct pair_set *align, int a, int b);

#endif /*MULTRIALIGN_H*/
