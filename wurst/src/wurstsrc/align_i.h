/*
 * 27 Aug 2001
 * rcsid = $Id: align_i.h,v 1.5 2014/07/04 12:37:53 torda Exp $;
 */
#ifndef ALIGN2_H
#define ALIGN2_H

struct score_mat;
struct coord;
struct pair_set *
score_mat_sum_full ( struct score_mat **rmat, struct score_mat *smat,
                     float pgap_open, float pgap_widen,
                     float qgap_open, float qgap_widen,
                     float *p_mult,   float *q_mult,
                     const int algn_type, const struct pair_set *bias_set);
  
size_t
get_pair_set_m (const struct pair_set *ps);

size_t
get_pair_set_n (const struct pair_set *ps);

char *
pair_set_get_strNum (const struct pair_set *ps, const size_t strNum,
                     const struct coord *c);

#endif  /* ALIGN2_H */
