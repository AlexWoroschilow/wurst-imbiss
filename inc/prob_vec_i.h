/*
 * 15 Feb 2005
 * This is the function interface for probability vectors.
 * The structure is defined in prob_vec.h.
 * One may only include this file after <stdlib.h>.
 * rcsid = $Id: prob_vec_i.h,v 1.2 2006/02/22 14:00:47 torda Exp $
 */
#ifndef PROB_VEC_I_H
#define PROB_VEC_I_H

struct prob_vec;
struct prob_vec  *
new_pvec (const size_t frag_len, const size_t prot_len,
          const size_t n_pvec, const size_t n_class);
void              prob_vec_unit_vec (struct prob_vec *p_v);
int               prob_vec_expand (struct prob_vec *p_vec);
void              prob_vec_destroy ( struct prob_vec *p_vec);
char             *prob_vec_info (struct prob_vec *pvec);
size_t            prob_vec_size(const struct prob_vec *pvec);
int               prob_vec_write (struct prob_vec *p_v, const char *fname);
struct prob_vec  *prob_vec_read (const char *fname);

#endif /* PROB_VEC_I_H */
