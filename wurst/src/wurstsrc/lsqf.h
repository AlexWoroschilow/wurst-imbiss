/* "$Id: lsqf.h,v 1.5 2014/07/04 11:58:48 torda Exp $" */
#ifndef LSQF_H 
#define LSQF_H

struct pair_set;

int get_rmsd(const struct pair_set *pairset, const struct coord *r1,
             const struct coord *r2, float *rmsd_ptr, int *count);

int coord_rmsd (const struct pair_set *pairset, const struct coord *coord1,
                struct coord *coord2, const int sub_flag, float *rmsd,
                struct coord **c1_new, struct coord **c2_new);

double
tm_score (const struct pair_set *pairset, const struct coord *coord1,
          const struct coord *coord2);

double
tm_score_s (const struct pair_set *pairset, const struct coord *coord1,
            const struct coord *coord2);

void coord_centre (struct coord *c);

#endif
