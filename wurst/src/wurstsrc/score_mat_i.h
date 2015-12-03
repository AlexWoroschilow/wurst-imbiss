/*
 * 1 March 2002
 * $Id: score_mat_i.h,v 1.3 2013/10/16 14:24:09 mosisch Exp $
 */
#ifndef SCORE_MAT_I_H
#define SCORE_MAT_I_H

struct score_mat;
struct seq;
struct pair_set;
struct prob_vec;

void score_mat_info (const struct score_mat *score_mat, float *min, float *max,
                    float *av, float *std_dev);
struct score_mat *score_mat_new (size_t n_rows, size_t n_cols);
void              score_mat_destroy (struct score_mat *s);

struct score_mat *
score_mat_add ( struct score_mat *mat1, struct score_mat *mat2,
                float scale, float shift);

struct score_mat *score_mat_scale (struct score_mat *mat1, const float scale);
struct score_mat *score_mat_shift (struct score_mat *mat1, const float shift);
char *score_mat_string (struct score_mat *smat, struct seq *s1, struct seq *s2);
struct score_mat *score_mat_read (const char *fname);
int score_mat_write (const struct score_mat *smat, const char *fname);
void score_mat_remove_alignment ( struct score_mat *smat, struct pair_set *p_s);

#endif /* SCORE_MAT_I_H */
