/*
 * 6 Sep 2011
 * Interface for reading autoclass sequence classification for Calpha and
 * Tau angles.
 */
#ifndef READ_AC_CALPHA_I_H
#define READ_AC_CALPHA_I_H

struct calpha_clssfcn;
struct prob_vec;
struct coord;
struct seq;

struct calpha_clssfcn *
ac_read_calpha (const char *fname, float t_error, float c_error, int corr_num);

float
calpha_clssfcn_return (struct calpha_clssfcn *ca_clssfcn,
                         int i_class, int pos, int att_num);

void
calpha_clssfcn_destroy (struct calpha_clssfcn *ca_clssfcn);

size_t
att_size_calpha (const struct calpha_clssfcn *ca_clssfcn);

size_t
ac_nclass_calpha (const struct calpha_clssfcn *ca_clssfcn);


struct prob_vec *
calpha_strct_2_prob_vec (struct coord *structure,
                         const struct calpha_clssfcn *ca,
                         const int norm);

void
calpha_getFragment(const size_t residue_num, const size_t n_att,
                   const size_t cor_n_att, struct coord *structure,
                   float *frag);

#endif /* READ_AC_CALPHA_I_H */
