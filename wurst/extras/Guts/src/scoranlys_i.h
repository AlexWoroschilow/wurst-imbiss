/* functions for analysing local scores for an alignment */

#ifndef SCORANLYS_I_H
#define SCORANLYS_I_H
struct scor_set;
struct pair_set;
struct score_mat;
struct coord;

struct scor_set *scor_set_rs (struct coord *coord, const float *P);
int mci_contact_rs( const char *fname, struct coord *coord, const float *P);
int mci_sc_n_rs( const char *fname, struct scor_set *scrs, struct coord *coord, const float *P);

#endif
