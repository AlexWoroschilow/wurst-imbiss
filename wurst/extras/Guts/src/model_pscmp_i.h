#ifndef MODEL_PSCMP_H
#define MODEL_PSCMP_H
struct coord **
ps_cmp_and_model (struct pair_set *a, struct score_mat *Score_a, struct coord *Sa,
		  struct pair_set *b, struct score_mat *Score_b, struct coord *Sb, 
		  char**frag_string, struct pair_set **cons_pairs);

#endif
