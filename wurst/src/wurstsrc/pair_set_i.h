/*
 * 25 September 2001
 * Interface to pair_set operations.
 * This does not show the internals of the structure.
 * rcsid = "$Id: pair_set_i.h,v 1.5 2013/10/16 14:24:09 mosisch Exp $"
 */

#ifndef PAIR_SET_I_H
#define PAIR_SET_I_H

struct pair_set;
struct seq_array;
struct seq;
char  *pair_set_string (struct pair_set *s, struct seq *s1, struct seq *s2);
char *pair_set_stringI(struct pair_set *pairSetString, size_t pos, struct seq *s);
char *pair_set_stringN(struct pair_set *pairSetString, size_t pos);
size_t pair_set_get_startpos(struct pair_set *s, size_t pos);
void  pair_set_dump (struct pair_set *s);
char  *multal_string (struct pair_set *pair_set);
int    pair_set_extend (struct pair_set *s, const size_t n0, const size_t n1,
                        const int long_or_short);
int
pair_set_coverage (struct pair_set *p_s, size_t s1, size_t s2,
                   char **pcover1, char **pcover2);
int
pair_set_gap (struct pair_set *p_s, float *open_cost, float *widen_cost,
              const float open_scale, const float widen_scale);
int    pair_set_score (struct pair_set *s, float *score, float *scr_smpl);
void   pair_set_destroy (struct pair_set *s);
struct pair_set *pair_set_xchange(struct pair_set *p_s);
void   pair_set_get_alignment_indices( struct pair_set *p_s, int sequencenumber, int *start, int *stop  );
int    pair_set_get_alignmentlength( struct pair_set *p_s );
int    pair_set_get_netto_alignmentlength( struct pair_set *p_s );
struct pair_set *pair_set_copy(struct pair_set *ps);
int    pair_set_get_residue( struct pair_set *p_s, int residueindex, int chainindex );

#endif
