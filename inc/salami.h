#ifndef SALAMI_H
#define SALAMI_H

/*
 *
 *	salami.h
 *  declarations for salami port (from Thomas Margraf's perl programm)
 * 
 *  @author Steve Hoffmann, shoffmann@zbh.uni-hamburg.de
 *  @company Center for Bioinformatics, Hamburg 
 *  @date 01/14/07 13:22:23 CET  
 *
 */

#include "intsequence.h"
#include "sufmatch.h"

struct salami_sequence {
	const char * sequence;
	unsigned int length;
};

struct salami_info {
	float rmsd, andrew_scr, frac_dme, z_scr, id, sw_cvr, sw_smpl_score, sw_score, sw_score_tot, nw_cvr, nw_smpl_score,
			nw_score, nw_score_tot;
	double tmscore;
	int sw_raw, nw_raw;
	Uint nw_length, sw_length;
};

struct score_mat;
struct coord;
struct pair_set;

void normalize_alt_scores(float *scrs, int len, float *mean, float *dev);
float* get_alt_scores(void *space, int num_scrs, struct score_mat *matrix, struct pair_set *set);
float get_dme_thresh(struct pair_set *set, struct coord *a, struct coord *b);
struct score_struct* get_scores(void *space, struct pair_set *set, struct coord *a, struct coord *b, void *to_use);

struct salami_sequence * salami_sequence_string(void *imbiss, IntSequence *sequence);
void salami_sequence_dump(struct salami_sequence * sequence);
struct salami_info* alignment_aacid(void *imbiss, void *space, const Matchtype *match, IntSequence **s, int len, void *info);
struct salami_info* alignment_wurst(void *config, void *space, const Matchtype *matchtype, IntSequence **sequences, int len, void *info);
#endif
