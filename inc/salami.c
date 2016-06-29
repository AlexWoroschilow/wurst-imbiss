/*
 *  salami.c
 *  a partial prot from Thomas Margraf's salami perl script
 *  to allmighty c w/ some adjustments
 * 
 *  @author Steve Hoffmann
 *  @email shoffmann@zbh.uni-hamburg.de
 *  @date 01/14/07 13:05:10 CET
 *  
 */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <getopt.h>
#include <sys/types.h>
#include <sys/times.h>
#include <time.h>
#include <string.h>
#include "memory.h"
#include "intsequence.h"
#include "sufmatch.h"
#include "depictseqs.h"
#include "stringutils.h"
#include "mathematics.h"
#include "salami.h"
#include "wurstimbiss.h"

/*WURST include*/
#include "read_seq_i.h"
#include "score_mat_i.h"
#include "prob_vec_i.h"
#include "prob_vec.h"
#include "coord_i.h"
#include "pair_set.h"
#include "pair_set_i.h"
#include "score_probvec.h"
#include "score_smat.h"
#include "coord.h"
#include "align_i.h"
#include "matrix.h"
#include "model.h"
#include "cmp_dmat_i.h"
#include "altscores.h"
#include "lsqf.h"
#include "seq.h"
#include "massert.h"
#include "read_seq_i.h"
#include "pair_set_p_i.h"

struct score_struct {
	float scr_tot;
	float cvr;
};

const float gapopen = 3.25;
const float gapwide = 0.8942;
const float sw1_pgap_open = 3.25;
const float sw1_pgap_widen = 0.8942;
const float nw1_pgap_open = 3.25;
const float nw1_pgap_widen = 0.8942;
const int N_AND_W = 0;
const int S_AND_W = 1;

/*--------------------------- normalize_alt_scores ---------------------------
 *    
 * a port from salami's normalize_alt_scores function.
 * 
 */

void normalize_alt_scores(float *scrs, int len, float *mean, float *dev) {
	int i;
	float temp;

	*mean = 0.0;
	*dev = 0.0;

	for (i = 0; i < len; i++) {
		*mean += scrs[i];
	}

	(*mean) /= len;

	for (i = 0; i < len; i++) {
		temp = scrs[i] - (*mean);
		*dev += temp * temp;
	}

	*dev /= (len - 1);
	*dev = sqrt((*dev));

	return;
}

/*------------------------------ get_alt_scores ------------------------------
 *    
 * a port from salami's get_alt_* function. Calculates scores on random paths 
 * through the scoring matrix parameters: number of paths/scores, scoring 
 * matrix, pair set of optimal path. 
 * 
 * @return scores
 */

float*
get_alt_scores(void *space, int num_scrs, struct score_mat *matrix, struct pair_set *set) {
	int i;
	float *scrs;

	scrs = ALLOCMEMORY(space, NULL, float, num_scrs);

	for (i = 0; i < num_scrs; i++) {
		scrs[i] = find_alt_path_score_simple(matrix, set);
	}

	return scrs;
}

/*------------------------------ get_dme_thresh ------------------------------
 *    
 * a port from salami's get_dme_thresh function
 * 
 */

float get_dme_thresh(struct pair_set *set, struct coord *a, struct coord *b) {
	struct coord *model;
	struct seq *seq_a;
	float frac = 0.0;

	seq_a = coord_get_seq(a);
	model = make_model(set, seq_a, b);

	if (model == NULL) {
		return 0.0;
	}

	if (coord_size(model) < 10) {
		return 0.0;
	}

	dme_thresh(&frac, a, model, 3.0);

	coord_destroy(model);
	seq_destroy(seq_a);

	return frac;
}

/*-------------------------------- get_scores --------------------------------
 *    
 * a port from salami's get_score function
 * 
 */

struct score_struct*
get_scores(void *space, struct pair_set *set, struct coord *a, struct coord *b, void *to_use) {
	int cover = 0, i;
	float cover_float;
	float smpl_score, scr_tot;
	char *pcover1;
	char *pcover2;
	struct seq* seq_a;
	struct score_struct *scores;
	/*float geo_gap = 0, nseq_gap = 0;*/
	float open_cost = 0, widen_cost = 0;

	scores = malloc(sizeof(struct score_struct));
	smpl_score = set->smpl_score;

	pair_set_coverage(set, coord_size(a), coord_size(b), &pcover1, &pcover2);
	for (i = 0; i < coord_size(a); i++) {
		if (pcover1[i] == '1')
			cover++;
	}

	free(pcover1);
	free(pcover2);

	seq_a = coord_get_seq(a);
	cover_float = (float) cover / (float) seq_size(seq_a);

	if (cover_float >= 0.05) {
		pair_set_gap(set, &open_cost, &widen_cost, 1, 1);
	}

	scr_tot = smpl_score + sw1_pgap_open * open_cost + sw1_pgap_widen * widen_cost;

	seq_destroy(seq_a);
	scores->scr_tot = scr_tot;
	scores->cvr = cover_float;

	return scores;
}

struct salami_sequence * salami_sequence_string(void *imbiss, IntSequence *sequence) {

	imbissinfo *imbissinfo = imbiss;
	massert((imbiss != NULL), "Imbiss info object can not be empty");

	char *binary = merge(merge(merge(imbissinfo->binarypath, "/"), sequence_code(sequence->url)), ".bin");

	struct coord *coordinates = coord_read(binary);
	massert((coordinates != NULL), "Coordinates object can not be null");
	free(binary);

	struct seq *sequence_wurst = coord_get_seq(coordinates);
	massert((sequence_wurst != NULL), "Coordinates object can not be null");

	struct salami_sequence * response = malloc(sizeof(struct salami_sequence));
	massert((response != NULL), "Salami sequence object can not be null");

	response->length = sequence_wurst->length;
	response->sequence = seq_print(sequence_wurst);

	coord_destroy(coordinates);
	seq_destroy(sequence_wurst);

	return response;
}

void salami_sequence_dump(struct salami_sequence * sequence) {
	unsigned long i;
	printf("Sequence length: %u \n", sequence->length);
	for (i = 0; i < sequence->length; i++) {
		printf("%d.", sequence->sequence[i]);
	}
	printf("\n");
}

struct salami_info* alignment_wurst(void *config, void *space, const Matchtype *match, IntSequence **sequences, int len, void *info) {
	float frac_dme, z_scr;
	float *altscores;
	unsigned int id;
	float mean, dev, rmsd, andrew_scr;
	double tmscore;
	stringset_t *in;
	struct coord *coord_a, *coord_b, *rmsd_a = NULL, *rmsd_b = NULL;
	struct seq *seq_a, *seq_b;
	struct prob_vec *pvec_a, *pvec_b;
	struct score_mat *matrix;
	struct score_struct *scores;
	struct pair_set *set_sw, *set_nw;
	struct score_mat *crap = NULL;
	struct salami_info *salami = NULL;

	imbissinfo *imbiss = (imbissinfo*) config;
	massert((imbiss != NULL), "Imbiss info object can not be empty");

	salami = malloc(sizeof(struct salami_info));
	massert((salami != NULL), "Can not allocate memory for salami object");

	in = (stringset_t*) info;

	char *binary = merge(merge(merge(imbiss->binarypath, "/"), sequence_code(sequences[match->id]->url)), ".bin");
	char *vector = merge(merge(merge(imbiss->vectorpath, "/"), sequence_code(sequences[match->id]->url)), ".vec");

	coord_a = coord_read(binary);
	coord_b = coord_read(in->strings[0].str);

	massert((coord_a != NULL), "Coordinates for sequence A can not be empty");
	massert((coord_b != NULL), "Coordinates for sequence B can not be empty");

	seq_a = coord_get_seq(coord_a);
	seq_b = coord_get_seq(coord_b);

	massert((seq_a != NULL), "Sequence A can not be empty");
	massert((seq_b != NULL), "Sequence B can not be empty");

	pvec_a = prob_vec_read(vector);
	pvec_b = prob_vec_read(in->strings[1].str);

	massert((pvec_a != NULL), "Vector for sequence A can not be empty");
	massert((pvec_b != NULL), "Vector for sequence B can not be empty");

	matrix = score_mat_new(seq_size(seq_a), seq_size(seq_b));
	score_pvec(matrix, pvec_a, pvec_b);

	prob_vec_destroy(pvec_a);
	prob_vec_destroy(pvec_b);

	FREEMEMORY(space, binary);
	FREEMEMORY(space, vector);

	/*nw align*/

	set_nw = score_mat_sum_full(&crap, matrix, sw1_pgap_open, sw1_pgap_widen, sw1_pgap_open, sw1_pgap_widen,
	NULL, NULL, N_AND_W, NULL);

	id = get_seq_id_simple(set_nw, seq_a, seq_b);
	scores = get_scores(space, set_nw, coord_a, coord_b, NULL);

	salami->id = (float)id / (float)set_nw->n;
	salami->nw_score = set_nw->score;
	salami->nw_length = (Uint)set_nw->n;
	salami->nw_smpl_score = set_nw->smpl_score;
	salami->nw_score_tot = scores->scr_tot;
	salami->nw_cvr = scores->cvr;
	salami->nw_raw = scores->cvr * seq_size(seq_a);

	score_mat_destroy(crap);
	pair_set_destroy(set_nw);
	free(scores);

	/*sw align*/

	set_sw = score_mat_sum_full(&crap, matrix, sw1_pgap_open, sw1_pgap_widen, sw1_pgap_open, sw1_pgap_widen,
	NULL, NULL, S_AND_W, NULL);
	score_mat_destroy(crap);
	scores = get_scores(space, set_sw, coord_a, coord_b, NULL);

	salami->sw_score = set_sw->score;
	salami->sw_length = (Uint)set_sw->n;
	salami->sw_smpl_score = set_sw->smpl_score;
	salami->sw_score_tot = scores->scr_tot;
	salami->sw_cvr = scores->cvr;
	salami->sw_raw = scores->cvr * seq_size(seq_a);

	frac_dme = get_dme_thresh(set_sw, coord_a, coord_b);
	altscores = get_alt_scores(space, 1000, matrix, set_sw);
	normalize_alt_scores(altscores, 1000, &mean, &dev);

	/*@todo: get scores from get_scores and calc z-scr here!*/
	z_scr = 0.0;
	if (dev != 0) {
		z_scr = (scores->scr_tot - mean) / dev;
	}

	/*rmsd*/
	coord_rmsd(set_sw, coord_a, coord_b, 0, &rmsd, &rmsd_a, &rmsd_b);
	tmscore = tm_score(set_sw, rmsd_a, rmsd_b);

	/*andrew score*/
	/*frac_dme (9) cvr (3)*/
	/*alignmentsize (6) = sw_cvr (3) * seq_size(coord_get_seq(coord_b))*/
	/*andrew_scr (2)= frac_dme(9) * (alignmentsize(6) / minsize) */

	andrew_scr = frac_dme * ((scores->cvr * seq_size(seq_a)) / (MIN(seq_size(seq_a), seq_size(seq_b))));

	salami->frac_dme = frac_dme;
	salami->z_scr = z_scr;
	salami->rmsd = rmsd;
	salami->tmscore = tmscore;
	salami->andrew_scr = andrew_scr;

	free(scores);
	FREEMEMORY(space, altscores);
	pair_set_destroy(set_sw);
	seq_destroy(seq_a);
	seq_destroy(seq_b);
	coord_destroy(coord_a);
	coord_destroy(coord_b);
	coord_destroy(rmsd_a);
	coord_destroy(rmsd_b);
	score_mat_destroy(matrix);

	return salami;
}

struct salami_info* alignment_aacid(void *config, void *space, const Matchtype *match, IntSequence **s, int len, void *info) {

	double zero_shift = 0.79;
	double gap_open = 8.94;
	double gap_widen = 0.88;

	imbissinfo *imbiss = (imbissinfo*) config;
	massert((imbiss != NULL), "Imbiss info object can not be empty");

	stringset_t* imbiss_info = (stringset_t*) info;
	struct salami_info *salami = malloc(sizeof(*salami));
	massert((salami != NULL), "Can not allocate memory for salami object");

	char *binary = merge(merge(merge(imbiss->binarypath, "/"), sequence_code(s[match->id]->url)), ".bin");

	struct coord *coord_a = coord_read(binary);
	struct coord *coord_b = coord_read(imbiss_info->strings[0].str);

	massert((coord_a != NULL), "Coordinates for sequence A can not be empty");
	massert((coord_b != NULL), "Coordinates for sequence B can not be empty");

	struct seq *seq_a = coord_get_seq(coord_a);
	struct seq *seq_b = coord_get_seq(coord_b);

	massert((seq_a != NULL), "Sequence A can not be empty");
	massert((seq_b != NULL), "Sequence B can not be empty");

	struct score_mat *matrix_score = score_mat_new(seq_size(seq_a), seq_size(seq_b));
	massert((matrix_score != NULL), "Substitution matrix can not be empty");
	struct score_mat *matrix_unknown = score_mat_shift(matrix_score, zero_shift);
	score_mat_destroy(matrix_unknown);

	score_smat(matrix_score, seq_a, seq_b, imbiss->matrix_substitition);

	struct score_mat *crap = NULL;
	struct pair_set *pair_set_nw = score_mat_sum_full(&crap, matrix_score, gap_open, gap_widen, gap_open, gap_widen,
	NULL, NULL, S_AND_W, NULL);

	unsigned id = get_seq_id_simple(pair_set_nw, seq_a, seq_b);
	struct score_struct *scores = get_scores(space, pair_set_nw, coord_a, coord_b, NULL);

	salami->id = (float)id / (float)pair_set_nw->n;
	salami->nw_score = pair_set_nw->score;
	salami->nw_length = (Uint)pair_set_nw->n;
	salami->nw_smpl_score = pair_set_nw->smpl_score;
	salami->nw_score_tot = scores->scr_tot;
	salami->nw_cvr = scores->cvr;
	salami->nw_raw = scores->cvr * seq_size(seq_a);

	score_mat_destroy(crap);
	pair_set_destroy(pair_set_nw);

	FREEMEMORY(space, binary);
	FREEMEMORY(space, scores);

	struct pair_set *pair_set_sw = score_mat_sum_full(&crap, matrix_score, gap_open, gap_widen, gap_open, gap_widen,
	NULL, NULL, S_AND_W, NULL);
	score_mat_destroy(crap);
	scores = get_scores(space, pair_set_sw, coord_a, coord_b, NULL);

	salami->sw_score = pair_set_sw->score;
	salami->sw_length = (Uint)pair_set_sw->n;
	salami->sw_smpl_score = pair_set_sw->smpl_score;
	salami->sw_score_tot = scores->scr_tot;
	salami->sw_cvr = scores->cvr;
	salami->sw_raw = scores->cvr * seq_size(seq_a);

	float frac_dme = get_dme_thresh(pair_set_sw, coord_a, coord_b);
	float *altscores = get_alt_scores(space, 1000, matrix_score, pair_set_sw);

	float dev = 0.0;
	float mean = 0.0;

	normalize_alt_scores(altscores, 1000, &mean, &dev);

	float z_scr = (dev != 0) ? ((scores->scr_tot - mean) / dev) : 0.0;

	float rmsd = 0.0;
	struct coord *rmsd_a = NULL;
	struct coord *rmsd_b = NULL;

	coord_rmsd(pair_set_sw, coord_a, coord_b, 0, &rmsd, &rmsd_a, &rmsd_b);
	double tmscore = tm_score(pair_set_sw, rmsd_a, rmsd_b);

	float andrew_scr = frac_dme * ((scores->cvr * seq_size(seq_a)) / (MIN(seq_size(seq_a), seq_size(seq_b))));

	salami->frac_dme = frac_dme;
	salami->z_scr = z_scr;
	salami->rmsd = rmsd;
	salami->tmscore = tmscore;
	salami->andrew_scr = andrew_scr;

	free(scores);
	FREEMEMORY(space, altscores);
	pair_set_destroy(pair_set_sw);
	seq_destroy(seq_a);
	seq_destroy(seq_b);
	coord_destroy(coord_a);
	coord_destroy(coord_b);
	coord_destroy(rmsd_a);
	coord_destroy(rmsd_b);
	score_mat_destroy(matrix_score);

	return salami;
}
