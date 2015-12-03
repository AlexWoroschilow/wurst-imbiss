
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
/*WURST include*/ 
 #include "read_seq_i.h"
 #include "score_mat_i.h"
 #include "prob_vec_i.h"
 #include "prob_vec.h"
 #include "coord_i.h"
 #include "pair_set.h"
 #include "pair_set_i.h"
 #include "score_probvec.h"
 #include "coord.h"
 #include "align_i.h"
 #include "matrix.h"
 #include "model.h"
 #include "cmp_dmat_i.h"
 #include "altscores.h"
 #include "lsqf.h"


 struct score_struct{
 	float scr_tot;
	float cvr;
   
 };


 const float gapopen = 3.25;
 const float gapwide = 0.8942; 
 const float sw1_pgap_open = 3.25;
 const float sw1_pgap_widen = 0.8942;
 const float nw1_pgap_open = 3.25;
 const float nw1_pgap_widen =0.8942 ;
 const int N_AND_W = 0;
 const int S_AND_W = 1;




/*--------------------------- normalize_alt_scores ---------------------------
 *    
 * a port from salami's normalize_alt_scores function.
 * 
 */
 
void
normalize_alt_scores (float *scrs, int len, float *mean, float *dev)
{
   int i;
   float temp;
   
   *mean = 0.0;
   *dev = 0.0;
   
   for(i=0; i < len; i++) {
       *mean += scrs[i];
   }

   (*mean) /= len;
   
   for(i=0; i < len; i++) {
       temp = scrs[i]-(*mean);
	   *dev += temp * temp;
   }

   *dev /= (len-1);
   *dev = sqrt((*dev));
   
   return ;
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
get_alt_scores (void *space, int num_scrs, struct score_mat  *matrix, 
				struct pair_set *set) 
{
    int i;
	float *scrs;
	
	scrs = ALLOCMEMORY(space, NULL, float, num_scrs);

	for (i=0; i < num_scrs; i++) {
		scrs[i]= find_alt_path_score_simple (matrix, set);
	}
	
	return scrs;
}


/*------------------------------ get_dme_thresh ------------------------------
 *    
 * a port from salami's get_dme_thresh function
 * 
 */
 
float
get_dme_thresh (struct pair_set *set, struct coord *a, struct coord *b)
{
    struct coord *model;
	struct seq *seq_a;
	float frac = 0.0;
	int thresh;
	
	seq_a = coord_get_seq(a);
	model = make_model(set, seq_a, b);
	
	if (model == NULL) {
		return 0.0;
	}
	
	if (coord_size(model) < 10) {
	  return 0.0;
	}
	
	thresh = dme_thresh(&frac, a, model, 3.0);
	
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
get_scores (void *space, struct pair_set *set, struct coord *a, struct coord *b, 
	        void *to_use)
{
   	int cover=0,i ;
	float cover_float;
	float score, smpl_score,
		  scr_tot;
    char *pcover1;
	char *pcover2;
	struct seq* seq_a;
	struct score_struct *scores;
	/*float geo_gap = 0, nseq_gap = 0;*/
	float open_cost = 0, widen_cost = 0;
	

	scores = malloc(sizeof(struct score_struct));	  
	score = set->score;
	smpl_score = set->smpl_score;
	
	pair_set_coverage(set, coord_size(a), coord_size(b), 
								&pcover1, &pcover2);
	for (i=0; i < coord_size(a); i++) {
		if(pcover1[i] == '1') cover++;
	}
	
	free(pcover1);
	free(pcover2);
	
	seq_a = coord_get_seq(a);
	cover_float = (float)cover/(float)seq_size(seq_a);
	
	if (cover_float >= 0.05) {
		pair_set_gap(set, &open_cost, &widen_cost, 1, 1);
	}

	scr_tot = smpl_score + sw1_pgap_open*open_cost + 
	  				sw1_pgap_widen*widen_cost;

	seq_destroy(seq_a);
	scores->scr_tot = scr_tot;
	scores->cvr = cover_float;
	
	return scores;
}


/*----------------------------- doWurstAlignment -----------------------------
 *    
 * invokes the Wurst library to calculate alignments and results
 *
 * 
 */
 
struct salami_info*
doWurstAlignment (void *space, Matchtype *m, IntSequence **s, int len, 
				  void *info)
{
    char *p_a;
    char *bin;
    char *vec;
	float frac_dme, z_scr;
	float *altscores;
	unsigned int id;
	float mean, dev, rmsd, andrew_scr;
    double tmscore;
	stringset_t *in;
    struct coord *coord_a, *coord_b, 
				 *rmsd_a = NULL, *rmsd_b=NULL;
	struct seq *seq_a, *seq_b;
	struct prob_vec *pvec_a, *pvec_b;
	struct score_mat  *matrix;
	struct score_struct *scores;
	struct pair_set *set_sw, *set_nw;
	struct score_mat *crap = NULL;
	struct salami_info *salami =NULL;

	salami = malloc(sizeof(struct salami_info));
	
	in = (stringset_t*) info;
	
	/*read coord and sequences*/
    p_a = malloc(100);	
    sprintf(p_a, "/smallfiles/public/no_backup/bm/pdb_all_bin/%5s.bin\0", s[m->id]->url+56);
    bin = p_a;
    coord_a = coord_read(p_a);
	coord_b = coord_read(in->strings[0].str);	
	/*coord_b = coord_read(bin);*/
	
	seq_a = coord_get_seq(coord_a);
	seq_b = coord_get_seq(coord_b);
	
	//FREEMEMORY(space, p_a);
	
	/* read prob vectors and generate score matrix */
	
    sprintf(p_a, "/smallfiles/public/no_backup/bm/pdb_all_vec_6mer_struct/%5s.vec\0", s[m->id]->url+56);
	pvec_a = prob_vec_read(p_a);
	pvec_b = prob_vec_read(in->strings[1].str);
		
	matrix = score_mat_new(seq_size(seq_a), seq_size(seq_b));
	score_pvec(matrix, pvec_a, pvec_b);
	
	prob_vec_destroy(pvec_a);
	prob_vec_destroy(pvec_b);

	FREEMEMORY(space, p_a);
	
	/*nw align*/
	
	set_nw = score_mat_sum_full(&crap, matrix, 
							 sw1_pgap_open, sw1_pgap_widen, 
							 sw1_pgap_open, sw1_pgap_widen, 
							 NULL, NULL, N_AND_W, NULL);
	
	id = get_seq_id_simple(set_nw, seq_a, seq_b);	
	scores = get_scores(space, set_nw, coord_a, coord_b, NULL);
	
	salami->id = (float) id/set_nw->n;
	salami->nw_score = set_nw->score;
	salami->nw_smpl_score = set_nw->smpl_score; 
	salami->nw_score_tot = scores->scr_tot;
	salami->nw_cvr = scores->cvr;
	salami->nw_raw = scores->cvr * seq_size(seq_a);
		
	score_mat_destroy(crap);
	pair_set_destroy(set_nw);
	free(scores);
		
	/*sw align*/
	
	set_sw = score_mat_sum_full(&crap, matrix, 
							 sw1_pgap_open, sw1_pgap_widen, 
							 sw1_pgap_open, sw1_pgap_widen, 
							 NULL, NULL, S_AND_W, NULL);
	score_mat_destroy(crap);
	scores=get_scores(space, set_sw, coord_a, coord_b, NULL);
	
	salami->sw_score = set_sw->score;
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
	
	andrew_scr = frac_dme * ((scores->cvr*seq_size(seq_a))/
								(MIN(seq_size(seq_a),seq_size(seq_b))));
	
	salami->frac_dme= frac_dme;
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

/*----------------------------- latexWurstAlignment -----------------------------
 *    
 * invokes the Wurst library to calculate alignments and results
 *
 * 
 */
 
void
latexWurstAlignment (void *space, Matchtype *m, IntSequence **s, int len, 
				  void *info)
{
    char *p_a;
    char *bin;
    char *vec;
	float frac_dme, z_scr;
	float *altscores;
	float mean, dev, rmsd, andrew_scr;
	stringset_t *in;
    struct coord *coord_a, *coord_b, 
				 *rmsd_a = NULL, *rmsd_b=NULL;
	struct seq *seq_a, *seq_b;
	struct prob_vec *pvec_a, *pvec_b;
	struct score_mat  *matrix;
	struct score_struct *scores;
	struct pair_set *set;
	struct score_mat *crap = NULL;
	
	in = (stringset_t*) info;
	
	p_a = ALLOCMEMORY(space, NULL, char, s[m->id]->urllen);
	strcpy(p_a, s[m->id]->url);
	
	p_a = concatdelim (space, p_a, "bin", s[m->id]->urllen, 3, '.');		
    bin = attachext(space, "/projects/bm/pdb_all_bin/", 25, in->strings[0].str, in->strings[0].len);
    vec = attachext(space, "/projects/bm/pdb_all_vec_6mer_struct/", 37, in->strings[1].str, in->strings[1].len);
    bin = attachext(space, bin, in->strings[0].len+25, ".bin", 4);       
    vec = attachext(space, vec, in->strings[1].len+37, ".vec", 4);
    //printf("#### bin: %s, vec: %s \n", bin, vec);
	
	coord_a = coord_read(p_a);
	/*coord_b = coord_read(in->strings[0].str);*/	
	coord_b = coord_read(bin);	
	FREEMEMORY(space, p_a);
	
	seq_a = coord_get_seq(coord_a);
	seq_b = coord_get_seq(coord_b);

	p_a = ALLOCMEMORY(space, NULL, char, s[m->id]->urllen);
	strcpy(p_a, s[m->id]->url);	
	p_a = concatdelim (space, p_a, "vec", s[m->id]->urllen, 3, '.');	


	pvec_a = prob_vec_read(p_a);
	/*pvec_b = prob_vec_read(in->strings[1].str);*/
	pvec_b = prob_vec_read(vec);
	
	matrix = score_mat_new(seq_size(seq_a), seq_size(seq_b));
	score_pvec(matrix, pvec_a, pvec_b);
	
	prob_vec_destroy(pvec_a);
	prob_vec_destroy(pvec_b);
	
	set = score_mat_sum_full(&crap, matrix, 
							 sw1_pgap_open, sw1_pgap_widen, 
							 sw1_pgap_open, sw1_pgap_widen, 
							 NULL, NULL, N_AND_W, NULL);
	
	scores=get_scores(space, set, coord_a, coord_b, NULL);
	pair_set_destroy(set);
	free(scores);
	score_mat_destroy(crap);
	
	set = score_mat_sum_full(&crap, matrix, 
							 sw1_pgap_open, sw1_pgap_widen, 
							 sw1_pgap_open, sw1_pgap_widen, 
							 NULL, NULL, S_AND_W, NULL);
	score_mat_destroy(crap);
	scores=get_scores(space, set, coord_a, coord_b, NULL);

	
	frac_dme = get_dme_thresh(set, coord_a, coord_b);
	altscores = get_alt_scores(space, 1000, matrix, set);
	normalize_alt_scores(altscores, 1000, &mean, &dev);
	
	/*@todo: get scores from get_scores and calc z-scr here!*/
	z_scr = 0.0;
	if (dev != 0) {
		z_scr = (scores->scr_tot - mean) / dev;
	}

    /*rmsd*/
	coord_rmsd(set, coord_a, coord_b, 0, &rmsd, &rmsd_a, &rmsd_b);

	/*andrew score*/
	/*frac_dme (9) cvr (3)*/
	/*alignmentsize (6) = sw_cvr (3) * seq_size(coord_get_seq(coord_b))*/
	/*andrew_scr (2)= frac_dme(9) * (alignmentsize(6) / minsize) */
	andrew_scr = frac_dme * ((scores->cvr*seq_size(seq_a))/
								(MIN(seq_size(seq_a),seq_size(seq_b))));
	
	printf("%.2f", rmsd);
	
	free(scores);
	FREEMEMORY(space, altscores);
	pair_set_destroy(set);
	FREEMEMORY(space, p_a);
	seq_destroy(seq_a);
	seq_destroy(seq_b);
	coord_destroy(coord_a);
	coord_destroy(coord_b);
	coord_destroy(rmsd_a);
	coord_destroy(rmsd_b);
	score_mat_destroy(matrix);

	return ;
}





