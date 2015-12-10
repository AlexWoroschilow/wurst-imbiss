/*
 *  wurstimbiss.c
 *  perform a boyer-moore exact string
 *  matching on probability sequences and score
 *  the matches w/ Wurst additionally
 *
 *  @author Steve Hoffmann
 *  @email shoffmann@zbh.uni-hamburg.de
 *  @date 12/07/06 22:55:11 CET
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
#include <ncurses.h>
#include <memory.h>
#include <stddef.h>
#include <assert.h>
#include <zlog.h>
#include <configini.h>

#include "inc/massert.h"
#include "lib/memory.h"
#include "inc/wurstimbiss.h"
#include "inc/blaststat.h"
#include "inc/imbissblast.h"
#include "inc/prob_vec.h"
#include "inc/prob_vec_i.h"
#include "inc/falphabet.h"
#include "inc/intsequence.h"

#include "lib/mathematics.h"
#include "lib/stringutils.h"
#include "lib/fileio.h"
#include "inc/createalphabet.h"
#include "inc/sufarray.h"
#include "inc/sufmatch.h"
#include "inc/mm.h"
#include "inc/encodeprobvec.h"
#include "inc/depictseqs.h"
#include "inc/salami.h"
#include "inc/probseqscore.h"
#include "inc/usage.h"
#include "lib/dpalign.h"

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
#include "lib/sort.h"
#include "gnuplot_i.h"

static int verbose_flag;

static struct option long_options[] = {
/* These options set a flag. */
{ "verbose", no_argument, &verbose_flag, 1 }, { "brief", no_argument, &verbose_flag, 0 }, { "sequences",
required_argument, 0, 's' }, { "query", required_argument, 0, 'q' }, { "matchlen", required_argument, 0, 'l' }, {
		"minmatches", required_argument, 0, 'c' }, { "depictsw", no_argument, 0, 'd' },
		{ "veggie", no_argument, 0, 'w' }, { "batchfile", required_argument, 0, 'b' }, { "percent", required_argument,
				0, 'p' }, { "subfile", required_argument, 0, 'x' }, { "reportfile", required_argument, 0, 'r' }, {
				"allscores",
				no_argument, 0, 'S' }, { "segmentscore", no_argument, 0, 'B' }, { "allalign", no_argument, 0, 'A' }, {
				"scorefilter", no_argument, 0, 'F' }, { "segfilter", no_argument, 0, 'G' }, { "match",
		required_argument, 0, 'M' }, { "mismatch", required_argument, 0, 'D' }, { "latex", no_argument, 0, 'L' }, {
				"maxmatches", required_argument, 0, 'm' }, { "help", no_argument, 0, 'h' }, { "gnuplot", no_argument, 0,
				'g' }, { 0, 0, 0, 0 } };

/*------------------------------- getProbChar --------------------------------
 *    
 * returns the position of character from a probability sequence over the 
 * given alphabet of size asize
 * 
 */

Uint getProbChar(void *alphabet, Uint aszize, void *p, Uint pos) {
	/*at this point of the development alphabet is a stub*/
	/*FAlphabet *a = (FAlphabet *)alphabet;*/
	IntSequence *seq = (IntSequence *) p;
	Uint ch = seq->sequence[pos];

	return ch;
}

int latexscores(void *space, Matchtype *m, IntSequence **s, Uint len, Uint match, void *info) {
	int sw;
	double explambda, E;

	/*FILE* fp;*/

	struct salami_info* salami;
	imbissinfo *imbiss;
	stringset_t *query;

	imbiss = (imbissinfo*) info;

	if (m->count <= imbiss->minseeds)
		return 0;
	if (match > imbiss->noofhits)
		return -1;

	query = tokensToStringset(space, "/.", s[m->id]->url, strlen(s[m->id]->url));

	explambda = exp(-imbiss->lambda * m->blast);
	E = imbiss->K * 2500000 * imbiss->substrlen * explambda;
	sw = floor(m->swscore);

	salami = doWurstAlignment(space, m, s, len, imbiss->query);

	printf("%s & %s & %.2f & %d & %.2f & %.2f & %.2f & %.2f\\\\ \n", query->strings[query->noofstrings - 1].str,
			s[m->id]->description, m->score, sw, salami->sw_score_tot, salami->id, log10(E), salami->rmsd);

	/* latexWurstAlignment(space, m, s, len, ((imbissinfo*)info)->query);
	 */
	FREEMEMORY(space, salami);
	destructStringset(space, query);

	return 1;
}

void consensus(void *space, Uint *consensus, Uint len, IntSequence *query, Uint seedlen, void* info) {

	Uint i, j;
	char *constr;
	double sum, norm;
	double *consensus_double;
	double *seedprobability;

	imbissinfo *imbiss;
	IntSequence *consensusSequence;
	gnuplot_ctrl *h;

	imbiss = (imbissinfo*) info;

	seedprobability = ALLOCMEMORY(space, NULL, double, len);
	for (i = 0; i < len; i++) {
		for (sum = 0.0, j = 0; j < seedlen; j++) {
			sum += imbiss->score[(Uint) query->sequence[i + j]];
		}
		seedprobability[i] = log10(imbiss->K * 2500000 * seedlen * exp(-(double) imbiss->lambda * sum));
	}

	consensusSequence = initSequence(space);
	consensusSequence->sequence = consensus;
	consensusSequence->length = len;

	constr = printSequence(space, consensusSequence, 60);
	printf("%s\n", constr);
	FREEMEMORY(space, constr);

	consensus_double = ALLOCMEMORY(space, NULL, double, len);

	h = gnuplot_init();
	gnuplot_setstyle(h, "lines");
	sum = 0;
	for (i = 0; i < len; i++) {
		sum += consensus[i];
	}

	for (i = 0; i < len; i++) {
		norm = consensus[i];
		norm = (norm > 0) ? norm : 1;
		consensus_double[i] = log10(norm);

		/*	log10(imbiss->K*3000000*len*exp(-(double)consensus[i]));
		 if (consensus_double[i] < -400) consensus_double[i]=-400;*/
	}

	gnuplot_cmd(h, "set title 'IMBISS - seed statistics' -28,0 font'Helvetica,15'");
	gnuplot_cmd(h, "set label '%s' at screen 0.12,0.92 font 'Helvetica,12'", imbiss->query->strings[0].str);
	gnuplot_cmd(h, "set label 'seed length: %d' at graph 0.05,0.95 font 'Helvetica, 12'", seedlen);
	gnuplot_set_xlabel(h, "query sequence residue");
	gnuplot_set_ylabel(h, "log");
	gnuplot_plot_x(h, consensus_double, len, "log(number of matches)");
	gnuplot_plot_x(h, seedprobability, len, "log(Kmn*e^{lambda*score(seed)})");

	FREEMEMORY(space, seedprobability);
	FREEMEMORY(space, consensus_double);
	FREEMEMORY(space, consensusSequence);
}

/*-------------------------------- allscores ---------------------------------
 *    
 * a handler function for ranked suffix matches
 * accepts pointer to the set of sequences, pointer to to ranked matches
 * the length of the match sequence and its index position
 *
 * returns an integer indicating if the calling function should 
 * stop (-1) or proceed (0) reporting matches.
 * 
 */

int allscores(void *space, Matchtype *m, IntSequence **s, Uint len, Uint match, void *info) {
	char *pic;
	float rmsd = -1;
	double explambda, E;
	FILE* fp;
	imbissinfo *imbiss;
	struct salami_info *salami;
	stringset_t *query;

	imbiss = (imbissinfo*) info;
	if (m->count <= imbiss->minseeds)
		return 0;
	if (match > imbiss->noofhits)
		return -1;

	/*report score stuff*/
	printf("[%d]: score: %f, count: %d\n", match, m->score, m->count);
	printf("%d\t%s\t%d\t", m->id, s[m->id]->url, m->count);
	pic = depictSequence(space, len, 20, m->pos, m->count, '*');
	printf("[%s]\n", pic);
	printf("%s\n", s[m->id]->description);

	if (imbiss->wurst) {
		salami = doWurstAlignment(space, m, s, len, imbiss->query);
		printf("sequence identity: %f\n", salami->id);
		printf("scr: %f (%f), scr_tot: %f, cvr: %f (raw: %d)\n", salami->sw_score, salami->sw_smpl_score,
				salami->sw_score_tot, salami->sw_cvr, salami->sw_raw);
		printf("frac_dme: %f, z_scr: %f, rmsd: %f, andrew_scr %f\n", salami->frac_dme, salami->z_scr, salami->rmsd,
				salami->andrew_scr);
		printf("tm_scr %f\n", salami->tmscore);
		FREEMEMORY(space, salami);
	}

	printf("gapless sw score: %f\n", m->swscore);

	/*report blast stuff*/
	printf("highest seed score (HSS): %f\n", m->blast);
	/*printf("lambda*S %19.16e\n", m->blast *((imbissinfo*)info)->lambda);*/
	explambda = exp(-imbiss->lambda * m->blast);
	/*printf("exp(-lambda*S): %19.16e\n", explambda);	*/
	E = ((imbissinfo*) info)->K * 2500000 * imbiss->substrlen * explambda;
	/*printf("E=Kmn * exp(-lambda*S): %19.16e\n", E);*/
	printf("log(HSS): %f\n", log10(E));
	printf("1-exp(-HSS): %19.16e\n", 1 - exp(-E));

	FREEMEMORY(space, pic);
	return 1;
}

/*----------------------------------- main -----------------------------------
 *    
 *  main function
 * 
 */

int main(int argc, char** argv) {
	Sint optindex, c;
	unsigned char depictsw = 0;
	unsigned char wurst = 1;

	Uint i, noofqueries = 0;
	Uint noofhits = 100;
	Uint substrlen = 10;
	Uint minseeds = 5;
	Uint maxmatches = 10000;
	char *vec, *bin;
	imbissinfo *imbiss;
	void *space = NULL;
	double *scores = NULL;

	int swscores[2] = { 3, -2 };
	char *inputfile = NULL;
	char *reportfile = NULL;

	int (*handler)(void *, Matchtype *, IntSequence **, Uint, Uint, void *) = allscores;

	double (*filter)(void *, Matchtype *, IntSequence *, IntSequence *, Uint *, Uint, Uint, void *) = swconstfilter;

	Matchtype* (*select)(void *, Matchtype *, Uint k, IntSequence *, IntSequence **, void *) = selectSW;

	stringset_t *queryurl, **queries = NULL;
	Suffixarray *suffix_array = NULL;
	IntSequence *input = NULL;
	FAlphabet *alphabet = NULL;
	PairSint *matches = NULL;

	time_t startsuf, endsuf;
	double difsuf, difmatch, difrank;

	Config *cfg = NULL;
	char file_batch[1024], file_sub[1024], file_abc[1024], file_seq[1024];
	assert(ConfigReadFile("wurstimbiss.conf", &cfg) == CONFIG_OK);
	ConfigReadString(cfg, "sources", "file_batch", file_batch, sizeof(file_batch), 0);
	ConfigReadString(cfg, "sources", "file_sub", file_sub, sizeof(file_sub), 0);
	ConfigReadString(cfg, "sources", "file_abc", file_abc, sizeof(file_abc), 0);
	ConfigReadString(cfg, "sources", "file_seq", file_seq, sizeof(file_seq), 0);

	ConfigFree(cfg);

	zlog_category_t *logger;
	assert(zlog_init("wurstimblog.conf") == CONFIG_OK);
	logger = zlog_get_category("wurstimbiss");
	zlog_info(logger, "File abc:\t%s", file_abc);
	zlog_info(logger, "File seq:\t%s", file_seq);
	zlog_info(logger, "File str:\t%s", file_batch);
	zlog_info(logger, "File sub:\t%s", file_sub);

	imbiss = ALLOCMEMORY(space, NULL, imbissinfo, 1);
	imbiss->reportfile = reportfile;
	imbiss->swscores = swscores;
	imbiss->noofhits = noofhits;
	imbiss->minseeds = minseeds;
	imbiss->wurst = wurst;

	time_t time_start, time_end;

	queries = readcsv(space, file_batch, "", &noofqueries);

	/*read alphabet*/
	zlog_debug(logger, "Load:\t%s", file_abc);
	alphabet = alphabet_load_csv(space, file_abc);

	Uint sequence_count = 0;
	zlog_debug(logger, "Load:\t%s", file_seq);

	/*load sequence database*/
	time(&time_start);
	IntSequence **sequences_wurst = sequence_load_csv(space, file_seq, "", &sequence_count, sequence_aacid_load);
	time(&time_end);

	zlog_debug(logger, "Time:\t pdb sequences loaded in %f sec", difftime(time_end, time_start));

	time(&time_start);
	suffix_array = suffix_array_init(space, sequences_wurst, sequence_count, NULL);
	time(&time_end);

	zlog_debug(logger, "Time:\t suffix array in %f sec", difftime(time_end, time_start));

	/*do search*/
	for (i = 0; i < noofqueries; i++) {

		/*get query form batchfile*/
		inputfile = SETSTR(queries[i], 0);

		input = sequence_aacid_load(space, inputfile);
		const char * sequence_printable = sequence_print(space, input, 60);
		zlog_debug(logger, "Sequence:\n %s", sequence_printable);

		time(&time_start);
		matches = sufSubstring(space, suffix_array, input->sequence, input->length, substrlen);
		time(&time_end);

		zlog_debug(logger, "Time:\t suffix array match in %f sec", difftime(time_end, time_start));

	}

	/*final cleanup*/
	for (i = 0; i < sequence_count; i++) {
		destructSequence(space, sequences_wurst[i]);
	}
	FREEMEMORY(space, sequences_wurst);
	suffix_array_destruct(space, suffix_array);

	zlog_fini();

	printf("Goodbye.\n");
	return EXIT_SUCCESS;
}

