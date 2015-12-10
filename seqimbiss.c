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
	if (m->count <= imbiss->minseeds) {
		return 0;
	}
	if (match > imbiss->noofhits) {
		return -1;
	}

	/*report score stuff*/
	printf("[%d]: score: %f, count: %d\n", match, m->score, m->count);
	printf("%d\t%s\t%d\t", m->id, s[m->id]->url, m->count);
	pic = depictSequence(space, len, 20, m->pos, m->count, '*');
	printf("[%s]\n", pic);
	printf("%s\n", s[m->id]->description);

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
	unsigned char wurst = 0;

	Uint i, noofqueries = 0;
	Uint minseeds = 5;
	Uint maxmatches = 10000;
	imbissinfo *imbiss;
	void *space = NULL;
	double *scores = NULL;

	int swscores[2] = { 3, -2 };
	char *reportfile = NULL;

	int (*handler)(void *, Matchtype *, IntSequence **, Uint, Uint, void *) = allscores;

	double (*filter)(void *, Matchtype *, IntSequence *, IntSequence *, Uint *, Uint, Uint, void *) = swconstfilter;

	Matchtype* (*select)(void *, Matchtype *, Uint k, IntSequence *, IntSequence **, void *) = selectSW;

	stringset_t *queryurl;
	Suffixarray *suffix_array = NULL;
	FAlphabet *alphabet = NULL;
	PairSint *matches = NULL;

	time_t startsuf, endsuf;
	double difsuf, difmatch, difrank;

	Config *cfg = NULL;
	Uint maximal_match, minimal_length;
	char file_batch[1024], file_sub[1024], file_abc[1024], file_seq[1024];
	assert(ConfigReadFile("wurstimbiss.conf", &cfg) == CONFIG_OK);
	ConfigReadString(cfg, "sources", "file_batch", file_batch, sizeof(file_batch), 0);
	ConfigReadString(cfg, "sources", "file_sub", file_sub, sizeof(file_sub), 0);
	ConfigReadString(cfg, "sources", "file_abc", file_abc, sizeof(file_abc), 0);
	ConfigReadString(cfg, "sources", "file_seq", file_seq, sizeof(file_seq), 0);
	ConfigReadUnsignedInt(cfg, "limits", "maximal_match", &maximal_match, 100);
	ConfigReadUnsignedInt(cfg, "limits", "minimal_length", &minimal_length, 10);
	ConfigFree(cfg);

	zlog_category_t *logger;
	assert(zlog_init("wurstimblog.conf") == CONFIG_OK);
	logger = zlog_get_category("wurstimbiss");
	zlog_info(logger, "File abc:\t%s", file_abc);
	zlog_info(logger, "File seq:\t%s", file_seq);
	zlog_info(logger, "File str:\t%s", file_batch);
	zlog_info(logger, "File sub:\t%s", file_sub);
	zlog_info(logger, "Max:\t%d matches from suffix array", maximal_match);
	zlog_info(logger, "Min:\t%d characters", minimal_length);

	imbiss = ALLOCMEMORY(space, NULL, imbissinfo, 1);
	imbiss->reportfile = reportfile;
	imbiss->swscores = swscores;
	imbiss->noofhits = maximal_match;
	imbiss->minseeds = minseeds;
	imbiss->wurst = wurst;

	time_t time_start, time_end;

	zlog_debug(logger, "Load:\t%s", file_abc);
	alphabet = alphabet_load_csv(space, file_abc);

	Uint sequence_count = 0;
	zlog_debug(logger, "Load:\t%s", file_seq);

	time(&time_start);
	IntSequence **sequences = sequence_load_csv(space, file_seq, "", &sequence_count, sequence_aacid_load);
	time(&time_end);

	zlog_debug(logger, "Time:\t pdb sequences loaded in %f sec", difftime(time_end, time_start));

	time(&time_start);
	suffix_array = suffix_array_init(space, sequences, sequence_count, NULL);
	time(&time_end);

	zlog_debug(logger, "Time:\t suffix array in %f sec", difftime(time_end, time_start));

	/*do search*/
	stringset_t ** queries = readcsv(space, file_batch, "", &noofqueries);
	for (i = 0; i < noofqueries; i++) {

		/*get query form batchfile*/
		char *inputfile = SETSTR(queries[i], 0);

		IntSequence *sequence = sequence_aacid_load(space, inputfile);
		sequence_dump_aacid(sequence);

		time(&time_start);
		matches = sufSubstring(space, suffix_array, sequence->sequence, sequence->length, minimal_length);
		time(&time_end);

		zlog_debug(logger, "Time:\t suffix array match in %f sec", difftime(time_end, time_start));

		char *vector = malloc(sizeof(char) * 66);
		sprintf(vector, "/smallfiles/public/no_backup/bm/pdb_all_vec_6mer_struct/%5s.vec\0", sequence->url + 56);
		//char *vector = scr_printf("/smallfiles/public/no_backup/bm/pdb_all_vec_6mer_struct/%5s.vec\0", sequence->url + 56);

		char *binary = malloc(sizeof(char) * 54);
		sprintf(binary, "/smallfiles/public/no_backup/bm/pdb_all_bin/%5s.bin\0", sequence->url + 56);
		//char *binary = scr_printf("/smallfiles/public/no_backup/bm/pdb_all_bin/%5s.bin\0", sequence->url + 56);

		queryurl = initStringset(space);
		addString(space, queryurl, binary, strlen(binary));
		addString(space, queryurl, vector, strlen(vector));

		imbiss->query = queryurl;
		imbiss->substrlen = minimal_length;
		imbiss->alphabet = alphabet;

		rankSufmatch(space, suffix_array, matches, (sequence->length - minimal_length), maxmatches, minimal_length,
				sequences, sequence_count, filter, select, handler, sequence, imbiss, scores, depictsw);

		destructSequence(space, sequence);

		FREEMEMORY(space, imbiss->score);
		FREEMEMORY(space, matches);
	}

	/*final cleanup*/
	for (i = 0; i < sequence_count; i++) {
		destructSequence(space, sequences[i]);
	}
	FREEMEMORY(space, sequences);
	suffix_array_destruct(space, suffix_array);

	zlog_fini();

	printf("Goodbye.\n");
	return EXIT_SUCCESS;
}

