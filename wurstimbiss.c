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
#include "inc/wurstimbiss.h"

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

/**
 * Build a output string with all parameters
 * this string should be a result of  work of this program
 */
const char * allscores_string(char * picture, IntSequence *sequence_a, IntSequence *sequence_b, Matchtype *matchtype,
		struct salami_info *salami) {
	char * response;

	if (picture == NULL && sequence_a == NULL && sequence_b == NULL && matchtype == NULL && salami == NULL) {
		return "picture;sequence a;sequence b;tm_score;andrew_score;matches;matches_score;matches_swscore;"
				"matches_blast;salami_id;salami_sw_score;salami_sw_smpl_score;salami_sw_score_tot;salami_sw_cvr;"
				"salami_sw_raw;salami_frac_dme;salami_z_scr;salami_rmsd;matchtype_id;description;";
	}

	massert((picture!=NULL), "Match picture can not be empty");
	massert((sequence_a!=NULL), "Sequence a can not be empty");
	massert((sequence_b!=NULL), "Sequence b can not be empty");
	massert((matchtype!=NULL), "Match type can not be empty");
	massert((salami!=NULL), "Salami object can not be empty");

	const char * sequence_a_code = sequence_code(sequence_a->url);
	const char * sequence_b_code = sequence_code(sequence_b->url);

	int strlen = snprintf(NULL, 0, "[%s];%s;%s;%f;%f;%u;%f;%f;%f;%f;%f;%f;%f;%f;%u;%f;%f;%f;%u;%s", picture,
			sequence_a_code, sequence_b_code, salami->tmscore, salami->andrew_scr, matchtype->count, matchtype->score,
			matchtype->swscore, matchtype->blast, salami->id, salami->sw_score, salami->sw_smpl_score,
			salami->sw_score_tot, salami->sw_cvr, salami->sw_raw, salami->frac_dme, salami->z_scr, salami->rmsd,
			matchtype->id, sequence_b->description);

	response = malloc((strlen + 1) * sizeof(*response));
	massert((response!=NULL), "Can not allocate memory for out string");

	snprintf(response, (strlen + 1), "[%s];%s;%s;%f;%f;%u;%f;%f;%f;%f;%f;%f;%f;%f;%u;%f;%f;%f;%u;%s", picture,
			sequence_a_code, sequence_b_code, salami->tmscore, salami->andrew_scr, matchtype->count, matchtype->score,
			matchtype->swscore, matchtype->blast, salami->id, salami->sw_score, salami->sw_smpl_score,
			salami->sw_score_tot, salami->sw_cvr, salami->sw_raw, salami->frac_dme, salami->z_scr, salami->rmsd,
			matchtype->id, sequence_b->description);

	return (const char *) response;
}

/*
 *
 * a handler function for ranked suffix matches
 * accepts pointer to the set of sequences, pointer to to ranked matches
 * the length of the match sequence and its index position
 *
 * returns an integer indicating if the calling function should
 * stop (-1) or proceed (0) reporting matches.
 *
 */
int allscores(void *space, IntSequence *sequence_a, Matchtype *matchtype, IntSequence **sequences, Uint len, Uint match,
		void *info) {
	imbissinfo *imbiss = (imbissinfo*) info;

	if (matchtype->count <= imbiss->minimal_seed) {
		return 0;
	}
	if (match > imbiss->maximal_hit) {
		return -1;
	}

	IntSequence *sequence_b = sequences[matchtype->id];

	char *picture = depictSequence(space, len, 20, matchtype->pos, matchtype->count, '*');
	massert((picture != NULL), "Picture object can not be null");

	struct salami_info *salami = alignment_wurst(info, space, matchtype, sequences, len, imbiss->query);
	massert((salami != NULL), "Salami alignment object can not be null");

	printf("CSV;%s\n", allscores_string(picture, sequence_a, sequence_b, matchtype, salami));

	FREEMEMORY(space, picture);
	FREEMEMORY(space, salami);

	return 1;
}

int main(int argc, char** argv) {
	unsigned char depictsw = 0;

	Uint i, noofqueries = 0;
	void *space = NULL;

	Config *cfg = NULL;
	double *scores = NULL;
	int swscores[2] = { 3, -2 };
	time_t time_start, time_end;

	imbissinfo *imbiss = ALLOCMEMORY(space, NULL, (*imbiss), 1);
	massert((imbiss != NULL), "Imbissinfo object can not be null");

	imbiss->score = NULL;
	imbiss->matrix_substitition = NULL;
	imbiss->consensus = NULL;
	imbiss->swscores = swscores;
	imbiss->handler = (imbissinfo_handler *) allscores;
	imbiss->filter = (imbissinfo_filter *) swconstfilter;
	imbiss->select = (imbissinfo_select *) selectSW;
	imbiss->file_configuration = getopt_configfile(argc, argv);

	assert(ConfigReadFile(imbiss->file_configuration, &cfg) == CONFIG_OK);
	ConfigReadString(cfg, "sources", "file_batch", imbiss->file_batch, sizeof(imbiss->file_batch), 0);
	ConfigReadString(cfg, "sources", "file_substitution", imbiss->file_substitution, sizeof(imbiss->file_substitution), 0);
	ConfigReadString(cfg, "sources", "file_alphabet", imbiss->file_alphabet, sizeof(imbiss->file_alphabet), 0);
	ConfigReadString(cfg, "sources", "file_sequences", imbiss->file_sequences, sizeof(imbiss->file_sequences), 0);
	ConfigReadString(cfg, "sources", "file_logconfig", imbiss->file_logconfig, sizeof(imbiss->file_logconfig), 0);
	ConfigReadString(cfg, "sources", "path_binary", imbiss->path_binary, sizeof(imbiss->path_binary), 0);
	ConfigReadString(cfg, "sources", "path_vector", imbiss->path_vector, sizeof(imbiss->path_vector), 0);


	ConfigReadUnsignedInt(cfg, "limits", "maximal_hit", &imbiss->maximal_hit, 10000);
	ConfigReadUnsignedInt(cfg, "limits", "maximal_match", &imbiss->maximal_match, 50);
	ConfigReadUnsignedInt(cfg, "limits", "minimal_seed", &imbiss->minimal_seed, 4);
	ConfigReadUnsignedInt(cfg, "limits", "minimal_length", &imbiss->minimal_length, 10);
	ConfigFree(cfg);

	assert(zlog_init((const char *)imbiss->file_logconfig) == CONFIG_OK);
	zlog_category_t *logger = zlog_get_category("wurstimbiss");
	zlog_info(logger, "File:\t configuration %s", imbiss->file_configuration);
	zlog_info(logger, "File:\t configuration for logger%s", imbiss->file_logconfig);
	zlog_info(logger, "File:\t alphabet %s", imbiss->file_alphabet);
	zlog_info(logger, "File:\t sequences %s", imbiss->file_sequences);
	zlog_info(logger, "File:\t structures %s", imbiss->file_batch);
	zlog_info(logger, "File:\t substitution matrix %s", imbiss->file_substitution);
	zlog_info(logger, "Max:\t%d matches from suffix array", imbiss->maximal_match);
	zlog_info(logger, "Min:\t%d characters", imbiss->minimal_length);
	zlog_info(logger, "Min:\t%d seeds", imbiss->minimal_seed);

	zlog_debug(logger, "Load:\t%s", imbiss->file_substitution);
	imbiss->matrix_substitition = sub_mat_read(imbiss->file_substitution);
	massert((imbiss->matrix_substitition != NULL), "Substitution matrix can not be empty");

	zlog_debug(logger, "Load:\t%s", imbiss->file_alphabet);
	imbiss->alphabet = alphabet_load_csv(space, imbiss->file_alphabet);
	massert((imbiss->alphabet != NULL), "Alphabet object can not be null");

	Uint sequence_count = 0;
	zlog_debug(logger, "Load:\t%s", imbiss->file_sequences);

	time(&time_start);
	IntSequence **sequences = sequence_load_csv(imbiss, space, imbiss->file_sequences, "", &sequence_count, sequence_salami_load);
	massert((sequences != NULL), "Sequence collection can not be empty");
	time(&time_end);

	zlog_debug(logger, "Time:\t pdb sequences loaded in %f sec", difftime(time_end, time_start));

	time(&time_start);
	zlog_debug(logger, "Build:\t suffix array");
	Suffixarray *suffix_array = suffix_array_init(space, sequences, sequence_count, NULL);
	massert((suffix_array != NULL), "Suffix array can not be empty");
	time(&time_end);

	zlog_debug(logger, "Time:\t suffix array in %f sec", difftime(time_end, time_start));

	/*do search*/
	stringset_t ** queries = readcsv(space, imbiss->file_batch, "", &noofqueries);
	massert((queries != NULL), "Queries collection can not be empty");

	printf("CSV;%s\n", allscores_string(NULL, NULL, NULL, NULL, NULL));
	for (i = 0; i < noofqueries; i++) {

		/*get query form batchfile*/
		char *inputfile = SETSTR(queries[i], 0);

		IntSequence *sequence = sequence_salami_load(imbiss, space, inputfile);
		massert((sequence != NULL), "Sequence object can not be null");

		time(&time_start);
		PairSint *matches = sufSubstring(space, suffix_array, sequence->sequence, sequence->length, imbiss->minimal_length);
		massert((matches != NULL), "Match collection can not be null");
		time(&time_end);

		zlog_debug(logger, "Time:\t suffix array match in %f sec", difftime(time_end, time_start));

		char *vector = merge(merge(merge(imbiss->path_vector, "/"), sequence_code(sequence->url)), ".vec");
		char *binary = merge(merge(merge(imbiss->path_binary, "/"), sequence_code(sequence->url)), ".bin");

		imbiss->query = initStringset(space);
		addString(space, imbiss->query, binary, strlen(binary));
		addString(space, imbiss->query, vector, strlen(vector));

		Uint matches_count = (sequence->length - imbiss->minimal_length);
		zlog_debug(logger, "Count:\t %u matches found", matches_count);

		const char * blast_statistic = getimbissblast(space, sequence, sequences, sequence_count, imbiss->alphabet, imbiss);
		zlog_debug(logger, "Statistic:\t %s", blast_statistic);

		rankSufmatch(space, suffix_array, matches, matches_count, imbiss->maximal_match, imbiss->minimal_length,
				sequences, sequence_count, imbiss->filter, imbiss->select, imbiss->handler, sequence, imbiss, scores,
				depictsw);

		destructSequence(space, sequence);

		FREEMEMORY(space, matches);
	}

	/*final cleanup*/
	for (i = 0; i < sequence_count; i++) {
		destructSequence(space, sequences[i]);
	}
	FREEMEMORY(space, sequences);
	suffix_array_destruct(space, suffix_array);

	zlog_fini();

	exit(EXIT_SUCCESS);
}

