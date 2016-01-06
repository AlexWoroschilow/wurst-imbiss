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
#include "inc/imbiss_common.h"
#include "inc/imbiss_getopt.h"
#include "inc/logger.h"

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
	imbiss->handler = (imbissinfo_handler) allscores_wurst;
	imbiss->loader = (imbissinfo_loader) sequence_salami_load;
	imbiss->filter = (imbissinfo_filter) swconstfilter;
	imbiss->select = (imbissinfo_select) selectSW;
	imbiss->file_configuration = getopt_configfile(argc, argv, "wurstimbiss.cnf");

	assert(ConfigReadFile(imbiss->file_configuration, &cfg) == CONFIG_OK);
	ConfigReadString(cfg, "common", "type", imbiss->type, sizeof(imbiss->type), "salami");
	ConfigReadString(cfg, "sources", "logger", imbiss->logger, sizeof(imbiss->logger), "logger.cnf");
	ConfigReadString(cfg, "sources", "structlist", imbiss->structlist, sizeof(imbiss->structlist), "/pdb_struct.lib");
	ConfigReadString(cfg, "sources", "substitution", imbiss->substitution, sizeof(imbiss->substitution), "blosum62.mat");
	ConfigReadString(cfg, "sources", "alphabet", imbiss->alphabetfile, sizeof(imbiss->alphabetfile), "alphabet_t500");
	ConfigReadString(cfg, "sources", "sequencelist", imbiss->sequencelist, sizeof(imbiss->sequencelist), "pdb_seq.list");
	ConfigReadString(cfg, "sources", "binarypath", imbiss->binarypath, sizeof(imbiss->binarypath), "./bin");
	ConfigReadString(cfg, "sources", "vectorpath", imbiss->vectorpath, sizeof(imbiss->vectorpath), "./vec");
	ConfigReadUnsignedInt(cfg, "limits", "maximal_hits", &imbiss->maximal_hit, 10000);
	ConfigReadUnsignedInt(cfg, "limits", "maximal_match", &imbiss->maximal_match, 50);
	ConfigReadUnsignedInt(cfg, "limits", "minimal_seed", &imbiss->minimal_seed, 4);
	ConfigReadUnsignedInt(cfg, "limits", "minimal_length", &imbiss->minimal_length, 10);
	ConfigFree(cfg);

	if (strcmp(imbiss->type, "aacid") == 0) {
		imbiss->handler = (imbissinfo_handler) allscores_aacid;
		imbiss->loader = (imbissinfo_loader) sequence_aacid_load;
		imbiss->filter = (imbissinfo_filter) swconstfilter;
		imbiss->select = (imbissinfo_select) selectSW;
	}

	logger_init((const char *) imbiss->logger, "wurstimbiss");
	logger_info("File:\t type %s", imbiss->type);
	logger_info("File:\t configuration %s", imbiss->file_configuration);
	logger_info("File:\t configuration for logger%s", imbiss->logger);
	logger_info("File:\t alphabet %s", imbiss->alphabetfile);
	logger_info("File:\t sequences %s", imbiss->sequencelist);
	logger_info("File:\t structures %s", imbiss->structlist);
	logger_info("File:\t substitution matrix %s", imbiss->substitution);
	logger_info("Max:\t %d matches from suffix array", imbiss->maximal_match);
	logger_info("Min:\t %d characters", imbiss->minimal_length);
	logger_info("Min:\t %d seeds", imbiss->minimal_seed);

	logger_info("Load:\t%s", imbiss->substitution);
	imbiss->matrix_substitition = sub_mat_read(imbiss->substitution);
	massert((imbiss->matrix_substitition != NULL), "Substitution matrix can not be empty");

	logger_info("Load:\t%s", imbiss->alphabetfile);
	imbiss->alphabet = alphabet_load_csv(space, imbiss->alphabetfile);
	massert((imbiss->alphabet != NULL), "Alphabet object can not be null");

	Uint sequence_count = 0;
	logger_info("Load:\t%s", imbiss->sequencelist);

	time(&time_start);
	IntSequence **sequences = sequence_load_csv(imbiss, space, imbiss->sequencelist, "", &sequence_count, imbiss->loader);
	massert((sequences != NULL), "Sequence collection can not be empty");
	time(&time_end);

	logger_info("Time:\t pdb sequences loaded in %f sec", difftime(time_end, time_start));

	time(&time_start);
	logger_info("Build:\t suffix array");
	Suffixarray *suffix_array = suffix_array_init(space, sequences, sequence_count, NULL);
	massert((suffix_array != NULL), "Suffix array can not be empty");
	time(&time_end);

	logger_info("Time:\t suffix array in %f sec", difftime(time_end, time_start));

	/*do search*/
	stringset_t ** queries = readcsv(space, imbiss->structlist, "", &noofqueries);
	massert((queries != NULL), "Queries collection can not be empty");

	printf("CSV;%s\n", allscores_string(NULL, NULL, NULL, NULL, NULL));
	for (i = 0; i < noofqueries; i++) {

		/*get query form batchfile*/
		char *inputfile = SETSTR(queries[i], 0);

		IntSequence *sequence = imbiss->loader(imbiss, space, inputfile);
		massert((sequence != NULL), "Sequence object can not be null");

		time(&time_start);
		PairSint *matches = sufSubstring(space, suffix_array, sequence->sequence, sequence->length, imbiss->minimal_length);
		massert((matches != NULL), "Match collection can not be null");
		time(&time_end);

		logger_info("Time:\t suffix array match in %f sec", difftime(time_end, time_start));

		char *vector = merge(merge(merge(imbiss->vectorpath, "/"), sequence_code(sequence->url)), ".vec");
		char *binary = merge(merge(merge(imbiss->binarypath, "/"), sequence_code(sequence->url)), ".bin");

		imbiss->query = initStringset(space);
		addString(space, imbiss->query, binary, strlen(binary));
		addString(space, imbiss->query, vector, strlen(vector));

		Uint matches_count = (sequence->length - imbiss->minimal_length);
		logger_info("Count:\t %u matches found", matches_count);

		if (strcmp(imbiss->type, "salami") == 0) {

			logger_info("Statistic:\t %s", getimbissblast( //
					space, //
					sequence, //
					sequences, //
					sequence_count, //
					imbiss->alphabet, //
					imbiss //
					));
		}

		time(&time_start);
		rankSufmatch( //
				space, //
				suffix_array, //
				matches, //
				matches_count, //
				imbiss->maximal_match, //
				imbiss->minimal_length, //
				sequences, //
				sequence_count, //
				imbiss->filter, //
				imbiss->select, //
				imbiss->handler, sequence, //
				imbiss, //
				scores, //
				depictsw //
				);
		time(&time_end);
		logger_info("Time:\t suffix match rank in %f sec", difftime(time_end, time_start));

		destructSequence(space, sequence);

		FREEMEMORY(space, matches);

		logger_info("Done:\t %f %", ((float )i / (float )noofqueries * 100.0));
	}

	/*final cleanup*/
	for (i = 0; i < sequence_count; i++) {
		destructSequence(space, sequences[i]);
	}
	FREEMEMORY(space, sequences);
	suffix_array_destruct(space, suffix_array);

	logger_destroy();

	exit(EXIT_SUCCESS);
}

