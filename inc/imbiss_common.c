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
#include <massert.h>

#include "salami.h"
#include "sufmatch.h"
#include "intsequence.h"
#include "wurstimbiss.h"
#include "depictseqs.h"
#include "logger.h"

/**
 * Build a output string with all parameters
 * this string should be a result of  work of this program
 */
const char * allscores_string(char * picture, IntSequence *sequence_a, IntSequence *sequence_b, const Matchtype *matchtype,
		struct salami_info *salami) {
	char * response;

	if (picture == NULL && sequence_a == NULL && sequence_b == NULL && matchtype == NULL && salami == NULL) {
		return "picture;sequence a;sequence b;tm_score;andrew_score;matches;matches_score;matches_swscore;"
				"matches_blast;salami_id;salami_sw_score;salami_sw_smpl_score;salami_sw_score_tot;salami_sw_cvr;"
				"salami_sw_raw;salami_frac_dme;salami_z_scr;salami_rmsd;matchtype_id;description;";
	}

	massert((picture != NULL), "Match picture can not be empty");
	massert((sequence_a != NULL), "Sequence a can not be empty");
	massert((sequence_b != NULL), "Sequence b can not be empty");
	massert((matchtype != NULL), "Match type can not be empty");
	massert((salami != NULL), "Salami object can not be empty");

	const char * sequence_a_code = sequence_code(sequence_a->url);
	const char * sequence_b_code = sequence_code(sequence_b->url);

	int strlen = snprintf(NULL, 0, "[%s];%s;%s;%f;%f;%u;%f;%f;%f;%f;%f;%f;%f;%f;%u;%f;%f;%f;%u;%s", picture,
			sequence_a_code, sequence_b_code, salami->tmscore, salami->andrew_scr, matchtype->count, matchtype->score,
			matchtype->swscore, matchtype->blast, salami->id, salami->sw_score, salami->sw_smpl_score,
			salami->sw_score_tot, salami->sw_cvr, salami->sw_raw, salami->frac_dme, salami->z_scr, salami->rmsd,
			matchtype->id, sequence_b->description);

	response = malloc((strlen + 1) * sizeof(*response));
	massert((response != NULL), "Can not allocate memory for out string");

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
int allscores_wurst(void *space, IntSequence *sequence_a, const Matchtype *matchtype, IntSequence **sequences, Uint len,
		Uint match, void *info) {
	time_t time_start, time_end;
	imbissinfo *imbiss = (imbissinfo*) info;
	massert((imbiss != NULL), "Imbiss info object can not be null");

	if (matchtype->count <= imbiss->minimal_seed) {
		return 0;
	}
	if (match > imbiss->maximal_hit) {
		return -1;
	}

	IntSequence *sequence_b = sequences[matchtype->id];

	char *picture = depictSequence(space, len, 20, matchtype->pos, matchtype->count, '*');
	massert((picture != NULL), "Picture object can not be null");

	time(&time_start);
	struct salami_info *salami = alignment_wurst(info, space, matchtype, sequences, len, imbiss->query);
	massert((salami != NULL), "Salami alignment object can not be null");
	time(&time_end);

	logger_info("Time:\t alignment_wurst in %f sec", difftime(time_end, time_start));

	printf("CSV;%s\n", allscores_string(picture, sequence_a, sequence_b, matchtype, salami));

	FREEMEMORY(space, picture);
	FREEMEMORY(space, salami);

	return 1;
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
int allscores_aacid(void *space, IntSequence *sequence_a, Matchtype *matchtype, IntSequence **sequences, Uint len,
		Uint match, void *info) {

	time_t time_start, time_end;
	imbissinfo *imbiss = (imbissinfo*) info;
	massert((imbiss != NULL), "Imbiss info object can not be null");

	if (matchtype->count <= imbiss->minimal_seed) {
		return 0;
	}
	if (match > imbiss->maximal_hit) {
		return -1;
	}

	IntSequence *sequence_b = sequences[matchtype->id];
	massert((sequence_b != NULL), "Sequence object can not be null");

	char *picture = depictSequence(space, len, 20, matchtype->pos, matchtype->count, '*');
	massert((picture != NULL), "Picture object can not be null");

	time(&time_start);
	struct salami_info *salami = alignment_aacid(info, space, matchtype, sequences, len, imbiss->query);
	massert((salami != NULL), "Salami alignment object can not be null");
	time(&time_end);

	logger_info("Time:\t allscores_aacid in %f sec", difftime(time_end, time_start));

	printf("CSV;%s\n", allscores_string(picture, sequence_a, sequence_b, matchtype, salami));

	FREEMEMORY(space, picture);
	return 1;
}

