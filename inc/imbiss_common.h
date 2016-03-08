#ifndef IMBISS_COMMON_H
#define IMBISS_COMMON_H

#include "salami.h"
#include "sufmatch.h"
#include "intsequence.h"
#include "wurstimbiss.h"


const char * allscores_string(char * picture, IntSequence *sequence_a, IntSequence *sequence_b, Matchtype *matchtype,
		struct salami_info *salami);

int allscores_wurst(void *space, IntSequence *sequence_a, Matchtype *matchtype, IntSequence **sequences, Uint len,
		Uint match, void *info);

int allscores_aacid(void *space, IntSequence *sequence_a, Matchtype *matchtype, IntSequence **sequences, Uint len,
		Uint match, void *info);

#endif
