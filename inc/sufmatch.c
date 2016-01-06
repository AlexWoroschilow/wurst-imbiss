/*
 *  sufmatch.c
 *  functions for matching in suffixarrays
 *
 *  @author Steve Hoffmann
 *  @email shoffmann@zbh.uni-hamburg.de
 *  @date 12/19/06 15:13:18 CET
 *  
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "basic-types.h"
#include "memory.h"
#include "mathematics.h"
#include "sufarray.h"
#include "sufmatch.h"
#include "mm.h"
#include "intsequence.h"
#include "list.h"
#include "depictseqs.h"
#include "dpalign.h"
#include "cantor.h"
#include "wurstimbiss.h"
#include "falphabet.h"
#include "massert.h"

/*------------------------------- suffixscore --------------------------------
 *    
 * scoring function for sw alignments of structure sequences 
 * 
 */

int suffixscore(symtype a, symtype b, void* info) {
	return (a == b) ? 3 : -2;
}

/*------------------------------- sufSubstring -------------------------------
 *    
 * retrieve the longest substring matches in a suffixarray
 * 
 */

PairSint*
sufSubstring(void *space, Suffixarray *arr, Uint *pattern, Uint len, Uint sublen) {
	Uint i;
	PairSint *res, d;

	if (len <= sublen) {
		return NULL;
	}

	res = ALLOCMEMORY(space, NULL, PairSint, len - sublen);
	for (i = 0; i < len - sublen; i++) {
		d = mmsearch(arr, &pattern[i], sublen, 0, 0, arr->numofsuffixes - 1);
		res[i].a = d.a;
		res[i].b = d.b;
	}

	return res;
}

/*------------------------------ reportSufmatch ------------------------------
 *    
 * returns a beautified match string
 * 
 */

void reportSufmatch(Suffixarray *a, PairSint *matches, Uint len, Uint threshold, IntSequence **s) {
	Uint i, j, idx;
	char *info;

	for (i = 0; i < len; i++) {
		if (matches[i].b >= ((matches[i].a) + threshold)) {
			/*valid matches*/
			for (j = matches[i].a; j <= matches[i].b; j++) {
				idx = getMultiSeqIndex(a->seq, a->suffixptr[a->suftab[j]]);
				info = s[idx]->description;
				printf("[%d]:\t %s\n", j, info);
			}
		}
	}

	return;
}

/*---------------------------------- cmp_* -----------------------------------
 *    
 * compare functions for clibs's qsort and bsearch
 * 1. cmp_suffixno : used in rankSufmatch to sort sequence numbers
 * 2. cmp_ranks    : used in rankSufmatch to sort sequence ranks
 * 3. cmp_ranks_ptr: used in rankSufmatchList to sort sequence ranks
 * 
 */

int cmp_blast(const void *a, const void *b) {
	Matchtype *first = (Matchtype*) a;
	Matchtype *second = (Matchtype*) b;
	double frac_first, frac_second;

	frac_first = (double) first->blast;
	frac_second = (double) second->blast;

	if (frac_first > frac_second)
		return 1;
	if (frac_first < frac_second)
		return -1;

	return 0;
}

int cmp_swscore(const void *a, const void *b) {
	Matchtype *first = (Matchtype*) a;
	Matchtype *second = (Matchtype*) b;
	double frac_first, frac_second;

	if (first->swscore == 0 && second->swscore == 0) {
		if (first->count > second->count)
			return 1;
		if (first->count < second->count)
			return -1;
	}

	frac_first = (double) first->swscore;
	frac_second = (double) second->swscore;

	if (frac_first > frac_second)
		return 1;
	if (frac_first < frac_second)
		return -1;

	return 0;
}

int cmp_score(const void *a, const void *b) {
	Matchtype *first = (Matchtype*) a;
	Matchtype *second = (Matchtype*) b;
	double frac_first, frac_second;

	if (first->score == 0 && second->score == 0) {
		if (first->count > second->count)
			return 1;
		if (first->count < second->count)
			return -1;
	}

	frac_first = (double) first->score;
	frac_second = (double) second->score;

	if (frac_first > frac_second)
		return 1;
	if (frac_first < frac_second)
		return -1;

	return 0;
}

int cmp_suffixno(const void *a, const void* b) {
	Matchtype *first = (Matchtype*) a;
	Matchtype *second = (Matchtype*) b;

	if (first->id > second->id)
		return 1;
	if (first->id < second->id)
		return -1;

	return 0;
}

int cmp_ranks(const void *a, const void* b) {
	Matchtype *first = (Matchtype*) a;
	Matchtype *second = (Matchtype*) b;

	if (first->count > second->count)
		return 1;
	if (first->count < second->count)
		return -1;

	return 0;
}

int cmp_rank_ptr(const void *a, const void* b) {
	Matchtype **first = (Matchtype**) a;
	Matchtype **second = (Matchtype**) b;

	if (first[0]->count > second[0]->count)
		return 1;
	if (first[0]->count < second[0]->count)
		return -1;

	return 0;
}

/*------------------------------ freeMatchtype -------------------------------
 *    
 * a function to delete Matchtypes in a list
 * 
 */

void freeMatchtype(void *space, void *data) {
	Matchtype *d = (Matchtype*) data;

	FREEMEMORY(space, d->pos);
	FREEMEMORY(space, data);
}

/*-------------------------------- getEntropy --------------------------------
 *    
 * calculates the entropy of a sequence, given probabilities.
 * 
 */

double getEntropy(void *space, Uint* sequence, Uint l, double* prob) {
	int i;
	double sum = 0;

	for (i = 0; i < l; i++) {
		sum += prob[sequence[i]] * log2(prob[sequence[i]]);
	}

	return sum;
}

/*local alignment*/
/*printf("max sw score: %f\n", occ[i-1].swscore);


 swres = swmatrix(space, queryseq->sequence, queryseq->length,
 s[occ[i-1].id]->sequence, s[occ[i-1].id]->length,
 -5, constscr, swscores);


 printf("max sw score: %d\n", swres[arraymax(swres,
 (queryseq->length+1)*(s[occ[i-1].id]->length+1))]);

 align = swgaplesstraceback(space, swres,
 queryseq->sequence, queryseq->length,
 s[occ[k].id]->sequence, s[occ[k].id]->length,
 //suffixscore,((imbissinfo*)info)->score
 -5,
 constscr, swscores,
 &alignsize);

 if (depictsw) {
 alignstr = printAlignment(space, align, alignsize,
 queryseq, s[occ[i-1].id], 80);
 printf("%s\n", alignstr);
 FREEMEMORY(space, alignstr);
 FREEMEMORY(space, align);
 }*/

Matchtype* selectBlastScoreSWconst(void *space, Matchtype *m, Uint k, IntSequence *a, IntSequence **s, void *info) {

	Uint l, i;
	int *swres;
	imbissinfo *imbiss;

	imbiss = (imbissinfo*) info;

	qsort(m, k, sizeof(Matchtype), cmp_blast);

	l = 0;
	for (i = k; i > 0 && l < 1000; i--) {
		if (m[i - 1].count >= imbiss->minimal_seed) {

			swres = swgapless(space, a->sequence, a->length, s[m[i - 1].id]->sequence, s[m[i - 1].id]->length, constscr, imbiss->swscores
			/*subscr, info*/
			);

			m[i - 1].swscore = swres[arraymax(swres, (a->length + 1) * (s[m[i - 1].id]->length + 1))];

			FREEMEMORY(space, swres);
		} else {
			m[i - 1].swscore = 0;
		}
		l++;
	}

	qsort(m, k, sizeof(Matchtype), cmp_swscore);
	return m;
}

Matchtype* selectScoreSWconst(void *space, Matchtype *m, Uint k, IntSequence *a, IntSequence **s, void *info) {

	Uint l, i;
	int *swres;
	imbissinfo *imbiss;

	imbiss = (imbissinfo*) info;

	qsort(m, k, sizeof(Matchtype), cmp_score);

	l = 0;
	for (i = k; i > 0 && l < 1000; i--) {
		if (m[i - 1].count >= imbiss->minimal_seed) {

			swres = swgapless(space, a->sequence, a->length, s[m[i - 1].id]->sequence, s[m[i - 1].id]->length, constscr, imbiss->swscores
			/*subscr, info*/
			);

			m[i - 1].swscore = swres[arraymax(swres, (a->length + 1) * (s[m[i - 1].id]->length + 1))];

			FREEMEMORY(space, swres);
		} else {
			m[i - 1].swscore = 0;
		}
		l++;
	}

	qsort(m, k, sizeof(Matchtype), cmp_swscore);
	return m;
}

Matchtype* selectSW(void *space, Matchtype *m, Uint k, IntSequence *a, IntSequence **s, void *info) {
	qsort(m, k, sizeof(Matchtype), cmp_swscore);
	return m;
}

Matchtype* selectBlastScore(void *space, Matchtype *m, Uint k, IntSequence *a, IntSequence **s, void* info) {
	qsort(m, k, sizeof(Matchtype), cmp_blast);
	return m;
}

Matchtype* selectScore(void *space, Matchtype *m, Uint k, IntSequence *a, IntSequence **s, void* info) {
	qsort(m, k, sizeof(Matchtype), cmp_score);
	return m;
}

double scorefilter(void *space, Matchtype *m, IntSequence *a, IntSequence *b, Uint *ptr, Uint len, Uint pos, void *info) {
	Uint l;
	double temp = 0;
	double sum = 0;
	imbissinfo *imbiss = (imbissinfo*) info;

	m->count++;
	m->pos = ALLOCMEMORY(space, m->pos, Uint, m->count);
	m->org = ALLOCMEMORY(space, m->org, Uint, m->count);
	m->pos[(m->count) - 1] = pos;
	m->org[(m->count) - 1] = pos;

	if (imbiss->score != NULL) {
		for (l = 0; l < len; l++) {
			temp = imbiss->score[*ptr];
			sum += temp;
			m->score += temp;
			ptr++;
		}
	}

	m->blast = m->blast > sum ? m->blast : sum;

	if (imbiss->consensus != NULL) {
		imbiss->consensus[pos] += (Uint) 1;
	}

	return sum > 0 ? sum : 0;
}

double swconstfilter(void *space, Matchtype *match, IntSequence *a, IntSequence *b, Uint *ptr, Uint len, Uint pos, void *info) {

	int *swres = NULL;
	imbissinfo *imbiss = NULL;
	double t;

	imbiss = (imbissinfo*) info;
	t = scorefilter(space, match, a, b, ptr, len, pos, info);

	if (match->count >= imbiss->minimal_seed) {
		swres = swgapless(space, a->sequence, a->length, b->sequence, b->length, constscr, imbiss->swscores);
		match->swscore = swres[arraymax(swres, (a->length + 1) * (b->length + 1))];

		FREEMEMORY(space, swres);
	}

	return t;
}

/*------------------------------ initMatchtype -------------------------------
 *    
 * init a Matchtype struct
 * 
 */

void initMatchtype(Matchtype *m, Uint id) {

	m->id = id;
	m->count = 0;
	m->pos = NULL;
	m->org = NULL;
	m->m = 0;
	m->blast = 0;
	m->score = 0;
	m->swscore = 0;

	return;
}

/*-------------------------------  rankSufmatch  ------------------------------
 *    
 * ranks matches 
 * given in an array of PairSint of length len. Sorting is done
 * by several calls to clib's qsort. For each item of the sorted
 * array a handler-function is invoked.
 * 
 */

void rankSufmatch(void *space, Suffixarray *suffix_array, PairSint *matches, Uint len, //
		Uint maxmatches, Uint substrlen, IntSequence **sequences, Uint noofseqs, //
		double (*filter)(void *, Matchtype *, IntSequence *, IntSequence *, Uint *, Uint, Uint, void *), // filter callback
		Matchtype* (*selector)(void *, Matchtype *, Uint, IntSequence *, IntSequence **, void *), // selector callback
		int (*handler)(void *, IntSequence *, Matchtype *, IntSequence **, Uint, Uint, void *), // handler callback
		IntSequence *sequence_a, void *info, double *scores, unsigned char depictsw)

{
	Uint *ptr = NULL;
	Matchtype *match = NULL;
	Matchtype *selected = NULL;
	int *hash_table = NULL;
	int i = 0, j = 0, k = 0;

	hash_table = ALLOCMEMORY(space, NULL, int, (noofseqs + 1));
	memset(hash_table, -1, sizeof(int) * (noofseqs + 1));

	for (i = 0; i < len; i++) {

		const PairSint * pair = &matches[i];
		if (pair->b < pair->a) {
			continue;
		}

		for (j = pair->a; j <= pair->b; j++) {
			const Uint index = getMultiSeqIndex(suffix_array->seq, suffix_array->suffixptr[suffix_array->suftab[j]]);
			const int r = hash_table[index];

			if (r == -1) {
				selected = ALLOCMEMORY(space, selected, Matchtype, k + 1);
				hash_table[index] = (&selected[k]) - (selected);
				initMatchtype(&selected[k], index);
				match = &selected[k];
				k++;
			} else {
				match = ((Matchtype*) (selected + r));
			}

			/*score the matches if no < maxmatches*/
			if ((pair->b - pair->a) >= maxmatches) {
				continue;
			}

			ptr = (suffix_array->suffixptr[suffix_array->suftab[j]]);
			IntSequence * sequence_b = (IntSequence *) sequences[match->id];

			massert((info != NULL), "Info can not be null");
			massert((sequence_a != NULL), "Sequence a can not be null");
			massert((sequence_b != NULL), "Sequence b can not be null");
			massert((match != NULL), "Match can not be null");
			massert((ptr != NULL), "Ptr can not be null");

			const int response = filter(space, match, sequence_a, sequence_b, ptr, substrlen, i, info);
			if (response == -1) {
				break;
			}

		}
	}

	selected = selector(space, selected, k, sequence_a, sequences, info);

	int length = 0;
	int response = 0;
	for (i = k; i > 0; i--) {
		Matchtype *match = (Matchtype *) &selected[i - 1];
		if ((response = handler(space, sequence_a, match, sequences, len, length, info))) {
			length++;
		}
		if (response == -1) {
			break;
		}
	}

	FREEMEMORY(space, hash_table);

	for (i = 0; i < k; i++) {
		FREEMEMORY(space, selected[i].pos);
		FREEMEMORY(space, selected[i].org);
	}

	FREEMEMORY(space, selected);
}

