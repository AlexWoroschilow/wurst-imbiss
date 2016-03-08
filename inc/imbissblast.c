/*
 *  imbissblast.c
 *  calculating blast statistics for imbiss
 *
 *  @author Steve Hoffmann
 *  @email shoffmann@zbh.uni-hamburg.de
 *  @date 04/03/07 14:28:01 CEST
 *  
 */

#include <stdio.h>
#include <stdlib.h>
#include "memory.h"
#include "basic-types.h"
#include "sort.h"
#include "imbissblast.h"
#include "falphabet.h"
#include "intsequence.h"
#include "mathematics.h"
#include "blaststat.h"

/*----------------------------- createsubmatrix ------------------------------
 *    
 * create a substitution matrix. Function accepts array of n scores
 * and a match matrix M^{nxn}
 * 
 */

double*
createsubmatrix(double *matrix, double *scr, Uint noofscores) {
	Uint i, j;

	for (i = 0; i < noofscores; i++) {
		for (j = 0; j < noofscores; j++) {
			MATRIX2D(matrix, noofscores, i, j)=
			MATRIX2D(matrix, noofscores, i, j) * scr[i];
			if(i==j) {
				MATRIX2D(matrix, noofscores, i, j) = scr[i];
			}
		}
	}

	return matrix;
}

int repsort(const void *a, const void *b) {
	double *first = (double *) a;
	double *second = (double *) b;

	if (*first < *second)
		return -1;
	if (*first == *second)
		return 0;
	if (*first > *second)
		return 1;

	return 0;

}

/*------------------------------ getimbissblast ------------------------------
 *    
 * getting lambda and K
 * 
 */

const char * getimbissblast(void *space, IntSequence *query, IntSequence **seqs, Uint noofseqs, FAlphabet *alphabet,
		imbissinfo *imbiss) {

	char * response;

	double *df, *sf, *scr;
	double avgsum = 0, inputscr = 0, lambda = 0, K = 0;
	Uint *sortind, i;
	evdparam *evd;

	/*frequency of query and database*/

	df = dbfreq(space, seqs, noofseqs, alphabet, 1);
	sf = seqfreq(space, query, alphabet);
	scr = logoddscr(space, df, sf, alphabet);

	for (i = 0; i < alphabet->domainsize; i++) {
		avgsum += df[i] * scr[i];
	}

	for (i = 0; i < query->length; i++) {
		inputscr += scr[query->sequence[i]];
	}

	sortind = quickSort(space, scr, alphabet->domainsize, cmp_dbl, NULL);

	evd = ALLOCMEMORY(space, NULL, evdparam, 1);
	evd->noofscores = alphabet->domainsize;
	evd->probs = df;
	evd->scores = scr;

	lambda = uniroot(0, 1, score_evd, 0.0000001, evd);
	FREEMEMORY(space, evd);

	K = relentropy(space, sortind, scr, alphabet->domainsize, df, lambda);
	if (K <= 0) {
		K = 1;
	}

	imbiss->score = scr;
	imbiss->H = 0;
	imbiss->K = K;
	imbiss->lambda = lambda;

	FREEMEMORY(space, sortind);
	FREEMEMORY(space, df);
	FREEMEMORY(space, sf);

	int strlen = snprintf(NULL, 0, "BLAST: E(score): %f, Inputscr: %f, Lambda: %19.16e, Check: %19.16e, K: %19.16e",
			avgsum, inputscr, lambda, checklambda(scr, alphabet->domainsize, df, avgsum, lambda), K);

	response = malloc((strlen + 1) * sizeof(*response));
	massert((response!=NULL), "Can not allocate memory for out string");

	snprintf(response, (strlen + 1), "BLAST: E(score): %f, Inputscr: %f, Lambda: %19.16e, Check: %19.16e, K: %19.16e",
			avgsum, inputscr, lambda, checklambda(scr, alphabet->domainsize, df, avgsum, lambda), K);

	return (const char *) response;
}

