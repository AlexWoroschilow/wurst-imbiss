#ifndef WURST_IMBISS_H
#define WURST_IMBISS_H

/*
 *
 *	wurstimbiss.h
 *  basic declarations
 * 
 *  @author Steve Hoffmann, shoffmann@zbh.uni-hamburg.de,
 *  @author Alex Woroschilow alex.woroschilow@gmail.com
 *
 *  @company Center for Bioinformatics, Hamburg 
 *  @date 03/25/07 01:24:22 CET  
 *
 */

#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "memory.h"
#include "basic-types.h"
#include "stringutils.h"
#include "falphabet.h"
#include "intsequence.h"
#include "sufmatch.h"
#include "massert.h"
#include "score_mat.h"

typedef int (*imbissinfo_handler)(void *, IntSequence *, Matchtype *, IntSequence **, Uint, Uint, void *);
typedef double (*imbissinfo_filter)(void *, Matchtype *, IntSequence *, IntSequence *, Uint *, Uint, Uint, void *);
typedef Matchtype* (*imbissinfo_select)(void *, Matchtype *, Uint k, IntSequence *, IntSequence **, void *);
typedef IntSequence* (*imbissinfo_loader)(void *imbiss, void *space, char *filename);

typedef struct {

	double *score;
	double *scrf;
	double *sf;
	double *df;
	double *sub;
	char *reportfile;

	int *swscores;

	unsigned char wurst;
	unsigned char depictsw;

	Uint substrlen;
	Uint *sortind;

	FAlphabet *alphabet;
	stringset_t *query;
	Uint *consensus;
	struct sub_mat * matrix_substitition;

	double lambda;
	double H;
	double K;
	imbissinfo_handler handler;
	imbissinfo_loader loader;
	imbissinfo_filter filter;
	imbissinfo_select select;

	Uint maximal_hit;
	Uint maximal_match;
	Uint minimal_seed;
	Uint minimal_length;

	const char * file_configuration;
	char type[32];
	char logger[1024];
	char structlist[1024];
	char substitution[1024];
	char alphabetfile[1024];
	char sequencelist[1024];
	char binarypath[1024];
	char vectorpath[1024];

} imbissinfo;

#endif
