
#ifndef WURST_IMBISS_H
#define WURST_IMBISS_H

/*
 *
 *	wurstimbiss.h
 *  basic declarations
 * 
 *  @author Steve Hoffmann, shoffmann@zbh.uni-hamburg.de
 *  @company Center for Bioinformatics, Hamburg 
 *  @date 03/25/07 01:24:22 CET  
 *
 */

#include "basic-types.h"
#include "stringutils.h"
#include "falphabet.h"

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
	Uint noofhits;	
	Uint minseeds;
	Uint *sortind;
	
	FAlphabet *alphabet;
	stringset_t *query;
	Uint *consensus;

	double lambda;
	double H;
	double K;

} imbissinfo;



#endif
