#ifndef IMBISSBLAST_H
#define IMBISSBLAST_H

/*
 *
 *	imbissblast.h
 *  declarations for imbissblast.c to get
 *  blast params lambda and K
 * 
 *  @author Steve Hoffmann, shoffmann@zbh.uni-hamburg.de
 *  @company Center for Bioinformatics, Hamburg 
 *  @date 04/03/07 14:34:18 CEST  
 *
 */
#include "intsequence.h"
#include "falphabet.h"
#include "wurstimbiss.h"

typedef struct {
	
  	double *scr;
	double H;
	double K;
	double lambda;

} imbissblast;

void getimbissblast(void *space, IntSequence *query, IntSequence **seqs, 
	Uint noofseqs, FAlphabet *alphabet, imbissinfo *imbiss); 

#endif
