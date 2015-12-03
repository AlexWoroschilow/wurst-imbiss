#ifndef BLAST_STAT_H
#define BLAST_STAT_H

/*
 *
 *	blaststat.h
 *  declarations for blast statistics
 * 
 *  @author Steve Hoffmann, shoffmann@zbh.uni-hamburg.de
 *  @company Center for Bioinformatics, Hamburg 
 *  @date 03/23/07 12:19:21 CET  
 *
 */

#include "falphabet.h"
#include "intsequence.h"

typedef struct {

  	double *probs;
	double *scores;
	Uint noofscores;

} evdparam;


double score_evd(double, void*);


double* logoddscr (void*, double*, double *, FAlphabet *);

double*
dbfreq (void *, IntSequence **, Uint, FAlphabet *,
		double);

double *
seqfreq (void *, IntSequence *, FAlphabet *);
 
double *
scorefreq ( void *, double *, Uint, FAlphabet *, 
			double *, double *);

double 
checklambda (double *, Uint, double*, double, double);

double
relentropy (void*, Uint *, double *, Uint, double *, double);

double
lambda (Uint *, double *, Uint, double *, double, double, double);

double
entropyK (Uint *, double *, Uint, double, double *, double, double);


#endif
