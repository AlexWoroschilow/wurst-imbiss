
/*
 *  blaststat.c
 *  an implementation of blast statistics
 *  completely re-written: no fancy shit implementation
 *  
 *  @author Steve Hoffmann
 *  @email shoffmann@zbh.uni-hamburg.de
 *  @date 02/26/07 21:18:34 CET
 *  
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "mathematics.h"
#include "memory.h"
#include "intsequence.h"
#include "falphabet.h"
#include "basic-types.h"
#include "blaststat.h"

#define BLAST_SUM_LIMIT 0.001 /*speedy*/
#define BLAST_ITER_MAX  100

/*--------------------------------  seqfreq  ---------------------------------
 *    
 * returns the frequency of query characters in alphabet A
 * Function accepts a query sequence and the corresponding
 * alphabet. Frequencies are returned in a vector.
 * 
 */
 
double *
seqfreq (void *space, IntSequence *seq, FAlphabet *alpha)
{
    Uint l, i;
	double *f, sum =0; 
	  
	l = alpha->domainsize;
	f = ALLOCMEMORY(space, NULL, double, l+1);
	memset(f, 0, sizeof(double)*(l+1));

	for(i=0; i < seq->length; i++) {
		f[seq->sequence[i]]++;
		sum++;
	}

	for(i=0; i < l; i++) {
		f[i]/=sum;
	}
	
	return f;
}


/*---------------------------------- dbfreq ----------------------------------
 *    
 * returns the frequencies in a set of sequences for each char in A 
 * ResFreqNormalize
 */
 
double*
dbfreq (void *space, IntSequence **set, Uint setsize, FAlphabet *alpha,
		double norm)
{
  	Uint i, j, l, sum=0;
	double *f;
	
	l = alpha->domainsize;
	f = ALLOCMEMORY(space, NULL, double, l+1);
	memset(f, 0, sizeof(double)*(l+1));
	
	for(i=0; i < l; i++) {
		f[i]=0;
	}

	for(i=0; i < setsize; i++) {
		for(j=0; j < set[i]->length; j++) {
			f[set[i]->sequence[j]]++;
			sum++;
		}
	}

	for(i=0; i < l; i++) {
		f[i] /= sum;
		f[i] *= norm;
	}

	
	return f;
}


#define BASIC_PROBABILITY 0.000001

/*-------------------------------- logoddscr ---------------------------------
 *    
 * calculations of log-odds scores (Karlin Altschul PNAS 89, example iii)
 * 
 */

double*
logoddscr (void* space, double* dbf, double *sf, FAlphabet *alpha)
{
	 Uint i;
	 double temp;
	 double *logodd=NULL;
 
	 logodd = ALLOCMEMORY(space, NULL, double, alpha->domainsize+1);
	 memset(logodd, 0, sizeof(double)*(alpha->domainsize+1));

	 for(i=0; i < alpha->domainsize; i++) {
		if (dbf[i] > 0) {
		    temp = (sf[i]+0.0000000001)/dbf[i];
			logodd[i] = log(temp)*10;	
		} else {
		  	logodd[i] = (double)0;
		}
	 }
	 
	return logodd;
}


void writescores(double *scr, Uint l, char* fname) {
	FILE *fp;
	Uint i;
	
	fp = fopen(fname, "w");
	if (fp == NULL) {
	    fprintf(stderr, "error: unable to open file.\n");
		return;
	}
	
	for(i=0; i < l; i++) {
		fprintf(fp, "%f\n", scr[i]);
	}

	fclose(fp);
}

/*-------------------------------- scorefreq ---------------------------------
 *    
 * calculates the frequencies of all given integer scores. 
 * Function accepts an array of scores and arrays containing
 * query and database frequencies for each char in alphabet A.
 * returns an array of frequencies. Additionally avg is filled.
 * 
 */
 
double *
scorefreq ( void *space, double *scr, Uint noofscr, FAlphabet *alpha, 
			double *sf, double *dbf)
{
	Uint i;
	double sum = 0;	
	double *f;

	f = ALLOCMEMORY(space, NULL, double, noofscr);
	memset(f, 0, sizeof(double)*noofscr);
	
	for(i=0; i < noofscr; i++) {
		f[i] = sf[i]*dbf[i];
		sum += sf[i]*dbf[i];
	}
	
	
	if (sum > 0.0001) {
		for(i=0; i < noofscr; i++) {
			f[i] /= sum;
		}
	}

  	return f;
}


double
checklambda (double *scr, Uint noofscores, double* sf, double avg, double e_lamb){
	Uint i;
    double val=0;
	
  	for(i=0; i < noofscores; i++) {
		val += sf[i]*exp((e_lamb*scr[i])); 
	}

	return val;
} 
		


/*-------------------------------- score_evd ---------------------------------
 *    
 * extreme value distribution for given scores (supplied in info)
 * 
 */
 
double
score_evd (double lambda, void* info)
{
  	Uint i;
	double sum=0;
	evdparam *p =(evdparam*) info;
	
	for(i=0; i < p->noofscores; i++) {
		sum += p->probs[i]*exp((double) lambda*p->scores[i]);
	}
	
	return (sum-1.0);
}

/*-------------------------------- relentropy --------------------------------
 *    
 * calculates the relative entropy of probabilities and target
 * frequencies. Subsequently, the parameter K is calculated.
 *
 * basics: 
 * scr is sorted - the observed minimum score is in scr[0]
 * probs contains the probabilities of scrs - eg. q[scr[0]]
 * 
 */
 
double
relentropy (void *space, Uint *sort, double *scr, Uint noofscr, double *q, double lambda)
{
  	double sum=0, 
		   av,
		   sumval, K, Ktmp;
	double *p, *P, *ptrP, *ptr1, *ptr2, *ptr1e, *resc, *scr2;
	int i, j, first, last, hi, lo, range=0;
	int fix=0;
	
  	if (lambda < 0)
	  return -1;

        /*sum the probabilities*/	
	for (i=0; i < noofscr; i++) {
		sum += q[i];
	}	

	/*normalize probabilites*/
	p = ALLOCMEMORY(space, NULL, double, noofscr);
	for (i=0; i < noofscr; i++) {
		p[i] = q[i] / sum;
	}

	/*resort*/
	resc  = ALLOCMEMORY(space, NULL, double, noofscr);
	scr2  = ALLOCMEMORY(space, NULL, double, noofscr);
	for(i=0; i < noofscr; i++) {
		resc[i] = p[i];
		scr2[i] = scr[sort[i]];
	}
	FREEMEMORY(space, p);
	p = resc;

	/*get average*/
	for (av=0.0, i=0; i < noofscr; i++) {
		av += p[i]*(scr2[i]+1)*(lambda*(scr2[i]+1));
	}

	range = ((int)scr2[noofscr-1]) - ((int) scr2[0]) ;
	
	/*get P*/
	fix += (int)1500*(int)range;

//	P = malloc(sizeof(double)* (size_t)fix);
	P = ALLOCMEMORY(NULL, NULL, double, fix);


	
	sumval=0;
	lo=0;
	hi=0;
	
	/*fancy sum /frac{1}{k=j} E[e^{\lambda*S_k}; S_k <0]+Prob(S_k>=0) }*/
	/*			^frac here		^calculated at the end of the loop     */
	
	for (*P = 1.0, sum = 1.0, j = 1;
			j <= 1500 && sum > 0.0000001; sumval += sum /= j++) {
	
	  first = last = range;
	  for (ptrP = P + ((hi += scr2[noofscr-1]) - (lo += scr2[0]));
		  ptrP >= P; *ptrP-- = sum) {
	  	ptr1   = ptrP - first;	
		ptr1e  = ptrP - last;
		ptr2   = p + first;
		
		for(sum = 0.0; ptr1 >= ptr1e;) {
		  sum += *ptr1-- * *ptr2++;
		}
		if (first) {
			--first;
		}
		if (((int) (ptrP-P)) <= range) {
			--last;
		}
	  }
	  
	  /*get \sum E[e^{\lambda*S_{k}}, S_k < 0]*/
	  for(sum=0.0, i = lo; scr2[i] < 0; ++i) {
	  	sum += *++ptrP * exp(lambda*scr2[i]);
	  }
	  for(; i <= hi; ++i) {
	  	sum += *++ptrP;
	  } 
	}

	Ktmp = (double) (exp(-2 * sumval));
	K = Ktmp / (av * (1.0 - exp(-lambda)));

	FREEMEMORY(space, p);
	FREEMEMORY(space, P);
	return K;
}
 
