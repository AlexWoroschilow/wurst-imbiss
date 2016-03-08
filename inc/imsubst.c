
/*
 *  imbsubst.c
 *  calculate the imbiss substitution matrix
 *
 *  @author Steve Hoffmann
 *  @email shoffmann@zbh.uni-hamburg.de
 *  @date 02/19/07 16:27:50 CET
 *  
 */

 #include <stdlib.h>
 #include <stdio.h>
 #include "memory.h"
 #include "basic-types.h"
 #include "prob_vec.h"
 #include "falphabet.h"
 #include "sort.h"
 #include "intsequence.h"
 #include "cantor.h"
 #include "encodeprobvec.h"
 #include "mathematics.h"
 #include "createalphabet.h"	



/*--------------------------------- avgpvec ----------------------------------
 *    
 * assign pvecs to character and sum up to get avg pvec 
 * 
 */

 double** 
 avgpvec (void *space, FAlphabet *alphabet, double** avg, 
	 	 struct prob_vec *pvec, Uint a, Uint b,
		 Uint (*assign)(void *, FAlphabet *, float *, Uint, void*),
		 void *info){
	
   Uint i, j, ndx;

	if (b==0)
		b = pvec->n_pvec;	
 	
 	for (i=a; i < b; i++) {	
	    ndx = assign(space, alphabet, pvec->mship[i], pvec->n_class, info);
	    for(j=0; j < pvec->n_class; j++) {
			avg[ndx][j] += pvec->mship[i][j];
		}
		avg[ndx][pvec->n_class]=avg[ndx][pvec->n_class]+1;
	}
	
	return avg;
 }


/*-------------------------------- submatrix ---------------------------------
 *    
 * return pairwise substitution matrix
 * 
 */
 
double**
submatrix (void *space, double **pvec, Uint length, Uint classes)
{
    double **m;
	Uint i,j,k;
	
	m = ALLOCMEMORY(space, NULL, double*, length+1);
	for (i=0; i < length+1; i++) {
		m[i]=ALLOCMEMORY(space, NULL, double, length+1);
		for(j=0; j < length+1; j++) {
			m[i][j]=0;
		}
	}

	for(i=0; i < length; i++) {
		for(j=0; j < length; j++) {
			for(k=0; k < classes; k++) {
				m[i][j]+=pvec[i][k]*pvec[j][k];
				m[i][length] += pvec[i][k]*pvec[j][k];
			}
		}
	}
	
	return m;
}

