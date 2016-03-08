#ifndef IMSUBST_H
#define IMSUBST_H

/*
 *
 *	imsubst.h
 *  declarations to calculate substitution matrices
 * 
 *  @author Steve Hoffmann, shoffmann@zbh.uni-hamburg.de
 *  @company Center for Bioinformatics, Hamburg 
 *  @date 02/19/07 18:39:11 CET  
 *
 */

 #include "intsequence.h"
 #include "falphabet.h"


 double** submatrix (void *, double **, Uint, Uint);

 double** avgpvec(void *, FAlphabet *, double **,
	 			 struct prob_vec *, Uint, Uint, 
				 Uint(*assign)(void*, FAlphabet *, float *, Uint, void *),
				 void *);

#endif
