#ifndef SUF_ARRAY_H
#define SUF_ARRAY_H

/*
 *
 *	sufarray.h
 *  declarations for enhanced suffix arrays
 *  for large integer alphabets
 * 
 *  @author Steve Hoffmann, shoffmann@zbh.uni-hamburg.de
 *  @company Center for Bioinformatics, Hamburg 
 *  @date 12/10/06 22:01:33 CET  
 *
 */

 #include "basic-types.h"
 #include "falphabet.h"
 #include "intsequence.h"
 #include "multiseq.h"


 typedef struct {
	MultiIntSeq  *seq;
	 
	Uint		numofsuffixes;
	Uint		**suffixptr;
	Uint 		*suftab;
	Uint		*inv_suftab;
 	
    Uint 		*bwttab; 	/*burrows-wheeler array*/
	Uint		*lcptab;
	/*
	 *  alternative: Abouelhoda et al.
		char 		*lcptab;    nB to store lcp values < 255
		PairUint	*llvtab;    array of 8B to store lcp val >=255
	*/
 	
	char		*chldtab;   /*a child table*/
 	Uint		*bcktab;    /*the bucket container*/
 
 } Suffixarray;

Suffixarray* constructSufArr(void *, IntSequence **, Uint, FAlphabet *);
void dumplcps (Suffixarray *arr);
void constructLcp (void *, Suffixarray *);
void destructSufArr (void *, Suffixarray *); 
PairSint testPairSint();

#endif

