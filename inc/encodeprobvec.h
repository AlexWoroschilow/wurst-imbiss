#ifndef ENCODE_PROB_VEC_H
#define ENCODE_PROB_VEC_H

/*
 * encodeprobvec.h
 * declarations for encoding probability sequences
 *
 * @author Steve Hoffmann
 * @date Tue 28 Nov 2006
 *
 */

 #include "intsequence.h"
 #include "falphabet.h"

 #define ALPHABET_OFFSET 1

 IntSequence *encode_prob_vec(void *, FAlphabet *, struct prob_vec *, 
 								Uint, Uint, 
 								Uint (*assign)(void *, FAlphabet *, float *, Uint, void *),
								void *);
	
 Uint cantorchar(void *, FAlphabet *, float *, Uint, void*);

 #endif
