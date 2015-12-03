
/*
 * encodeprobvec.c
 * encodes probability vectors to sequences
 *
 * @author Steve Hoffmann
 * @date Mon 27 Nov 2006
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

 IntSequence *encode_prob_vec(void *space, FAlphabet *alphabet, 
	 							struct prob_vec *pvec, Uint a, Uint b,
 								Uint (*assign)(void*, FAlphabet *, float *, 
									Uint, void *),
								void *info) 
 {
	Uint i;
	IntSequence *seq;

	if (b==0)
		b = pvec->n_pvec;

	seq = initSequence(space);	
	seq->sequence = ALLOCMEMORY(space, NULL, Uint, (b-a));
	seq->length = b-a;
	
	for (i=a; i < b; i++) {	
		seq->sequence[i-a]= assign(space, alphabet, pvec->mship[i], 
								pvec->n_class, info);
	}
	return seq;
 }


 Uint cantorchar(void *space, FAlphabet *alphabet, float *v, Uint size, 
	 				void *info) 
 {
	Uint ch, l, *sorted;
	vector_t a, *recoded = (vector_t*)info;

	INITVECTOR(&a);
	sorted = quickSort(space, v, size, cmp_flt, NULL);
	
	APPENDVEC(space, &a, sorted[size-1]+ALPHABET_OFFSET);
	if ( v[sorted[size-2]] > 0.0001)
		APPENDVEC(space, &a, sorted[size-2]+ALPHABET_OFFSET);
	else
		APPENDVEC(space, &a, sorted[size-1]+ALPHABET_OFFSET);
		
	ch = codeCantor(&a);	
	l = lookupChar(alphabet, ch);
	
	/*in case lookup for ch failed*/
	if (l > alphabet->mapsize) {
		/*notify via info*/
		APPENDVEC(space, recoded, VECTOR(&a,1));
		
		VECTOR(&a,1) = sorted[size-1]+ALPHABET_OFFSET;
		ch = codeCantor(&a);
		l = lookupChar(alphabet, ch);
		
	} else {
		APPENDVEC(space, recoded, 0);
	}

	FREEMEMORY(space, sorted);
	FREEMEMORY(space, a.elements);
	return alphabet->characters[l];
 }
 
