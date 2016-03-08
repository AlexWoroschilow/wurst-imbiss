
/*
 *  multiseq.c
 *  some functions to handle multiseqs
 *
 *  @author Steve Hoffmann
 *  @email shoffmann@zbh.uni-hamburg.de
 *  @date 12/15/06 11:42:53 CET
 *  
 */

 #include <stdio.h>
 #include <stdlib.h>
 #include "basic-types.h"
 #include "memory.h"
 #include "intsequence.h"
 #include "sort.h"
 #include "multiseq.h"
 
/*---------------------------- concatIntSequences ----------------------------
 *    
 * concatenates IntSequences using a given Uint delimiter
 * and stores them in a MultiIntSeq container.
 * 
 */
 
MultiIntSeq *
concatIntSequences (void *space, IntSequence **s, Uint len, 
					Uint delim, Uint sentinel)
{
    Uint *buf=NULL, i, j, k=0, 
		 totallength=0, 
		 *markpos;
	MultiIntSeq *mseq;

	mseq = ALLOCMEMORY(space, NULL, MultiIntSeq, 1);
	markpos = ALLOCMEMORY(space, NULL, Uint, len);

		
	for(i=0; i < len; i++) {
	  	totallength += (s[i]->length+1);
		buf = ALLOCMEMORY(space, buf, Uint, totallength);
		if (buf==NULL) printf("allocation failed: exiting\n");

		for(j=0; j < s[i]->length; j++) {
			buf[k] = s[i]->sequence[j];
			k++;
		}
		/*separate sequences or finalize*/
		if (i == (len-1)) {
		  buf[k] = sentinel;
		  markpos[i] = k;
		  k++;
		} else {
		  buf[k] = delim;
		  markpos[i] = k;
		  k++;
		}
	}
	mseq->totallength = totallength;
	mseq->numofsequences = len;
	mseq->sequences = buf;
	mseq->markpos = markpos;
	
	return mseq;
}


/*----------------------------- destructMultiSeq -----------------------------
 *    
 * destructs a MultiSeq structure
 * 
 */

void
destructMultiSeq (void *space, MultiIntSeq *mseq)
{
    
	FREEMEMORY(space, mseq->sequences);
	FREEMEMORY(space, mseq->markpos);
	FREEMEMORY(space, mseq);
	return ;
}


/*------------------------------- cmp_markpos --------------------------------
 *    
 * compare function for getMultiSeqIndex
 * 
 */
 
Uint
cmp_markpos (Uint a, void *data, void *key, void *info)
{
    Uint *d = (Uint*) data;
	Uint *k = (Uint*) key;
	
	if (d[a] > *k) {
		if (a > 0) 
		{
			if (d[a-1] < *k) 
			{
				return 0;
			} 
			else 
			{
				return 1;
			}
		} 
		else 
		{
			return 0;
		}
	}
	
    if (d[a] < *k) return 2;

	return 0;
}

/*-------------------------- getMultiSeqDescription --------------------------
 *    
 * returns the index of a sequence in multiseq addressed by a pointer
 * 
 */
 
Uint
getMultiSeqIndex (MultiIntSeq *mseq, Uint *ptr)
{	
	Uint pos, i;
	
	pos = (ptr - mseq->sequences);
 
	i=binarySearch(mseq->markpos, mseq->numofsequences, &pos, 
			       cmp_markpos, NULL);

	return i;
}

/*---------------------------- getMultiSeqRelPos -----------------------------
 *    
 * returns the relative position of a pointer to mulitseq
 * with respect to the addressed sequence.
 * 
 */
 
Uint
getMultiSeqRelPos (MultiIntSeq *mseq, Uint *ptr)
{
	return 0;
}


/*------------------------------- dumpMultiSeq -------------------------------
 *    
 * dumps a multiseq to the screen
 * 
 */

void
dumpMultiSeq (MultiIntSeq *mseq)
{
  	Uint i;

	for(i=0; i < mseq->totallength; i++) {
		printf("%d-", mseq->sequences[i]);	
	}

	printf("\n");
	return ;
}

