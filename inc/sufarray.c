
/*
 *  sufarray.c
 *  implementations for enhanced suffix arrays
 *  for large integer alphabets
 *
 *  @author Steve Hoffmann
 *  @email shoffmann@zbh.uni-hamburg.de
 *  @date 12/11/06 14:56:57 CET
 *  
 */
 
 #include <stdlib.h>
 #include <stdio.h>
 #include "basic-types.h"
 #include "memory.h"
 #include "mathematics.h"
 #include "sufarray.h"
 #include "intsequence.h"
 #include "falphabet.h"
 #include "sort.h"

/* ------------------------------ cmpIntSequence ----------------------------
 *    
 * function to compare IntSequences for mulitkey sort (sort.c)
 * 
 */
 
Uint
cmpIntSequence (Uint a, Uint b, Uint depth, void *data, void *info)
{
    Uint **s = (Uint**) data;	
	
	/*quick fix to meet end of multiintsequence criterion*/
	if(s[b][0] == 4001) {
		if (s[a][depth] == 4001) {
			return 0;
		}
		return 1;
	}
	
	/*real comparison*/
	if (s[a][depth] > s[b][depth]) return 1;
	if (s[a][depth] < s[b][depth]) return 2;
	
	return 0;
}


/*------------------------------- getSuffixPtr -------------------------------
 *    
 * returns an array of pointers for all "real" suffixes in mseq.
 * 
 */

Uint**
getSuffixPtr(void *space, MultiIntSeq *mseq, Uint delim, Uint sentinel) 
{
	Uint i,r,k=0;
	Uint **ptr;

	ptr = ALLOCMEMORY(space, NULL, Uint*, 
					(mseq->totallength)-(mseq->numofsequences)+1);
	
	for(i=0; i < mseq->totallength; i++) {
	  	if((r=mseq->sequences[i]) != delim){
			ptr[k]= &mseq->sequences[i];	  
			k++;
		}
	}
	
	return ptr;
}
  

/* ---------------------------- constructSufArr -----------------------------
 *    
 * constructs a suffix array from an (unsigned) integer sequence
 * should be working in O(n log n). It uses a multi-key quicksort.
 * 
 */

Suffixarray*
constructSufArr(void *space, IntSequence **s, Uint len, FAlphabet* alphabet)
{
 	Uint i, numofsuffixes,
		 *sorted, 
		 *inv_suftab,
		 **suffixptr;
	MultiIntSeq *mseq;
	Suffixarray *arr;
	
	mseq = concatIntSequences(space, s, len, 4000, 4001);
	numofsuffixes = (mseq->totallength - mseq->numofsequences)+1; 
	
	inv_suftab = ALLOCMEMORY(space, NULL, Uint , numofsuffixes);
	arr = ALLOCMEMORY(space, NULL, Suffixarray, 1);
	
	suffixptr = getSuffixPtr(space, mseq, 4000, 4001);
	sorted = quickSortMultikey (space, suffixptr, numofsuffixes, 
	cmpIntSequence, numofsuffixes-1, NULL);
		
	for (i=0; i < numofsuffixes; i++) {
		inv_suftab[sorted[i]]=i;
	}

	arr->seq = mseq;
	arr->numofsuffixes = numofsuffixes;
	arr->suffixptr = suffixptr;
	arr->suftab = sorted;
	arr->inv_suftab = inv_suftab;
	
	return arr;
}


/*------------------------------ destructSufArr ------------------------------
 *    
 * destruct a suffix array.
 * 
 */
void
destructSufArr (void *space, Suffixarray *arr)
{
    FREEMEMORY(space, arr->suftab);
	if (arr->lcptab != NULL)
	FREEMEMORY(space, arr->lcptab);
	if (arr->inv_suftab != NULL)
	FREEMEMORY(space, arr->inv_suftab);
	if (arr->suffixptr != NULL)
	FREEMEMORY(space, arr->suffixptr);
	if (arr->seq != NULL)
	destructMultiSeq(space, arr->seq);
	FREEMEMORY(space, arr);
	
	return ;
}

/*------------------------------ computeLcpTab -------------------------------
 *    
 * computes the lcp tab from suftab and inv_suftab in O(n).
 * 
 */
 
void
constructLcp (void *space, Suffixarray *arr)
{
  	Uint i, j, k;
	Uint s,t;
	int l=0;	
	
	arr->lcptab = ALLOCMEMORY(space, NULL, Uint, arr->numofsuffixes);
	
	for(i=0; i < arr->numofsuffixes; i++) {
	  	j = arr->inv_suftab[i];
		
		if (j > 0) {
		  	k = arr->suftab[j-1];
			s = arr->suffixptr[k]-arr->seq->sequences;
			t = arr->suffixptr[i]-arr->seq->sequences;
			
			l=l-1;
			if (l < 0) l=0;
						
			while ((s+l < arr->seq->totallength) && (t+l < arr->seq->totallength) && (arr->seq->sequences[s+l] == arr->seq->sequences[t+l])){ 
				l++;
			}
			arr->lcptab[j] = l;
		}
	}
	return;
}


/*-------------------------------- dumpSufArr --------------------------------
 *    
 * dumps a suffix array to a screen
 * 
 */
 
void
dumpSufArr (Suffixarray *arr)
{
  	Uint i;
	
  	for(i=0; i < arr->numofsuffixes; i++) {
		printf("%d \t %d \t %d \t %d \t %d \t %d\n", i, 
								arr->suftab[i], 
								arr->lcptab[i],
								arr->inv_suftab[i], 
								*arr->suffixptr[arr->suftab[i]],
								arr->seq->sequences[arr->suffixptr[arr->suftab[i]]-arr->seq->sequences]);
	}
	
	return;
}

void
dumplcps(Suffixarray *arr) {
	Uint i, j, s, t;


	for(i=0; i < arr->numofsuffixes; i++) {
		if (arr->lcptab[i] > 0) {
		  	s = arr->suffixptr[arr->suftab[i-1]]-arr->seq->sequences;
			t = arr->suffixptr[arr->suftab[i]]-arr->seq->sequences;
		    printf("lcp of suffix %d and %d has length %d\t:\n", i-1, i, arr->lcptab[i]);
			for(j=0; j <= arr->lcptab[i]; j++) printf(" %d ", arr->seq->sequences[s+j]);
			printf("\n");
			for(j=0; j <= arr->lcptab[i]; j++) printf(" %d ", arr->seq->sequences[t+j]);
			printf("\n");
		}
	}
}



