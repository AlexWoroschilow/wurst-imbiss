
/*
 *  dpalign.c
 *  basic routines for dynamic programming
 *
 *  @author Steve Hoffmann
 *  @email shoffmann@zbh.uni-hamburg.de
 *  @date 02/15/07 17:46:54 CET
 *  
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "basic-types.h"
#include "mathematics.h"
#include "dpalign.h"

/*---------------------------------- edist -----------------------------------
 *    
 * evaluation of the edit distance in O(n) space
 * the function accepts two sequences of symtype and
 * an m-dimensional substitution matrix 
 * 
 */

Uint edist(void *space, symtype *sa, Uint lena, symtype *sb, Uint lenb, 
				Uint indel, Uint *sub, Uint m) {

	Uint minlen, maxlen;
	Uint i, j, r1=0, r2=0;
	Uint pen;
	Uint* col;
	symtype *min;
	symtype *max;

	minlen = (MIN(lena, lenb)) + 1;
	maxlen = (MAX(lena, lenb)) + 1;
	col = ALLOCMEMORY(space, NULL, Uint, minlen);
	
	if(minlen-1 == lena) {
		min = sa;
		max = sb;
	} else {
		min = sb;
		max = sa;
	}
	
	for(i=0; i < maxlen; i++) {
		for(j=0; j < minlen; j++) {
			if (i==0) {
			  col[j]=j; 
			} else {
				if (j==0) {
					r1 = col[j]++; 
				} else {
					r2 = col[j];
					if (sub == NULL) {
						pen = (min[j-1]==max[i-1]) ? 0 : 1;
					} else {
						pen = (min[j-1]==max[i-1]) ? 0 : 
					  		MATRIX2D(sub, m, min[j-1], max[i-1]);
					}
					col[j] = MIN ((col[j-1]+indel),  
						     (MIN ((r1 + pen), (col[j]+indel))));
					r1 = r2;	
				}
			}
			printf("%d\n",col[j]);
		}
		printf("--\n");
	}
	
	r1 = col[minlen-1];
	FREEMEMORY(space, col);
	return r1;
}


/*--------------------------------- constscr ----------------------------------
 *    
 * a function that assigns constant scores for matches and mismatches
 * given in info[0] and info[1], respectively.
 * 
 */
 
int
constscr (symtype a, symtype b, void *info)
{	
  	int* scores;
	
  	scores = (int*) info;
	if(a == b) return scores[0];
	
  	return scores[1];
}


/*--------------------------------- swgapless ----------------------------------
 *    
 * smith-waterman local similarity alignment w/o gaps
 * returns a matrix of size (m+1)*(n+1) where m is length of given sequence a
 * and n the length of sequence b, respectively. Function expects
 * a function to calculate a substitution score
 * 
 */
 
int*
swgapless (void *space, symtype *a, Uint m, symtype *b, Uint n,
			Sint (*sub)(symtype, symtype, void *), void *nfo)
{
  	int i, j, cols, rows, size;
	int *L;
	
	rows = m+1;
	cols = n+1;

	size = rows*cols;
	L = ALLOCMEMORY(space, NULL, int, size);
	L = memset(L, 0, sizeof(int)*size);
	
	for(i=1; i < m+1; i++) {
		for(j=1; j < n+1; j++) {
		    
			MATRIX2D(L, cols, i, j) =
			     MAX(0,
				   MATRIX2D(L, cols, (i-1), (j-1)) + sub(a[i-1], b[j-1], nfo)
				);
		}
	}
	
	return L;
}





/*--------------------------------- swalign ----------------------------------
 *    
 * smith-waterman local similarity alignment
 * returns a matrix of size (m+1)*(n+1) where m is length of given sequence a
 * and n the length of sequence b, respectively. Function expects
 * a function to calculate a substitution score
 * 
 */
 
int*
swmatrix (void *space, symtype *a, Uint m, symtype *b, Uint n, int indel,
			Sint (*sub)(symtype, symtype, void *), void *nfo)
{
  	int i, j, cols, rows, size;
	int *L;
	
	rows = m+1;
	cols = n+1;

	size = rows*cols;
	L = ALLOCMEMORY(space, NULL, int, size);
	L = memset(L, 0, sizeof(int)*size);
	
	for(i=1; i < m+1; i++) {
		for(j=1; j < n+1; j++) {
		    
			MATRIX2D(L, cols, i, j) = 
			  MAX4(0,
				   MATRIX2D(L, cols, (i-1), j) + indel ,	
				   MATRIX2D(L, cols, i, (j-1)) + indel , 
				   MATRIX2D(L, cols, (i-1), (j-1)) + sub(a[i-1], b[j-1], nfo)
			  );
		}
	}
	
	return L;
}


/*------------------------------- swtraceback --------------------------------
 *    
 * traceback to find optimal local alignment path
 * 
 */
 
int*
swtraceback (void *space, int *M,  
			 symtype *a, Uint m, symtype *b, Uint n, int indel,
			 Sint (*sub)(symtype, symtype, void *), void *nfo, int* alignsize)
{
	Uint i, j, ncol, cur, start;
	int *align = NULL;
	
	*alignsize = 0;
	ncol = (n+1);
	start = arraymax(M, (m+1)*ncol);
	i = start / ncol;
	j = start % ncol;
	
	while(i > 0 && j > 0) {

 	  	cur = MATRIX2D(M, ncol, i, j);
		if (cur==0)
			return align;
		

		align = ALLOCMEMORY(space, align, int, *alignsize+2);
		align[*alignsize]=i;
		align[*alignsize+1]=j;
		*alignsize +=2;
	
		
		if (MATRIX2D(M, ncol, i-1, j) + indel == cur){
			i--;	
		} else {
			 if (MATRIX2D(M, ncol, i, j-1) + indel == cur) {
			  	j--;
			 } else {
			 	if (MATRIX2D(M, ncol, i-1, j-1)+sub(a[i-1], b[j-1], nfo) 
						== cur){
			  		i--; j--;
			  	}
			 }
		 }
	 }

	return align;
}

/*------------------------------- swgaplesstraceback -----------------------------
 *    
 * traceback to find optimal local alignment path
 * 
 */
 
int*
swgaplesstraceback (void *space, int *M,  
			 symtype *a, Uint m, symtype *b, Uint n, 
			 Sint (*sub)(symtype, symtype, void *), void *nfo, int* alignsize)
{
	Uint i, j, ncol, cur, start;
	int *align = NULL;
	
	*alignsize = 0;
	ncol = (n+1);
	start = arraymax(M, (m+1)*ncol);
	i = start / ncol;
	j = start % ncol;
	
	while(i > 0 && j > 0) {

 	  	cur = MATRIX2D(M, ncol, i, j);
		if (cur==0)
			return align;

		align = ALLOCMEMORY(space, align, int, *alignsize+2);
		align[*alignsize]=i;
		align[*alignsize+1]=j;
		*alignsize +=2;

		i--; j--;

	 }

	return align;
}



