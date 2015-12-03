
/*
 *  mm.c
 *  implementations for searches manber myers style
 *  on enhanced suffix arrays
 *
 *  @author Steve Hoffmann
 *  @email shoffmann@zbh.uni-hamburg.de
 *  @date 12/22/06 19:13:27 CET
 *  
 */

 #include <stdlib.h>
 #include <stdio.h>
 #include "basic-types.h"
 #include "memory.h"
 #include "mathematics.h"
 #include "sufarray.h"
 #include "intsequence.h"
 #include "mm.h"

/*---------------------------------- mmleft ----------------------------------
 *    
 * part of the manber-myers pattern search algorithm
 * 
 */
 
int
mmleft(Suffixarray *arr, Uint *pattern, Uint len, int h, int l, int r) 
{
  	PairSint lbound, rbound, ibound;
	int mid;

	lbound = mmcompare (arr, pattern, len, l, h);
	if(lbound.b <= 0) {
		return l;
	}
	rbound = mmcompare (arr, pattern, len, r, h);
	if(rbound.b > 0) {
		return r+1;
	}	
	while (r > l+1) {
		mid = (l+r)/2;
		ibound = mmcompare (arr, pattern, len, mid, MIN(lbound.a, rbound.a));
		if(ibound.b <= 0) {
			rbound.a = ibound.a;
			r = mid;
		} else {
			lbound.a = ibound.a;
			l = mid;
		}
	}
	return r;
}

/*--------------------------------- mmright ----------------------------------
 *    
 * part of the manber-myers pattern search algorithm
 * 
 */
 
int
mmright (Suffixarray *arr, Uint *pattern, Uint len, int h, int l, int r)
{
  	PairSint lbound, rbound, ibound;
	int mid;
	
	lbound = mmcompare(arr, pattern, len, l, h);
	if (lbound.b < 0) {
		return -1; 
	}
	rbound = mmcompare(arr, pattern, len, r, h);
	if (rbound.b >= 0) {
		return r;
	}
	
	while (r > l+1) {
		mid = (l+r)/2;
		ibound = mmcompare(arr, pattern, len, mid, MIN(lbound.a, rbound.a));
		if (ibound.b >= 0 ) {
			lbound.a = ibound.a;
			l = mid;
		} else {
			rbound.a = ibound.a;
			r = mid;
		}
	}	
	return l;
}

/*-------------------------------- mmcompare ---------------------------------
 *    
 * manber 'n' myers compare for suffixarrays
 * 
 */
 
PairSint
mmcompare (Suffixarray *arr, Uint *pattern, Uint len, int idx, int start)
{
    Uint *sufptr;
	int t, margin;
	PairSint res;
	
	sufptr = arr->suffixptr[arr->suftab[idx]];
	t = start;
						/*length = totallength - relativeposition*/
	margin = MIN((len),(arr->seq->totallength-(sufptr-arr->seq->sequences)));
	
	while (t < margin) {
		if(pattern[t] < sufptr[t]) {
			res.a=t;
			res.b=-1;
			
			return res;
		} else {
			if(pattern[t] > sufptr[t]) {
				res.a=t;
				res.b=1;
				
				return res;
			} else {
				t++;
			}
		}
	}
	
	if (t == len) {
		  res.a=t;
		  res.b=0;
		  return res;
	} 
	
	res.a=t;
	res.b=-1;
	
	return res;
}


/*--------------------------------- mmsearch ---------------------------------
 *    
 * search manber'n'myers style
 * 
 */
 
PairSint
mmsearch (Suffixarray *arr, Uint *pattern, Uint len, Uint h, Uint l, Uint r)
{
    int p, q;
	PairSint ps;

	p = mmleft(arr, pattern, len, h, l, r);
	q = mmright(arr, pattern, len, h, l, r);

	if (p <= q) {
	    ps.a = p;
		ps.b = q;
		
		return ps;
	}
	
	ps.a = 1;
	ps.b = 0;
	
	return ps;
}


