
/*
 * boyermoore.c
 * implementation of exact string matching
 * Boyer-Moore style for arbitrary alphabets
 *
 * @author Steve Hoffmann
 * @date Sat 2 Nov 2006
 *
 */

 #include <string.h>
 #include "basic-types.h"
 #include "mathematics.h"
 #include "memory.h"

 Uint getASCIIcode(void * alphabet, Uint asize, void *p, Uint pos) {
	char *pat = (char *)p;
	Uint ch = (Uint) pat[pos];

	return ch;
 }

 /* ------------------------- badChar --------------------------------
  *
  * returns an array for bad character heuristics
  * for char c bad[c] returns the last occurence of c in pattern p
  *
  */
 Sint * 
 badChar(void *space, void *alphabet, Uint asize, 
 		 void *p, Uint plen,
 		 Uint (*charpos) (void *, Uint, void *, Uint)) 
 {
	Uint i; 
	Sint *bad;
	
	bad = ALLOCMEMORY(space, NULL, Sint, asize);
	for(i=0; i < asize; i++) bad[i]=-1;
	for(i=0; i < plen; i++) 
		bad[charpos(alphabet,asize,p,i)]=i;
	
	return bad;
 }

 /* ------------------------- goodSuffix -----------------------------
  *
  * returns an array for good suffix heuristics
  * for each position i in pattern p of length l, s[i] indicates 
  * the suffix shifting distance. 
  *
  */
 Sint *
 goodSuffix(void *space, void *alphabet, Uint asize,
 			void *p, Uint l, 
			Uint (*charpos) (void *, Uint, void *, Uint)) 
 {
	Uint i,j=l+1;
	Sint *s, *f;
	
	f = ALLOCMEMORY(space, NULL, Sint, l+1);
	s = ALLOCMEMORY(space, NULL, Sint, l+1);
	memset(s, 0, (l+1)*sizeof(Sint));
	memset(f, 0, (l+1)*sizeof(Sint));
	
	f[l]=l+1;
	for(i=l; i > 0; i--) {
		while (j <= l && 
			charpos(alphabet, asize, p, i-1) !=
			charpos(alphabet, asize, p, j-1) ) {
				if(s[j]==0) {
					s[j]=j-i;
				}
				j = f[j];
		}
		f[i-1]=j-1;
		j--;
	}

	j = f[0];
	for(i=0; i <= l; i++) {
		if (s[i] == 0) s[i]=j;
		if (i==j) j = f[j];
	}

	FREEMEMORY(space, f);
	return s;
 }

 /* ------------------------- boyerMoore -----------------------------
  *
  * returns a vector v indicating all matches of pattern p and
  * template t of length tlen and plen respectively.
  *
  */
 vector_t *
 boyerMoore(void *space, void *alphabet, Uint asize, 
 			void *t, Uint tlen,
			void *p,  Uint plen,
			Uint (*charpos) (void *, Uint, void *, Uint)) 
 {
 	Sint i=0, j;
	Sint *bad, *suf;
	vector_t *res=NULL;
	
	bad = badChar(space, alphabet, asize, p, plen, charpos);
	suf = goodSuffix(space, alphabet, asize, p, plen, charpos);
	res = ALLOCMEMORY(space, NULL, vector_t, 1);
	INITVECTOR(res);
	
 	while (i<=tlen-plen) {
		j = plen-1;
		while(j>0 && 
				charpos(alphabet, asize, p, j) ==
				charpos(alphabet, asize, t, i+j)) j--;
				
		if(j==0 && 
			charpos(alphabet, asize, p, 0) ==
			charpos(alphabet, asize, t, i)) {
			
			APPENDVEC(space, res, i);	
			i+=suf[0];
		} else {
			i+=MAX(suf[j+1],(j-bad[charpos(alphabet, asize, t, i+j)]));
		}
	}
	
	FREEMEMORY(space, bad);
	FREEMEMORY(space, suf);
	return res;
 }
	

