
/*
 *
 *	hash.c
 *  implementations for various hash functions
 * 
 *  @author Steve Hoffmann, shoffmann@zbh.uni-hamburg.de
 *  @company Center for Bioinformatics, Hamburg 
 *  @date 04/03/07 18:58:53 CEST  
 *
 */

#include "hash.h"

/*-------------------------------- knuthhash ---------------------------------
 *    
 * implementation of the *poor* knuth hash function
 * 
 */
 
unsigned long 
knuthhash(unsigned char *key, int len) {
	unsigned long hash;
	for (hash=len; len--;) {
		hash = ((hash << 5)^(hash>>27))^*key++;
	}
	return hash;
}


/*------------------------------- jenkinshash --------------------------------
 *    
 * implementation of the *widely appreciated* 
 * jenkins one-at-a-time hash function
 * 
 */
 
unsigned long
jenkinshash (unsigned char *key, int len)
{
     unsigned long hash = 0;
	 int i;
	
	 for (i = 0; i < len; i++) {
	   hash += key[i];
	   hash += (hash << 10);
	   hash ^= (hash >> 6);
	}
	hash += (hash << 3);
	hash ^= (hash >> 11);
	hash += (hash << 15);
	
	return hash;
}
