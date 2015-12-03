
/*
 *  depictseqs.c
 *  simple implementations to visualize a sequence
 *  and its properties.
 *
 *  @author Steve Hoffmann
 *  @email shoffmann@zbh.uni-hamburg.de
 *  @date 12/25/06 15:53:21 CET
 *  
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "memory.h"
#include "depictseqs.h"
#include "basic-types.h"


/*------------------------------ depictSequence ------------------------------
 *    
 * depicts and annotates a (virtual) sequence.
 * 
 */
 
char*
depictSequence (void *space, Uint seqlen, Uint maplen, 
				Uint *depict, Uint len, char ch)
{
  	Uint i, pos;
	double scale;
	char *seq;
	
	if (maplen > seqlen) maplen = seqlen;
	scale = ((float) maplen)/((float)seqlen);
	
	seq = ALLOCMEMORY(space, NULL, char, maplen+1);
	memset(seq, '-', maplen);
	seq[maplen]='\0';

	for (i=0; i < len; i++) {
		pos = depict[i]*scale;
		seq[pos]=ch;
	}
	
	return seq;
}


