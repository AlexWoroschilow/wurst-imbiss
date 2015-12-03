
/*
 * falphabet.c
 * implmentations for a flexible alphabet
 *
 */

 #include <stdlib.h>
 #include "basic-types.h"
 #include "falphabet.h"
 #include "sort.h"
 #include "memory.h"

 Uint lookupChar(FAlphabet* alphabet, Uint mapped) {
	Uint ch;
	
	ch = binarySearch(alphabet->mapdomain, alphabet->mapsize, &mapped, cmp_int_bin, NULL);
	return ch;
 	
 }

 void destructAlphabet (void *space, FAlphabet *alphabet) {
	
   	FREEMEMORY(space, alphabet->characters);
	FREEMEMORY(space, alphabet->mapdomain);
 	FREEMEMORY(space, alphabet);
 }
