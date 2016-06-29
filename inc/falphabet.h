 #ifndef FALPHABET_H
 #define FALPHABET_H

/*
 * alphabet.h
 * declarations for a flexible alphabet
 *
 */

 #include "basic-types.h"

 typedef struct {

	Uint *characters,
		 *mapdomain;
	
	Uint domainsize,
		 mapsize,

		 mappedwildcards,
		 undefsymbol,
		 *symbolmap;

 } FAlphabet;


 /*from mapdomain to character*/
 Uint lookupChar(FAlphabet *, Uint);	
 void destructAlphabet(void *space, FAlphabet *);
 

#endif
