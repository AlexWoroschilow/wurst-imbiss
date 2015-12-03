
/*
 * createalphabet.c
 * implementing functions to create and manipulate 
 * flexible alphabets
 *
 * @author Steve Hoffmann
 * @date Sat 25 Nov 2006
 *
 */

 #include <stdlib.h>
 #include <stdio.h>
 #include "basic-types.h"
 #include "stringutils.h"
 #include "fileio.h"
 #include "sort.h"
 #include "falphabet.h"
 #include "memory.h"

 FAlphabet* loadCSValphabet(void *space, char* filename) {
	Uint size, i, mapto;
	char* trimmed;
	FAlphabet *alphabet;
	stringset_t **csv;
	
	csv=readcsv(space, filename,",", &size);
	 
	alphabet = ALLOCMEMORY(space, NULL, FAlphabet, 1);
	alphabet->characters = ALLOCMEMORY(space, NULL, Uint, size-1);
	alphabet->mapdomain = ALLOCMEMORY(space, NULL, Uint, size-1);
	alphabet->domainsize = size-1;
	alphabet->mapsize = size-1;
	 
	for(i=1; i < size; i++) {
		switch(csv[i]->noofstrings) {
		case 1:
			alphabet->characters[i-1]=i;
			alphabet->mapdomain[i-1]=atoi(csv[i]->strings[0].str);
			break;
		case 2:
			trimmed = strtrimquote(space, csv[i]->strings[0].str, &csv[i]->strings[0].len);
			if (trimmed) {
				mapto=atoi(trimmed);
				FREEMEMORY(space, trimmed);
			} else {
				mapto = 0;
			}
			
			if(mapto == 0) mapto = i;
			alphabet->characters[i-1]=mapto;
			alphabet->mapdomain[i-1]=atoi(csv[i]->strings[1].str);
			break;

		default:
			fprintf(stderr, "malformed csv alphabet: %s. Exit forced.\n", filename);
			exit(EXIT_FAILURE);
		}
		destructStringset(space, csv[i]);
	}	

	/*remind this piece - the csv header*/
	destructStringset(space, csv[0]);
	FREEMEMORY(space, csv);
	return alphabet;
 }

 void sortMapdomain(void *space, FAlphabet *alphabet) {
	Uint i, *sorted, *mapdomain, *characters, 
		size = alphabet->mapsize;
	
	mapdomain = ALLOCMEMORY(space, NULL, Uint, size);
	characters= ALLOCMEMORY(space, NULL, Uint, size);
	
	sorted = quickSort(space, alphabet->mapdomain, size, cmp_int, NULL);
	for (i=0; i < size; i++) {
		
		characters[i] = alphabet->characters[sorted[i]];
		mapdomain[i] = alphabet->mapdomain[sorted[i]];
	}

	FREEMEMORY(space, alphabet->mapdomain);
	FREEMEMORY(space, alphabet->characters);
	FREEMEMORY(space, sorted);

	alphabet->characters=characters;
	alphabet->mapdomain=mapdomain;
 }


 void dumpAlphabet(FAlphabet *alphabet) {
	int i;
	
	printf("dumping alphabet: \n");
	for(i=0; i<alphabet->domainsize; i++) {
		printf("[%d]: %d\n", alphabet->characters[i], alphabet->mapdomain[i]);
	}
 }
 
