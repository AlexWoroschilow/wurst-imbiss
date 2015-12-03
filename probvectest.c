
/*
 * alphabet_test.c
 * testing alphabet routiness
 *
 */

 #include <stdlib.h>
 #include <stdio.h>
 #include <string.h>
 #include "lib/memory.h"
 #include "inc/prob_vec.h"
 #include "inc/prob_vec_i.h"
 #include "inc/falphabet.h"
 #include "inc/intsequence.h"
 #include "lib/mathematics.h"
 #include "lib/stringutils.h"
 #include "lib/fileio.h"
 #include "inc/createalphabet.h"
 #include "lib/vtprogressbar.h"
 #include "inc/encodeprobvec.h"
 
 const char *prog_bug="";
 const char *null_point="";

/*----------------------------------- main -----------------------------------
 *    
 * the main function
 * 
 */
 
 int 
 main(int argc, char** argv) 
 {
	Uint noofvecs, i;
	vector_t info;
#ifdef MEMMAN_H
	Spacetable spacetab;
#endif 
	void *space = NULL;
	char *vecext="vec";
	char *seqext="seq";	
	struct prob_vec *p_vec;
	IntSequence *sequence;
	FAlphabet *alphabet;	
	stringset_t** fn;
	
#ifdef MEMMAN_H
	initmemoryblocks(&spacetab, 1000);
	space = &spacetab;
#endif
	
	fn=readcsv(space, "listshort.pvecs", ".", &noofvecs);
	alphabet = loadCSValphabet(space, "alphabet_t100.csv");
	sortMapdomain(space, alphabet);
   
	for(i=0; i<noofvecs; i++) 
	{	
	  
		INITVECTOR(&info);	
	  	SETSTR(fn[i],0) = concatdelim(space, SETSTR(fn[i],0), vecext,
							  SETSTRLEN(fn[i],0), 3, '.');
	   
		p_vec = prob_vec_read (SETSTR(fn[i],0));	
			
		if (p_vec->mship == NULL) 
			prob_vec_expand(p_vec);
			
		sequence = encode_prob_vec(space, alphabet, p_vec, 0, 0, 
									cantorchar, &info);
		sequence->info = (Uint*) info.elements;
		sequence->descrlen = SETSTRLEN(fn[i],0);
		COPYSTR(space, sequence->description, SETSTR(fn[i],0), 
							 SETSTRLEN(fn[i],0));
		SETSTR(fn[i],0) = concatdelim(space, SETSTR(fn[i],0), seqext,
							 SETSTRLEN(fn[i],0), 3, '.');	


		saveSequence(sequence, SETSTR(fn[i],0));
		
		destructSequence (space, sequence);		
		prob_vec_destroy (p_vec);			 
		destructStringset (space, fn[i]);
		progressBarVT("probability vectors converted", noofvecs-1, i, 25);	
	}
	
	printf("\nexit.\n");
	FREEMEMORY(space, fn);
	destructAlphabet(space, alphabet);
		
	return EXIT_SUCCESS;
 }


