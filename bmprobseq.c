
/*
 *  bmprobseq.c
 *  perform a boyer-moore exact string
 *  matching on probability sequences
 *
 *  @author Steve Hoffmann
 *  @email shoffmann@zbh.uni-hamburg.de
 *  @date 12/07/06 22:55:11 CET
 *  
 */

 #include <stdlib.h>
 #include <stdio.h>
 #include <getopt.h>
 #include "lib/memory.h"
 #include "inc/prob_vec.h"
 #include "inc/prob_vec_i.h"
 #include "inc/falphabet.h"
 #include "inc/intsequence.h"
 #include "lib/mathematics.h"
 #include "lib/stringutils.h"
 #include "lib/fileio.h"
 #include "inc/createalphabet.h"
 #include "lib/boyermoore.h"
 #include "inc/encodeprobvec.h"
 
 const char *prog_bug="";
 const char *null_point="";

 static int verbose_flag;
 static struct option long_options[] =
 {
	/* These options set a flag. */
	{"verbose", no_argument,       &verbose_flag, 1},
	{"brief",   no_argument,       &verbose_flag, 0},	
	{"file",    required_argument, 0, 'f'},
	{0, 0, 0, 0}
 };


/*------------------------------- getProbChar --------------------------------
 *    
 * returns the position of character from a probability sequence over the 
 * given alphabet of size asize
 * 
 */
 
Uint
getProbChar (void *alphabet, Uint aszize, void *p, Uint pos)
{
  	/*at this point of the development alphabet is a stub*/
  	/*FAlphabet *a = (FAlphabet *)alphabet;*/
	IntSequence *seq = (IntSequence *)p;
	Uint ch = seq->sequence[pos];
	
	return ch;
}

/*---------------------------------- usage -----------------------------------
 *    
 * the usage function
 * 
 */
 
void
usage (char *name)
{
	printf("usage: %s -f <filename>\n", name);
  	return ;
}

/*----------------------------------- main -----------------------------------
 *    
 *  main function
 * 
 */
 
int
main (int argc, char** argv)
{	
 	Sint optindex, c;
	Uint i, matchcount=0, noofvecs=0;
#ifdef MEMMAN_H	
	Spacetable spacetab;
#endif	
	/*const Uint pattern[] = {156,144,132,142};*/
    const Uint pattern[] = {1,1,1,1,1,1,1,1,1,1,1,1,1};
	void *space = NULL;
	char *buf;
	char *pveclistfile = NULL;
	char *alphabetfile = NULL;
	vector_t *res = NULL;
	stringset_t **fn;
	IntSequence *tmplseq;
	IntSequence *patseq;
	FAlphabet *alphabet = NULL;
	
#ifdef MEMMAN_H 	
	initmemoryblocks(&spacetab, 100000);
	space = &spacetab;
#endif

    /*146  179  156  144  132  142  196*/
    patseq = initSequence(space);
	patseq->sequence = (Uint*) &pattern;
	patseq->length = 13;
	
    while(1) 
	{
		c=getopt_long(argc, argv, "f:a:", long_options, &optindex);
		if (c==-1) break;
		
		switch(c) {
			case 'f':
				pveclistfile = optarg;	
				break;
			case 'a':
				alphabetfile = optarg;	
				break;
			default:
				usage(argv[0]);
				exit (EXIT_FAILURE);
		}
	}
	if (pveclistfile==NULL) {
		usage(argv[0]);
		exit (EXIT_FAILURE);
	}
	

	fn=readcsv(space, pveclistfile, "", &noofvecs);
	if (alphabetfile != NULL) {
		alphabet = loadCSValphabet(space, alphabetfile);
		sortMapdomain(space, alphabet);
    }

	for(i=0; i<noofvecs; i++) 
	{		  
		tmplseq = loadSequence(space, SETSTR(fn[i],0));	
		res = boyerMoore(space, alphabet, 1408,
						tmplseq, tmplseq->length, patseq, patseq->length, getProbChar);
	
		if (LENGTHVEC(res)>0) {
	   		matchcount++;
	   		dumpVector(res);
	   		buf = printSequence(space, tmplseq, 80);
			printf("%s\n", buf);		
			FREEMEMORY(space, buf);
		}

		destructSequence(space, tmplseq);	
		destructVector(space, res);	
		destructStringset(space, fn[i]);
	}
	
	FREEMEMORY(space, patseq);
	FREEMEMORY(space, fn);
	
#ifdef MEMMAN_H
	activeblocks(space);
#endif
	
	printf("%d matches found in %d sequences.\n", matchcount, noofvecs);
	printf("Goodbye.\n");
	
	return EXIT_SUCCESS;
}

