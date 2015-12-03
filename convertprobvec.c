
/*
 * convertprobvec.c
 * 
 * converting probability vectors to
 * probaility sequences
 *
 * @author Steve Hoffmann
 * @date 
 * 
 */

 #include <stdlib.h>
 #include <stdio.h>
 #include <string.h>
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
 #include "lib/vtprogressbar.h"
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
	{"alphabet", required_argument, 0, 'a'},
	{"outpath", required_argument, 0, 'o'},
	{0, 0, 0, 0}
 };


/*---------------------------------- usage -----------------------------------
 *    
 * the usage function
 * 
 */
 
void
usage (char *name)
{
    printf("usage: %s -f <pveclistfile> -a <alphabetfile>\n", name);
	return ;
}

/*----------------------------------- main -----------------------------------
 *    
 * the main function
 * 
 */
 
 int 
 main(int argc, char** argv) 
 {
	Uint noofvecs, i;
	Sint optindex, c;
	vector_t info;
#ifdef MEMMAN_H	
	Spacetable spacetab;
#endif
	void *space = NULL;
	char *url = NULL;
	char *outpath = NULL;
	char *pveclistfile = NULL;
	char *alphabetfile = NULL;
	char *vecext="vec";
	char *seqext="seq";	
	struct prob_vec *p_vec;
	IntSequence *sequence;
	FAlphabet *alphabet;	
	stringset_t *tok;
	stringset_t **fn;
	
#ifdef MEMMAN_H 
	initmemoryblocks(&spacetab, 1000);
	space = &spacetab;
#endif

    while(1) {
		c=getopt_long(argc, argv, "f:a:o:", long_options, &optindex);
		if (c==-1) break;
		
		switch(c) {
			case 'f':
				pveclistfile = optarg;	
				break;
			case 'a':
				alphabetfile = optarg;	
				break;
			case 'o':
				outpath = optarg;
				break;
			default:
				usage(argv[0]);
				exit (EXIT_FAILURE);
		}

	}

	if (pveclistfile==NULL || alphabetfile == NULL) {
		usage(argv[0]);
		exit (EXIT_FAILURE);
	}
	
	fn=readcsv(space, pveclistfile, ".", &noofvecs);
	alphabet = loadCSValphabet(space, alphabetfile);
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
		sequence->namelen= strlen(alphabetfile);
		
		COPYSTR(space, sequence->alphabetname, alphabetfile, 
							 strlen(alphabetfile));

		/*this is a potential security risk*/
		if (p_vec->compnd_len > 0) {
		  sequence->descrlen = p_vec->compnd_len-1; 
		  COPYSTR(space, sequence->description, p_vec->compnd, 
							 p_vec->compnd_len-1);
		} else {
		  sequence->descrlen = 14;			
		  COPYSTR(space, sequence->description, "descriptor n/a", 14);
		}

		sequence->urllen = SETSTRLEN(fn[i],0);
		COPYSTR(space, sequence->url, SETSTR(fn[i],0), 
						   SETSTRLEN(fn[i],0));
				
		SETSTR(fn[i],0) = concatdelim(space, SETSTR(fn[i],0), seqext,
							SETSTRLEN(fn[i],0), 3, '.');	
		SETSTRLEN(fn[i],0) += 4;
		
		if (outpath) {
			tok = tokensToStringset(space, "/", SETSTR(fn[i],0), 
					SETSTRLEN(fn[i],0));
			
			COPYSTR(space, url, outpath, strlen(outpath));
	
			url = concat(space, url, SETSTR(tok, tok->noofstrings-1), 
				strlen(url), SETSTRLEN(tok, tok->noofstrings-1));
		
			saveSequence(sequence, url);
			
			destructStringset(space, tok);
			FREEMEMORY(space, url);
			url = NULL;
		
		} else {
	
			saveSequence(sequence, SETSTR(fn[i],0));
		}
		
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


