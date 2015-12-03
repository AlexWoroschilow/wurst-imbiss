
/*
 *  createsubmatrix.c
 *  create a substitution matrix
 *
 *  @author Steve Hoffmann
 *  @email shoffmann@zbh.uni-hamburg.de
 *  @date 02/19/07 17:49:02 CET
 *  
 */

 #include <stdlib.h>
 #include <stdio.h>
 #include <string.h>
 #include <getopt.h>
 #include <math.h>
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
 #include "inc/imsubst.h"
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
	Uint noofvecs, i, j;
	Sint optindex, c;
	vector_t info;
#ifdef MEMMAN_H	
	Spacetable spacetab;
#endif
	void *space = NULL;
	char *pveclistfile = NULL;
	char *alphabetfile = NULL;
	char *vecext="vec";	
	unsigned char scores=0;
	struct prob_vec *p_vec;	
	FAlphabet *alphabet;	
	double **avg;
	stringset_t** fn;
	Uint sum=0;
	
#ifdef MEMMAN_H 
	initmemoryblocks(&spacetab, 1000);
	space = &spacetab;
#endif

    while(1) {
		c=getopt_long(argc, argv, "f:a:", long_options, &optindex);
		if (c==-1) break;
		
		switch(c) {
		    case 's':
			  	scores = 1;
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

	if (pveclistfile==NULL || alphabetfile == NULL) {
		usage(argv[0]);
		exit (EXIT_FAILURE);
	}
	
	fn=readcsv(space, pveclistfile, ".", &noofvecs);
	alphabet = loadCSValphabet(space, alphabetfile);
	sortMapdomain(space, alphabet);
   
	avg = ALLOCMEMORY(space, NULL, double*, 1408);
	for (i=0; i < 1408; i++) {
		avg[i] = ALLOCMEMORY(space, NULL, double, 309);
		for(j=0; j < 309; j++) {
			avg[i][j] = 0;
		}
	}
	
	for (i=0; i<noofvecs; i++) 
	{	
	  
		INITVECTOR(&info);	
	  	SETSTR(fn[i],0) = concatdelim(space, SETSTR(fn[i],0), vecext,
			SETSTRLEN(fn[i],0), 3, '.');
	   
		p_vec = prob_vec_read (SETSTR(fn[i],0));	
	
		if (p_vec->mship == NULL) 
			prob_vec_expand(p_vec);
			
		avg = avgpvec (space, alphabet, avg, p_vec, 0, 0, 
			cantorchar, &info);
			
		prob_vec_destroy (p_vec);			 
		destructStringset (space, fn[i]);
		progressBarVT("probability vectors scanned", noofvecs-1, i, 25);	
	}

	for(i=0; i < 1408; i++) {
	//printf("%f;",  log10(m[i][j]*308));}	
		//printf("\n");
	}


	for(i=0; i < 1408; i++) {
		sum+=avg[i][308];
	}
	
	
	for(i=0; i < 1408; i++) {
		if (scores)
	  	   	printf("%f\n", 
				(log10(avg[i][308]/sum)*(-1)));
		else
			printf("%f\n",
				(avg[i][308]/sum));
	}
	

	printf("\nexit.\n");
	FREEMEMORY(space, fn);
	destructAlphabet(space, alphabet);
		
	return EXIT_SUCCESS;
 }


