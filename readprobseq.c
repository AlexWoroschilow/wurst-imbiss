
/*
 *  readprobseq.c
 *  read a probility sequence
 *
 *  @author Steve Hoffmann
 *  @email shoffmann@zbh.uni-hamburg.de
 *  @date 12/05/06 22:55:11 CET
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
	char *filename;
	IntSequence *seq;
	
	c=getopt_long(argc, argv, "f:", long_options, &optindex);
	
	switch(c) {
		case 'f':
			filename = optarg;
			break;
		default:
			usage(argv[0]);
			exit (EXIT_FAILURE);
	}
  	
	seq = loadSequence(NULL, filename);
	printf("%s\n", printSequence(NULL, seq, 80));	
	return EXIT_SUCCESS;
}
