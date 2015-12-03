
/*
 *
 *  usage.h
 *  the usage message goes into here
 * 
 *  @author Steve Hoffmann, shoffmann@zbh.uni-hamburg.de
 *  @company Center for Bioinformatics, Hamburg 
 *  @date 04/13/07 19:40:26 CEST  
 *
 */

 #include <stdio.h>
 #include <stdlib.h>

 void usage(char *name) {

 fprintf(stderr,"usage:\n"); 
 fprintf(stderr," %s [-AFSdwhM:D:] [-b batchfile] [-s subfile] [-r reportfile]\n", name);
 fprintf(stderr,"                 -q queryfile  -s dbfile  -a alphabet  \n");
 fprintf(stderr,"                 -l matchlen -c minmatches \n"); 
 fprintf(stderr,"\nIndex mediated bayesian structure search (IMBISS) with Wurst \n");
 fprintf(stderr,"\n");
 fprintf(stderr," [INPUT]\n");
 fprintf(stderr," -q, --query        a file containing the query sequence\n");
 fprintf(stderr," -s, --sequences    a file containing the paths to all sequences of the db\n");
 fprintf(stderr," -a, --alphabet     a csv file containing the alphabet \n");
 fprintf(stderr," -x, --subfile      a csv file containing a substitution matrix (see descr.)\n");
 fprintf(stderr," -l, --matchlen     length of a match\n");
 fprintf(stderr," -p, --percent      use percentage of sequence length as seed\n");
 fprintf(stderr," -c, --minmatches   minimum number of matches\n");
 fprintf(stderr," -b, --batchfile    a csv file with queries\n");
 fprintf(stderr," [OPTIONS]\n");
 fprintf(stderr," -n, --maxresults   max. number of results returned\n");
 fprintf(stderr," -S, --allscores    use total scores for ranking only\n");
 fprintf(stderr," -B, --segmentscore use highest seed scores for ranking\n");
 fprintf(stderr," -F, --scorefilter  use total scores as filter\n");
 fprintf(stderr," -G, --segfilter    use high. seed scores as filter\n");
 fprintf(stderr," -A, --allalign     use local alignment scores only (default)\n");
 fprintf(stderr," -M, --match        alignment score for matches (default: 3)\n");
 fprintf(stderr," -D, --mismatch     alignment score for mismatches (default: -2)\n");
 fprintf(stderr," -d, --depictsw     sw alignments will be show for each match\n");
 fprintf(stderr," -w, --veggie       no wurst\n");
 fprintf(stderr," -r, --reportfile   a filename for csv reports\n");
 fprintf(stderr," -g, --gnuplot      draw matching statistics\n");
 fprintf(stderr," -h, --help         this message\n");
 fprintf(stderr," [ENVIRONMENT]\n");
 fprintf(stderr," set DISPLAY to redirect gnuplot output\n");
 fprintf(stderr," [REFERENCES]\n");
 fprintf(stderr,"  IMBISS is free software distributed under the terms of the GNU Public Licence\n");
 fprintf(stderr,"  (C) 2007 Steve Hoffmann (the Torda-Group, Centre for Bioinformatics, Hamburg) \n");
 fprintf(stderr," [BUGS]\n");
 fprintf(stderr,"  Please report bugs to <steve[dot]hoffmann[at]gmx[dot]de>\n");
 fprintf(stderr,"\n");

 exit(EXIT_FAILURE);
}



