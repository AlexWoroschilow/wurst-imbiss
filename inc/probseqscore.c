
/*
 *  probseqscore.c
 *  probability structure sequence
 *  scoring
 *
 *  @author Steve Hoffmann
 *  @email shoffmann@zbh.uni-hamburg.de
 *  @date 01/26/07 23:02:33 CET
 *  
 */
 #include "stdio.h" 
 #include "intsequence.h"
 #include "sufmatch.h"

 double
 directscore(void *space, Matchtype *m, IntSequence **seq, int len) {
 	Uint i, sum = 0;
	double score;
	
	for (i=0; i < m->count; i++) {

		sum += seq[m->id]->sequence[m->pos[i]];
	}
	
	score = (double) sum / len;
	return score;
 }


