 #ifndef MULTI_SEQ_H
 #define MULTI_SEQ_H

/*
 *	multiseq.h
 *  declarations for a datastructure containing
 *  multiple integer sequences
 * 
 *  @author Steve Hoffmann, shoffmann@zbh.uni-hamburg.de
 *  @company Center for Bioinformatics, Hamburg 
 *  @date 12/11/06 15:09:15 CET  
 *
 */

 #include "basic-types.h" 
 #include "intsequence.h"

 typedef struct {
	void *ref; 	
	BOOL refersFile;

 } SeqReference;


 typedef struct {
 	Uint numofsequences;
	Uint totallength;
	Uint *markpos;		/*markpos[i] is the position of a*/
						/*separator character between S_i and S_i+1*/
	Uint *sequences; 	/*array of concatenated sequences*/
    SeqReference *ref;  /*ref[i] points to the original sequence*/
						/*that starts at position markpos[i]*/
 } MultiIntSeq;


 void dumpMultiSeq (MultiIntSeq *);
 MultiIntSeq* concatIntSequences(void *, IntSequence **, Uint, Uint, Uint);
 void destructMultiSeq(void*, MultiIntSeq *);
 Uint getMultiSeqIndex(MultiIntSeq *, Uint *);
 Uint getMultiSeqRelPos(MultiIntSeq *, Uint *);

 #endif
