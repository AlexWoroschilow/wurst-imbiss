 #ifndef INTSEQUENCE_H
 #define INTSEQUENCE_H

/*
 * probsequence.h
 * declaration of a probability sequence
 * and functions working on it
 *
 * @author Steve Hoffmann
 * @date Mon 27 Nov 2006
 *
 */

 #include "basic-types.h"


 typedef struct {
	Uint descrlen;			
    Uint namelen;
	Uint urllen;
	
	char *description; 		/*a description*/
    char *alphabetname;		/*the name of the corresponding alphabet*/
	char *url;				/*the name of the sequences url*/
	
	Uint *sequence,			/*the sequence itself*/
		 *info;				/*additional information*/
	Uint length;
	
 } IntSequence;
 
 void destructSequence(void *, IntSequence *);
 IntSequence* initSequence(void *);
 char* printSequence(void *, IntSequence *, Uint);
 void dumpSequence(IntSequence *s);
 void saveSequence (IntSequence *s, char *filename); 
 IntSequence* loadSequence (void *space, char *filename);
 char * printAlignment (void *, int *, Uint, IntSequence *, IntSequence *, 
	 Uint);
 IntSequence** createSequenceHash(void *, Uint);

 
 #endif
