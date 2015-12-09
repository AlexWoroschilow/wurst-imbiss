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

	char *description; /*a description*/
	char *alphabetname; /*the name of the corresponding alphabet*/
	char *url; /*the name of the sequences url*/

	Uint *sequence, /*the sequence itself*/
	*info; /*additional information*/
	Uint length;

} IntSequence;

void destructSequence(void *, IntSequence *);
IntSequence* initSequence(void *);
char* printSequence(void *, IntSequence *, Uint);
void dumpSequence(IntSequence *s);
void saveSequence(IntSequence *s, char *filename);
IntSequence* loadSequence(void *space, char *filename);
char * printAlignment(void *, int *, Uint, IntSequence *, IntSequence *, Uint);
IntSequence** createSequenceHash(void *, Uint);

/**
 * Convert sequence to printable format
 */
char * sequence_userfriendly(void *space, IntSequence *s, Uint cols);

/**
 * Load sequence from given file using wurst
 * the normal amino-acids based sequence should be loaded
 */
IntSequence* sequence_load_pdb(void *space, char *filename);

/**
 * Load sequence from given file using wurst
 * the structure-alphabet based sequence should be loaded
 */
IntSequence* sequence_load_wurst(void *space, char *filename);

/**
 * Load a sequence list from csv file
 * using some load handler like 'sequence_load_wurst'
 * or 'sequence_load_pdb'
 */
IntSequence ** sequence_load_csv(void *space, char* filename, char *delimeter, Uint *linecount,
		IntSequence* (*loader)(void *space, char *filename));

#endif
