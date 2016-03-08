#ifndef FILEIO_H
#define FILEIO_H

/*
 * fileio.h
 * declarations for file io
 *
 * @author Steve Hoffmann
 * @date Sat 25 Nov 2006
 *
 */

#ifndef ALLOCMEMORY
	#include "memory.h"
#endif

#include "stringutils.h"

char* readfile(void *, char *, Uint*);
stringset_t **readcsv(void *, char *, char*, Uint *);
void writeY(char *, double  *, Uint);
void writeXYUint(char *filename, Uint *X, Uint *Y, Uint len);


#endif
