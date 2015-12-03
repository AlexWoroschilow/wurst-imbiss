#ifndef BOYER_MOORE_H
#define BOYER_MOORE_H

/*
 * boyermoore.h
 * declarations for a boyer moore over
 * flexible alphabets 
 *
 * @author Steve Hoffmann
 * @date Sun 23 Nov 2007
 *
 */
 
 #include "basic-types.h"
 #include "mathematics.h"

 Uint getASCIIcode(void *alphabet, Uint, void *, Uint);
 vector_t *boyerMoore(void *, void *, Uint, void *, Uint, void *, 
 					  Uint, Uint (*charpos)(void *, Uint, void *, Uint));


#endif

