 
#ifndef DPALIGN_H
#define DPALIGN_H

/*
 *
 *	dpalign.h
 *  declarations for routines for dynamic programming  
 * 
 *  @author Steve Hoffmann, shoffmann@zbh.uni-hamburg.de
 *  @company Center for Bioinformatics, Hamburg 
 *  @date 03/19/07 00:19:57 CET  
 *
 */

#include "basic-types.h"

typedef Uint symtype;

Uint edist(void *, symtype *, Uint, symtype *, Uint, Uint, Uint *, Uint);
int constscr (symtype, symtype, void *);
int* swmatrix (void *, symtype*, Uint, symtype*, Uint, int,
			   Sint (*sub)(symtype, symtype, void *), void *);

int* swgapless (void *, symtype*, Uint, symtype*, Uint, 
			   Sint (*sub)(symtype, symtype, void *), void *);

int* swgaplesstraceback (void *, int *,  
			 symtype *, Uint, symtype *, Uint, 
			 Sint (*sub)(symtype, symtype, void *), void *, int*);

#endif
