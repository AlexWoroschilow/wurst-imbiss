#ifndef CREATEALPHABET_H
#define CREATEALPHABET_H

#include "basic-types.h"
#include "falphabet.h"

/*
 * createalphabet.h
 * reading and creating of flexible alphabets
 *
 * @author Steve Hoffmann
 * @date Sat 25 Nov 2006
 *
 */

 FAlphabet *loadCSValphabet(void *, char*); 
 void sortMapdomain(void *, FAlphabet *);
 void dumpAlphabet(FAlphabet*); 


 #endif
