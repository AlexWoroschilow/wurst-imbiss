#ifndef CANTOR_H
#define CANTOR_H

/*
 * cantor.h
 * declarations for coding cantor style
 *
 * @author Steve Hoffmann
 * @date Thu 23 Nov 2006
 *
 */

 #include "memory.h"
 #include "mathematics.h"
 #include "basic-types.h"

 Uint codeCantor(vector_t *);
 vector_t *decodeCantor(void *, Uint, Uint);

 #endif
