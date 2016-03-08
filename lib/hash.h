#ifndef HASH_H
#define HASH_H

/*
 *
 *	hash.h
 *  declarations for various hash functions 
 * 
 *  @author Steve Hoffmann, shoffmann@zbh.uni-hamburg.de
 *  @company Center for Bioinformatics, Hamburg 
 *  @date 04/03/07 19:21:09 CEST  
 *
 */

unsigned long 
knuthhash (unsigned char *, int);
	
unsigned long
jenkinshash (unsigned char *, int); 

#endif
