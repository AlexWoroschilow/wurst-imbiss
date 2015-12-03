/*
 * bit_set.h
 *
 *  Created on: Sep 8, 2011
 *      Author: ibondarenko
 */

#ifndef BIT_SET_H_
#define BIT_SET_H_

size_t
bit_set_getsize(size_t size);

char*
bit_set_calloc(size_t size);

void
bit_set_free(char* bs);

void
bit_set_bitset(char* bs, size_t pos);

unsigned char
bit_set_ifbitset(char* bs, size_t pos);

void
bit_set_clean(char* bs, size_t pos);

#endif /* BIT_SET_H_ */
