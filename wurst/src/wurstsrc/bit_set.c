/*
 * bit_set.c
 *
 *  Created on: Oct 13, 2011
 *      Author: ibondarenko
 */

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include "e_malloc.h"
#include "bit_set.h"

size_t
bit_set_getsize(size_t size) {

    return (((size % CHAR_BIT) != 0) ? size / CHAR_BIT + 2 : size / CHAR_BIT
            + 1);

}

char*
bit_set_calloc(size_t size) {

    return ((char*) E_MALLOC (bit_set_getsize(size) ));

}

void
bit_set_free(char* bs) {

    free_if_not_null(bs);
    bs = 0;
}

void bit_set_bitset(char* bs, size_t pos) {

    char XXXmask = 1;
    XXXmask <<= ((pos) % CHAR_BIT);
    bs[(pos) / CHAR_BIT] |= XXXmask;

}

unsigned char bit_set_ifbitset(char* bs, size_t pos) {

    char XXXmask = 1;
    XXXmask <<= ((pos) % CHAR_BIT);
    return ((bs[(pos) / CHAR_BIT] & XXXmask) == XXXmask);

}

void bit_set_clean(char* bs, size_t pos) {

    char XXXmask = 1;
    XXXmask <<= ((pos) % CHAR_BIT);
    bs[(pos) / CHAR_BIT] ^= XXXmask;

}

