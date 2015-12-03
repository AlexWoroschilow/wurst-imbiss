/* pair_settrim_i.h
 *  8 July 2004
 * this will probably be most useful when you only want the
 * best region of an alignment.
 */

#ifndef PAIR_SET_TRIM_I_H
#define PAIR_SET_TRIM_I_H

struct pair_set;

struct pair_set *
pair_set_trim(struct pair_set *p_s, size_t start, size_t end);

#endif
