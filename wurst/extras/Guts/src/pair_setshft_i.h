/* pair_setshft_i.h
 *  8 July 2004
 *  this is useful when mapping a threading of a sub-sequence 
 *  to the total alignment set for the parent sequence.
 */

#ifndef PAIR_SET_SHFT_I_H
#define PAIR_SET_SHFT_I_H

struct pair_set;

int pair_set_shift(struct pair_set *p_s, int shft);
  
#endif
