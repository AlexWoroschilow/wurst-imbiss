/* 13 June 2004
   pair_set_xchange
   this is not a normal threading operation!
*/

   
#include <stdio.h>

#include "e_malloc.h"
#include "mprintf.h"
#include "pair_set.h"
#include "pair_setxch_i.h"
/* ---------------- pair_set_xchange  --------------------------
    interchanges the position of the pairs of residue numbers in
    an alignment, returning a new pairset.
 */
struct pair_set *pair_set_xchange (struct pair_set *p_s) {
  struct pair_set *nw;
  nw = NULL;
  if (p_s!=NULL) {
      if ((nw = E_MALLOC(sizeof(*nw)))==NULL)
	  err_printf("pair_set_xchange","No memory for pair_set structure.\n");
      nw->pairs = NULL;
      nw->n = 0;
      nw->score = nw->smpl_score = 0.0;
      if (p_s->n>0) {
	  struct apair *p,*n,*plast;
	  if ((nw->pairs = E_MALLOC(sizeof(*nw->pairs)*p_s->n))==NULL)
	      err_printf("pair_set_xchange","No memory for pair_set pairs.\n");
	  plast = ((p = p_s->pairs)+p_s->n);
	  nw->n = p_s->n;
	  n = nw->pairs;
	  while (p<plast) {
	      n->a = p->b;
	      n->b = p->a;
	      p++; n++;
	  };
	  nw->score = p_s->smpl_score;
	  nw->smpl_score = p_s->score;
      }
      
  }
  return nw;
}

