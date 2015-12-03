/* pair_settrim.c
   Routines to extract a sub-alignment for a range of
   sequence or structure.
*/

#include <stdio.h>

#include "e_malloc.h"
#include "mprintf.h"
#include "pair_set.h"
#include "pair_settrim_i.h"

/* ---------------- pair_set_trim  --------------------------
 * extract a range of pairs aligning to a specified range of
 * sequence, and return them as a new pairset.
 * end=0 and we define the range from start->end of sequence
 * start = 0 and we define the range from start->end
 */
/* perhaps need two more options for specifying structure range */
struct pair_set *
pair_set_trim(struct pair_set *p_s, size_t start, size_t end) {
    struct pair_set *n;
    size_t strt, ends;
    n = NULL;
    
    
    if (p_s != NULL) {
        /* find all pairs in range */
        /* make array and copy */
        if (!start)
            err_printf("pair_set_trim", 
                       "Broken script - start is in range 1..seq_length\n");
        
        if (end && (end<=start)) 
            err_printf("pair_set_trim", 
                       "Broken script - end is less than start\n");
        
        strt = start-1;
        ends = end-1;
        if ((n = E_MALLOC(sizeof(*n)))==NULL)
            err_printf("pair_set_trim","No memory for pair_set structure.\n");
        
        n->pairs = NULL;
        n->n = 0;
        n->score = n->smpl_score = 0.0;
        if (p_s->n>0) {
            struct apair *p,*en,*p1, *p2;
            p=p_s->pairs;
            en=p+p_s->n;
            while ((p<en) && (p->a<(int)strt))
                p++;
            if (p<en) {
                p1 = p;
                while ((p<en) && (p->a<=(int)ends)) {
                    n->n++;
                    p++;
                }
                p2 = p;
            } else {
                err_printf("pair_set_trim", 
                           "start is longer than aligned sequence!\n");
            };
            
            if ((n->pairs = E_MALLOC(sizeof(*n->pairs)*n->n))==NULL)
                err_printf("pair_set_xchange",
                           "No memory for pair_set pairs.\n");
            p = p1;
            en = n->pairs;
            while (p<p2) {
                *en = *p;
                p++;
                en++;
            }
            /* don't modify score at all at the moment */
        }

    }
    return n;
}
