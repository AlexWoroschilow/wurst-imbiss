/* pair_setshft.c
   Routines for renumbering and pruning pairsets
   Used by certain scripts for Wurst::Modeller
*/

#include <stdio.h>

#include "e_malloc.h"
#include "mprintf.h"
#include "pair_set.h"
#include "pair_setshft_i.h"

/* ---------------- pair_set_shift  --------------------------
 * Adds a signed constant to the sequence entry of a pair_set 
 * object.
 */

int pair_set_shift ( struct pair_set *p_s, int shft) {
    struct apair *pb, *pe;
    if (p_s != NULL) {
        if (shft!=0) {
            pb = p_s->pairs;
            pe = pb+p_s->n;
            while (pb<pe) {
                if ((pb->a!=GAP_INDEX) && ((pb->a+=shft) < 0))
                    err_printf("pair_set_shift","Warning - pair_set entry has gone negative! Probably a bug.\n");
                pb++;
            }
        }
    } else {
        return 0;
    }
    return 1;
}
