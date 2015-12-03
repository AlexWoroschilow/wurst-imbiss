/* wurstguts_i.h
 * 
 * 18th June 2004
 */

#ifndef WURSTGUTS_I_H
#define WURSTGUTS_I_H
struct scor_set;
struct pair_set;
struct score_mat;
struct coord;

/* functions to delete sequence or structure stretches */

#include "coord_deletion_i.h"

/* functions for analysing local scores for an alignment */

#include "scoranlys_i.h"

/* mci-output for interaction matrices
   #include "mcigraph_io_i.h"
   this is the interface for general float** matrices
*/

/* pairset functions */

#include "pair_setxch_i.h"
#include "pair_setshft_i.h"
#include "pair_settrim_i.h"
/* seqprof functions */
#include "seqprof_trim_i.h"


/* comparison of two alignments for consistency */

#include "model_pscmp_i.h"


#endif
