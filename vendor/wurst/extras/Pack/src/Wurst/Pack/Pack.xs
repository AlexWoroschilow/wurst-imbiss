#include "EXTERN.h"
#include "perl.h"
#include "XSUB.h"

#include "ppport.h"

#include "packw_obj_i.h"

#include "const-c.inc"

# /* this really ought to be an include taken from wurst
#  * (XStypemap.h)
#  */

typedef struct sec_s_data Sec_s_data;
typedef struct coord      Coord;
typedef struct seqprof    Seqprof;


MODULE = Wurst::Pack            PACKAGE = Wurst::Pack		

INCLUDE: const-xs.inc

# /* These were originally in Wurst.xs - handcoded xs. */
#/* ---------------- coord_pack        -------------------------
# * Pack a coord structure into a binary scalar
# */

void
coord_pack (topack)
        Coord *topack

        PREINIT:
        void *packed;

        PPCODE:
        packed = coord_pack(topack);
        if (packed!=NULL) {
                EXTEND(SP, 1);
                PUSHs(sv_2mortal(newSVpv((char *)packed, (*(size_t *)packed))));
                free(packed);
        } else {
                XSRETURN_UNDEF;
        }

#/* ---------------- coord_unpack        ------------------------
# * UnPack a coord structure from a binary scalar
# */

Coord *
coord_unpack (packed)
        char *packed;

#/* -------------------- pack_ss_data -------------------------- 
#*/
void
sec_s_pack (topack)
        Sec_s_data *topack

        PREINIT:
        void *packed;

        PPCODE:
        packed = sec_s_pack(topack);
        if (packed!=NULL) {
                EXTEND(SP, 1);
                PUSHs(sv_2mortal(newSVpv((char *)packed, (*(size_t *)packed))));
                free(packed);
        } else {
                XSRETURN_UNDEF;
        }

#/* -------------------- unpack_sec_s -------------------------- 
#*/

Sec_s_data *
sec_s_unpack (packed)
           char *packed;


#/* -------------------- pack/unpack checkpoints ---------------
# */
Seqprof *
seqprof_unpack( packed)
           char *packed;

void
seqprof_pack (topack)
        Seqprof *topack

        PREINIT:
        void *packed;

        PPCODE:
        packed = seqprof_pack(topack);
        if (packed!=NULL) {   
                EXTEND(SP, 1);
                PUSHs(sv_2mortal(newSVpv((char *)packed, (*(size_t *)packed))));
                free(packed);
        } else {
                XSRETURN_UNDEF;
        }

