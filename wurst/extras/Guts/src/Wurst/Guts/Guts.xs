#include "EXTERN.h"
#include "perl.h"
#include "XSUB.h"

#include "ppport.h"

typedef struct coord      Coord;
typedef struct seq        Seq;
typedef struct seq_array  Seq_array;
typedef struct sub_mat    Sub_mat;
typedef struct score_mat  Score_mat;
typedef struct pair_set   Pair_set;
typedef void              Param;
typedef struct FXParam    FXParam;
typedef float             RSParam;
typedef struct seqprof    Seqprof;
typedef struct scor_set   Scor_set;
#include "wurstguts_i.h"

#include "const-c.inc"

MODULE = Wurst::Guts		PACKAGE = Wurst::Guts		



INCLUDE: const-xs.inc

#/* ---------------- coord_segment   --------------------------
# */


Coord *
coord_segment (oc, start, end)
        Coord *oc
        size_t start
        size_t end;


#/* ---------------- coord_deletion   --------------------------
# */


Coord *
coord_deletion (oc, start, slength,qlength)
        Coord *oc
        size_t start
        size_t slength
        size_t qlength;

#/* -------------------- seq_deletion --------------------
# * Undocumented call used by coord_deletion
# */

Seq *
seq_deletion (oseq, start, slength, qlength)
    Seq *oseq
   size_t start
   size_t slength
   size_t qlength;
    

#/* ---------------- ps_cmp_and_model --------------------------
# */


void
ps_cmp_and_model ( a, score_a, Sa, b, score_b, Sb)
    Pair_set *a
    Score_mat *score_a
    Coord *Sa
    Pair_set *b
    Score_mat *score_b
    Coord *Sb

# /* Found this magic from the responses to: 
# * http://archive.develooper.com/perl-xs@perl.org/msg00514.html
# * Specifically :
# * http://archive.develooper.com/perl-xs@perl.org/msg00516.html
# * { No thanks to the other reply from Hisrch where he misses
# *   the point entirely! }
# * For us - we are returning blessed object pointers- hence the
# * _pv on 'sv_setref_**' and the "CoordPtr" entry. */

    PREINIT:
        SV *model1;
        SV *model2;
        SV *model3;
        SV *model4;
        SV *conserved_pairs;
        char *f_string;
        struct pair_set *cons_pset;
        struct coord **models;
    PPCODE:
        models = ps_cmp_and_model (a, score_a, Sa, b, score_b, Sb, &f_string, &cons_pset);

        if (models!=NULL) {
            model1 = sv_newmortal();
            model2 = sv_newmortal();
            model3 = sv_newmortal();
            model4 = sv_newmortal();
            sv_setref_pv(model1, "CoordPtr",(void *) models[0]);
            sv_setref_pv(model2, "CoordPtr",(void *) models[1]);
            sv_setref_pv(model3, "CoordPtr",(void *) models[2]);
            sv_setref_pv(model4, "CoordPtr",(void *) models[3]);
            conserved_pairs = sv_newmortal();    
            sv_setref_pv(conserved_pairs,"Pair_setPtr",(void *) cons_pset);
            EXTEND(SP, 6);
            PUSHs(model1);
            PUSHs(model2);
            PUSHs(model3);
            PUSHs(model4);
            PUSHs(sv_2mortal(newSVpv (f_string,0)));
            PUSHs(conserved_pairs);
            free(f_string);
            free(models);
        } else {
            XSRETURN_UNDEF;
        }


#/* ---------------- pair_set_xchange --------------------------
# */

Pair_set *
pair_set_xchange (toinvert)
    Pair_set *toinvert;


#/* ---------------- scor_set_rs --------------------
# */


Scor_set *
scor_set_rs( coord, params )
        Coord *coord
        RSParam *params;

#/* ---------------- mci_contact_rs --------------------
# */

        
int
mci_contact_rs( fname, coord, params )
        const char *fname
        Coord *coord
        RSParam *params;

#/* ---------------- mci_sc_n_rs --------------------
# */                 

int
mci_sc_n_rs( fname, lscores, coord, params )
        const char *fname
        Scor_set *lscores
        Coord *coord
        RSParam *params;

#/* ----pair_set_shift ---------------------------------
# */

int pair_set_shift ( pairset, offset)
        Pair_set *pairset
        int offset;
        
#/* ----pair_set_trim ----------------------------------
# */
Pair_set *
pair_set_trim( pairset, start, end)
        Pair_set *pairset
        size_t start
        size_t end;


#/* ----seqprof_trim------------------------------------
# */

Seqprof *
seqprof_trim ( profile, start, end)
        Seqprof *profile
        size_t  start
        size_t  end;

#/* ----------------------------------------------------
Seqprof *
seqprof_merge ( profile1, profile2);
        Seqprof *profile1
        Seqprof *profile2;

# */
#/* ----------------------------------------------------
Coord *
coord_merge ( coord1, coord2);
        Coord *coord1
        Coord *coord2;

# */
#/* ----------------------------------------------------
# */
#/* ----------------------------------------------------
# */
#/* ----------------------------------------------------
# */
