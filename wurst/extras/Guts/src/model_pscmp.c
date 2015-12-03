/* ---------------- ps_cmp_and_model  --------------------------
 * give me a model of consistent aligned residues from sequence
 to partisan's structure (!seq_targ) or the native structure of all
 aligned residues

 constructor and
 helper functions have been duplicated (ie score_smpl) for
 speed.
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "coord.h"
#include "coord_i.h"
#include "e_malloc.h"
#include "model.h"
#include "mprintf.h"
#include "pair_set.h"
#include "read_seq_i.h"
#include "seq.h"
#include "score_mat.h"
#include "model_pscmp_i.h"

/*
 * Return an empty model structure
 */
static struct coord *
model_dummy(struct coord *template) {
    const char *this_sub = "model_dummy";
    struct coord *mdl;
    if ((mdl = coord_template (template,0))==NULL) {
        /* Allocates the actual coord struct */
        err_printf(this_sub, "No Memory for model\n");
        free(mdl);
        return NULL;
    }
    /* Set size to zero */
    mdl->seq = NULL;
    mdl->size = 0;

    return mdl;
}

static struct coord *
model_fromAnnotAlign(char *a_code, int *A_pos,
                     int szSeq, int m,
                     struct coord *seqStr, struct coord *targStr,
                     int seq_targ) {

    struct coord *mdl, *cfme;
    int i,j;
    const char *this_sub = "model_fromAnnotAlign";

    cfme = (seq_targ) ? targStr : seqStr;
    if ((mdl = coord_template (cfme,m))==NULL) {
        /* Allocates the actual coord struct */
        err_printf(this_sub, "No Memory for model\n");
        free(mdl);
        return NULL;
    }

    mdl->seq = seq_copy (seqStr->seq);
    for ( i=0, j=0; i<szSeq; i++) {
        if ((a_code[i]=='C') || ((!seq_targ) && (A_pos[i]>=0))) {
            size_t ndx_a = i;  /* sequence index */
            size_t ndx = ( seq_targ ) ? A_pos[i] : i;    /* structure index */
            mdl->rp_ca[j] = cfme->rp_ca[ndx];
            mdl->rp_cb[j] = cfme->rp_cb[ndx];
            mdl->rp_n[j]  = cfme->rp_n[ndx];
            mdl->rp_c[j]  = cfme->rp_c[ndx];
            mdl->rp_o[j]  = cfme->rp_o[ndx];

            mdl->icode[j] = ' ';
            if (cfme->psi)
                mdl->psi[j] = cfme->psi[ndx];
            if (cfme->sec_typ)
                mdl->sec_typ[j] = cfme->sec_typ[ndx];

            mdl->orig[j]  = i+1;  /* Model gets number from sequence */
            mdl->seq->seq[j] = seqStr->seq->seq[ndx_a]; /* Edit model's coord sequence */
            j++;
        }
    }
    mdl->seq->seq[j]='\0';
    mdl->seq->length=j;
    mdl->seq->seq = E_REALLOC(mdl->seq->seq, (j+1) * sizeof (mdl->seq->seq[0]));

    return mdl;
}


/* from align.c
   score_smpl() walks through a score matrix, based on the numbers in a pair_set
   and fills in the smpl_score element of the structure
*/
static void
smpl_score (struct pair_set *pair_set, const struct score_mat *smat)
{
    /* We account for the extra row and column (see align.c c.a. line 435)
       in the score matrix, here.
    */
    struct apair *p, *plast;
    float **mat = smat->mat;
    float smpl = 0.0;
    p = pair_set->pairs;
    plast = p + pair_set->n;
    
    for ( ; p < plast; p++) {
        int nc, nr;
        nc = p->a;
        nr = p->b;
        if (nc == GAP_INDEX)
            continue;
        if (nr == GAP_INDEX)
            continue;
        smpl += mat[nc+1][nr+1];
    }
    pair_set->smpl_score = smpl;
}

/*
 *   test consistency of sequence alignment (onto potential homolog) against original structure
 *   Make model of (consistent pairs of q_a on Sa and Sb, consistent pairs of q_b on Sb and Sa)
 *   currently only q_a on Sb and q_b on Sa are made.

*/

/* make a fragment string */
static void add_fragment(char **frag_set, int f_Sa, int f_Sb, int f_len) {
    char frag_buffer[25]; /* scratch for forming each fragment */
    if (f_len==1)
        sprintf(frag_buffer, "(%i,%i)",f_Sa, f_Sb);
    else
        sprintf(frag_buffer, "(%i,%i,%i)",f_Sa, f_Sb,f_len);
    *frag_set = E_REALLOC(*frag_set, strlen(frag_buffer)+strlen(*frag_set)+1);
    strcat(*frag_set, frag_buffer);
}

struct coord **
ps_cmp_and_model (struct pair_set *a, struct score_mat *Score_a, struct coord *Sa,
                  struct pair_set *b, struct score_mat *Score_b, struct coord *Sb, 
                  char **frag_string, struct pair_set **cons_pairs)
{
    /*
      a is always sequence of a to structure Sb
      b is always sequence of b to structure Sa
      We compute a consistency measure for the moment.
    */
    const char *this_sub = "pair_set_compare";
    

    char *frag_set;
    int *A_algn, *B_algn;
    char *A_code, *B_code;
    size_t i,l,m,Al_a,Al_b;
    size_t f_Sa,f_Sb,f_len;
    size_t szA,szB;
    struct coord **model;
    struct pair_set *cnsvd;
    struct apair *cnsv,*cons_set,*cs_p,*ccnsv,*bcnsv;
    struct Sippl_Res {
        int *s_distr[2];
        int low_s,high_s;
        int n_incorrect[2];
        int n_extra[2];
        int n_missed[2];
        int n_correct;
        int ng_correct;
        int Lc,Lp;
        double S_bar[2];
    } r;


    /* Codify Alignments to sequence/coordinate numbers
     */

    szA = Sa->seq->length;
    szB = Sb->seq->length;
    frag_set = E_MALLOC((sizeof(char)*3));
    A_algn=(int *)E_MALLOC(sizeof(*A_algn)*4*szA);
    B_algn=(int *)E_MALLOC(sizeof(*B_algn)*4*szB);
    A_code=(char *)E_MALLOC(sizeof(*A_code)*(1+szA));
    B_code=(char *)E_MALLOC(sizeof(*A_code)*(1+szB));

    frag_set[0]='\0';
    for (i=0;i<(4*szA); i++) {
        A_algn[i]=-2;
        if (i<szA)
            A_code[i]=' ';
        if (i>(2*szA))
            A_algn[i]=-1; /* Gap Mapping for pair_set pairs */
    }
    for (i=0;i<(4*szB); i++) {
        B_algn[i]=-2;
        if (i<szB)
            B_code[i]=' ';
        if (i>(2*szB))
            B_algn[i]=-1; /* Gap Mapping for pair_set pairs */
    }

    B_code[szB]=A_code[szA] = '\0';
    Al_a=Al_b=0;
    r.ng_correct = r.n_correct=0;
    for (i=0;i<a->n;i++) 
        if (a->pairs[i].a!=GAP_INDEX) {
            l=a->pairs[i].a;
            if (a->pairs[i].b!=GAP_INDEX) {
                Al_a++;
                m = a->pairs[i].b;
                A_algn[l] = m;
                B_algn[m + szB] = l;
            } else
                A_algn[l] = GAP_INDEX;
        } else if (a->pairs[i].b!=GAP_INDEX) {
            B_algn[a->pairs[i].b+szB] = GAP_INDEX;
        }
    
    l=m=0;
    for (i=0;i<b->n;i++)
        if (b->pairs[i].a!=GAP_INDEX) {
            m = b->pairs[i].a;
            if (b->pairs[i].b!=GAP_INDEX) {
                l = b->pairs[i].b;
                Al_b++;
                B_algn[m] = l;
                A_algn[l+szA] = m;
            } else
                B_algn[m] = GAP_INDEX;
        } else if (b->pairs[i].b!=GAP_INDEX) {
            A_algn[b->pairs[i].b+szA] = GAP_INDEX;
        }

    /* Score Consistency.
     */
    f_Sa = f_Sb = f_len = 0;
    cons_set = E_MALLOC(sizeof(*cnsv)*(szA+szB));
    cnsv = E_MALLOC(sizeof(*cnsv)*szA*2);

    cs_p = cons_set;
    ccnsv = cnsv;
    bcnsv = cnsv+szA;
    
    for (l=0;l<szA; l++) {
        if (A_algn[l]==A_algn[l+szA]) {
            if ((A_algn[l]==-2)||(A_algn[l]==GAP_INDEX)) {
                r.ng_correct++;
                A_code[l] = (A_algn[l]==GAP_INDEX) ? 'g':'c';
                if (f_len) {
                    add_fragment(&frag_set, f_Sa, f_Sb, f_len);
                    f_len=0;
                }
            }
            else {
                if (f_len == 0) {
                    f_len = 1;
                    f_Sa = l;
                    f_Sb = A_algn[l];
                } else {
                    if ((l-f_Sa) != (A_algn[l]-f_Sb)) {
                        add_fragment(&frag_set, f_Sa, f_Sb, f_len);
                        f_len=0;
                        f_Sa= l;
                        f_Sb = A_algn[l];
                    }
                    f_len++;
                }
                r.n_correct++;
                A_code[l] = 'C';
                cs_p->a = l;
                cs_p->b = A_algn[l];
                cs_p++;
                ccnsv->a = l;
                ccnsv->b = m = A_algn[l]; /* Assume correct pairs contained in matrix? */
                bcnsv->b = l;
                bcnsv->a = m;
                ccnsv++;
                bcnsv++;
            }
        } else {
            if (f_len) {
                /* Finish fragment */
                add_fragment(&frag_set, f_Sa, f_Sb, f_len);
                f_len=0;
            }
            
            /* Dont Care ? for the mo.
             *  put a truth table in for the code applied to these
             *  conditions eventually
             
             if ((A_algn[l]<0) || (A_algn[l+szA]<0)) {
             char ca,cb;
             switch (A_algn[l]) {
             case GAP_INDEX:
             ca = 'A';
             break;
             case -2:
             ca = 'U';
             break;
             default:
             ca = '.';
             }
             switch (A_algn[l+szA]) {
             case GAP_INDEX:
             cb = (ca=='U') ? 'n' : 'A';
             break;
             case -2:
             cb = (ca=='A') ? 'N' : 'U';
             break;
             default:
             cb = (ca=='.') ? ''.';
             }
             A_code[l] = (ca=='.') ? ((cb=='A') : ' ';
             ((A_algn[l]==GAP_INDEX)
             && (A_algn[l+szA]==GAP_INDEX))
             ? 'G' : ' ';
             A_code[l] = ((A_algn[l]==-2) &&
             (A_algn[l+szA] == -2))
             ? 'U' : A_code[l];

             A_code[l] = ((A_algn[l]==-2) &&
             (A_algn[l+szA] == -2))
             ? 'U' : A_code[l];

             A_code[l] = '
             } else {
             **** Shift Error Between ALignments ***
             A_code[l] = (A_algn[l]<0) ? 's' : 'S';
             }
            */
            A_code[l] = ((A_algn[l]>=0) && (A_algn[l+szA]>=0))
                ? 's' : 'D';
        }
    }
    
    /* Final fragment */
    if (f_len) {
        add_fragment(&frag_set, f_Sa, f_Sb, f_len);
        f_len=0;
    }
    
    /* Make conserved pairset */
    cnsvd = NULL;
    cnsvd = (struct pair_set *) E_MALLOC(sizeof(*cnsvd));
    if (r.n_correct)
        cons_set = E_REALLOC(cons_set,(sizeof(*cons_set)*r.n_correct));
    
    cnsvd->pairs = cons_set;
    cnsvd->n = r.n_correct;
    cnsvd->smpl_score = cnsvd->score = 0.0;
    *cons_pairs = cnsvd;
    if (r.n_correct == 0) {
        cnsvd->pairs = NULL;
        free(cons_set);
        free(cnsv);
    } else {
        cnsvd->pairs = cnsv+szA;
        smpl_score(cnsvd,Score_b);
        cnsvd->score = cnsvd->smpl_score;
        cnsvd->pairs = cnsv;
        smpl_score(cnsvd, Score_a);
        free(cnsv);
        cnsvd->pairs = cons_set;
    }
    
    for (m=0;m<szB; m++)
        if (B_algn[m]==B_algn[m+szB]) {
            if ((B_algn[m]==-2)||(B_algn[m]==GAP_INDEX)) {
                r.ng_correct++;
                B_code[m] = 'c';
            }
            else {
                r.n_correct++;
                B_code[m] = 'C';
            }
        } else {
            char c;
            /* Dont Care ? for the mo. - etc */
            c = (abs(B_algn[m]-B_algn[m+szB])<10)
                ? ((char) ((int) '0')+abs((B_algn[m]-B_algn[m+szB])))
                : 'S';
            B_code[m] = ((B_algn[m]>=0) && (B_algn[m+szB]>=0))
                ? c : ((B_algn[m]<0) ? ((B_algn[m]==-2) ? 'U' : 'I')
                       : ((B_algn[m+szB]== -2) ? 'u' : 'i'));
        }
    
    
    /* Make models and compute coverages */
    /* Models consist of common C residues from both structures */
    if (r.n_correct % 2) {
        err_printf(this_sub,
                   "confusion: count of consistent pairs is not a whole number!\nACode:%s\nBCode:%s\n",
                   A_code, B_code);
        return NULL;
    }
    model = E_MALLOC(sizeof(*model)*4);

    m = r.n_correct / 2; /* Size of Models */

    /* Could be clever and reuse make_model, but we just copy the code instead.*/
    /* Allocate space for coordinates, then the sequence */
    
    model[0] = model_fromAnnotAlign(A_code,A_algn, szA, Al_a, Sa, Sb, 0); /* Sa orig  */
    model[1] = (!m) ? model_dummy(Sa)
        : model_fromAnnotAlign(A_code,A_algn, szA, m, Sa, Sb, 1); /* Sa align consistent */
    model[2] = (!m) ? model_dummy(Sb)
        : model_fromAnnotAlign(B_code,B_algn, szB, m, Sb, Sa, 1); /* Sb align consistent */
    model[3] = model_fromAnnotAlign(B_code,B_algn, szB, Al_b, Sb, Sa, 0); /* Sb orig  */

    /* Ugly informationals  - not niced as yet */
    /*    mfprintf(stdout, "Correct %i\nC-Gaps %i\nA:\n%s\n",r.n_correct,
     *      (r.ng_correct/2), seq_print(Sa->seq));
     *
     * mfprintf(stdout, "%s\nB:\n%s\n%s\n", A_code, seq_print(Sb->seq), B_code);
     */

    free(A_code);
    free(B_code);
    free(A_algn);
    free(B_algn);
    
    *frag_string = frag_set;
    return model;
}

