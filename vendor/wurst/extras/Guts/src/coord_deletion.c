/* coord_deletion.c
   Simple deletions to structures
   J.B. Procter 2002 April
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "e_malloc.h"
#include "mgc_num.h"
#include "mprintf.h"
#include "str.h"
#include "read_seq_i.h"
#include "seq.h"
#include "coord.h"
#include "coord_deletion_i.h"
#include "coord_i.h"
#include "seqprof_trim_i.h"

#if !defined (lint) && !defined (DONT_SEE_RCS)
 static const char *rcsid =
 "$Id: coord_deletion.c,v 1.1 2007/09/28 16:57:16 mmundry Exp $";
#endif /* !defined (lint) && !defined (DONT_SEE_RCS) */

struct coord* coord_merge ( struct coord *c1, struct coord *c2) {

    struct coord *nc, *cc;
    struct seq *nseq;
    size_t i;
    if (c1==NULL) 
        err_printf("coord_merge","c1 is undefined!\n");
    if (c2==NULL)
        err_printf("coord_merge","c2 is undefined!\n");
    
    nseq=seq_merge(c1->seq, c2->seq);
    if (!nseq) {
        err_printf("coord_merge","couldn't merge sequences\n");
        return(NULL);
    }

    nc=coord_template(c1, c1->size+c2->size);
    nc->seq = nseq;
    if (!nc) {
        err_printf("coord_merge","Couldn't make new coord object\n");
        return(NULL);
    }
    if (c1->units!=c2->units) {
        cc = coord_template(c2, c2->size);
        if (c1->units==nm) {
            coord_a_2_nm(cc);
        } else {
            if (c1->units==angstrom) {
                coord_nm_2_a(cc);
            }
        }
    } else {
        cc = c2;
    }
    /* merge c1 and cc 
       nc takes on all characteristics (chainid, etc)
    */
    nc->units=c1->units;
    for (i = 0; i < ACQ_SIZ; i++)
        nc->pdb_acq[i] = c1->pdb_acq[i];
    nc->chain=c1->chain; 
    /* and all coordinates present in both structures */
    
    if ((c1->sec_typ) && (c2->sec_typ)) {
        memcpy (nc->sec_typ, c1->sec_typ, c1->size * sizeof (*nc->sec_typ));
        memcpy (nc->sec_typ+c1->size, c2->sec_typ, c2->size * sizeof (*nc->sec_typ));
    }
    if (c1->rp_ca && c2->rp_ca) {
        memcpy (nc->rp_ca, c1->rp_ca, c1->size * sizeof (*nc->rp_ca));
        memcpy (nc->rp_ca+c1->size, c2->rp_ca, c2->size * sizeof (*nc->rp_ca));
    }
    if (c1->rp_cb && c2->rp_cb) {
        memcpy (nc->rp_cb, c1->rp_cb, c1->size * sizeof (*c1->rp_cb));
        memcpy (nc->rp_cb+c1->size, c2->rp_cb, c2->size * sizeof (*c2->rp_cb));
    }
    if (c1->rp_n && c2->rp_n) {
        memcpy (nc->rp_n, c1->rp_n, c1->size * sizeof (*c1->rp_n));
        memcpy (nc->rp_n+c1->size, c2->rp_n, c2->size * sizeof (*nc->rp_n));
    }
    if (c1->rp_c && c2->rp_c) {
        memcpy (nc->rp_c, c1->rp_c, c1->size * sizeof (*c1->rp_c));
        memcpy (nc->rp_c+c1->size, c2->rp_c, c2->size * sizeof (*nc->rp_c));
    }

    if (c1->rp_o && c2->rp_o) {
        memcpy (nc->rp_o, c1->rp_o, c1->size * sizeof (*nc->rp_o));
        memcpy (nc->rp_o+c1->size, c2->rp_o, c2->size*sizeof(*nc->rp_o));
    }
    if (c1->orig && c2->orig) {
        short ofset, *ncp,*nce;
        /* pain in the neck, here */
        memcpy (nc->orig, c1->orig, c1->size * sizeof (*c1->orig));
        memcpy(nc->orig+c1->size, c2->orig, c2->size*sizeof(*nc->orig));
        /* We really don't want to do this every time... 
        ofset = nc->orig[c1->size-1];
        ncp=nc->orig+c1->size;
        nce=nc->orig+nc->size;
        do {
            *ncp+=ofset;
        } while (++ncp<nce);
        */
    }
    if (c1->icode && c2->icode) {
        memcpy (nc->icode, c1->icode, c1->size * sizeof(*c1->icode));
        memcpy (nc->icode+c1->size, c2->icode, c2->size * sizeof(*c2->icode));
    }
    if (c1->psi && c2->psi) {
        nc->psi = E_MALLOC (nc->size * sizeof (*c1->psi));
        memcpy (nc->psi, c1->psi, c1->size * sizeof (*c1->psi));
        memcpy (nc->psi+c1->size, c2->psi, c2->size * sizeof (*c2->psi));
    }
        
    if (c2->units!=cc->units)
        coord_destroy(cc);
    return(nc);
}

struct coord *
coord_segment ( struct coord *c, size_t start, size_t end )
{
    struct coord *nc;
    struct seq *nseq;
    nc = NULL;
    nseq = NULL;
    if (c) {
        size_t n,i;
        /* range check */
        
        if (!end || (end<start)) 
            err_printf("coord_trim", 
                       "Broken script - end is < than start\n");
        
        if (start>=c->size)
            err_printf("coord_trim", 
                       "Broken script - start is beyond coord size (-ve position perhaps?)\n");
        
        if (end>c->size)  
            err_printf("coord_trim", 
                       "Broken script - end is beyond coord size\n");
        
        /* make new sequence */
        nseq = E_MALLOC(sizeof(*nseq));
        nseq->comment = NULL; /* we don't keep this */
        n = (nseq->length = end-start);
        nseq->seq = E_MALLOC(sizeof(*nseq->seq)*(nseq->length+1));
        nseq->format = c->seq->format;
        memcpy(nseq->seq, c->seq->seq+start, nseq->length);
        nseq->seq[nseq->length] = '\0';
        
        /* and new coord - following routines above */
        nc = coord_template (c, n);
        if (!nc) {
            err_printf("coord_trim", "---- NO COORDINAT SPACE ----");
            return(NULL);
        }
        nc->units=c->units;
        for (i = 0; i < ACQ_SIZ; i++)
            nc->pdb_acq[i] = c->pdb_acq[i];
        nc->chain=c->chain; /*not strictly true*/
        if (c->sec_typ)
            memcpy (nc->sec_typ, c->sec_typ+start, n * sizeof (*c->sec_typ));
        if (c->rp_ca)
            memcpy (nc->rp_ca, c->rp_ca+start, n * sizeof (*c->rp_ca));
        if (c->rp_cb)
            memcpy (nc->rp_cb, c->rp_cb+start, n * sizeof (*c->rp_cb));
        if (c->rp_n)
            memcpy (nc->rp_n, c->rp_n+start, n * sizeof (*c->rp_n));
        if (c->rp_c)
            memcpy (nc->rp_c, c->rp_c+start, n * sizeof (*c->rp_c));
        if (c->rp_o)
            memcpy (nc->rp_o, c->rp_o+start, n * sizeof (*c->rp_o));
        if (c->orig) {
            memcpy (nc->orig, c->orig+start, n * sizeof (*c->orig));
        }
        if (c->icode)
            memcpy (nc->icode, c->icode+start, n * sizeof(*c->icode));
        if (c->psi) {
            nc->psi = E_MALLOC (n * sizeof (*c->psi));
            memcpy (nc->psi, c->psi+start, n * sizeof (*c->psi));
        }
        nc->seq=nseq;
    }
    return (nc);
}
 
/* ---------------- seq_deletion ------------------------------
 * perl: seq_deletion seq start length
 *  returns  a new sequence with appropriately modified set
 *  Operation :
 *    for a deletion at s of length l,
 *    1.remove residues s to s+l-1 from the sequence
 */
struct seq *
seq_deletion ( struct seq *oseq, size_t start,
	       size_t sl, size_t ql)
{
    const char *this_sub="seq_deletion";
    /* From seq_copy */
    size_t n;
    struct seq *dst;

    if (oseq==NULL) {
        err_printf (this_sub, "---- NULL SEQUENCE ----\n");
        return(NULL);
    }
    if ((start>=oseq->length) ||
	((start+sl+ql)>oseq->length)) {
        err_printf (this_sub, "Unbounded Deletion of (%i, s%i, q%i) in sequence of %i residues.\n",start,sl,ql,oseq->length);
        return (NULL);
    }
    n = oseq->length-sl-ql;
    dst = E_MALLOC (sizeof (*oseq));
    if (dst!=NULL) {
        seq_ini (dst);
        if (oseq->seq) {
            dst->seq = E_MALLOC(sizeof(oseq->seq[0])*(n+1));
            if (start)
                memcpy(dst->seq, oseq->seq, (start)*sizeof(oseq->seq[0]));
            if (n-start)
                memcpy(dst->seq+start, oseq->seq+(start+ql), (1+n-start)*sizeof(oseq->seq[0]));
            dst->seq[n]='\0'; /* excessive ensureing */
            if (oseq->comment)
                dst->comment = save_str (oseq->comment);
        }
        dst->length = n;
        dst->format = oseq->format;
        return dst;
    } else
        err_printf (this_sub, "---- NO SEQUENCE MEMORY ----");
    return NULL;
}
/* ---------------- coord_deletion ----------------------------
 *
 * perl: coord_deletion coord start length
 *  returns  a new (coord) with appropriately modified atom and sequence set
 *  Operation :
 *    for a deletion at s of length l,
 *    1.remove residues s to s+l-1 from the sequence entry of the
 *    structure.
 *    2.remove coordinates at the end of the set for l residues
 */
struct coord *
coord_deletion ( struct coord *oc, size_t start, size_t sl, size_t ql)
{
    const char *this_sub = "coord_deletion";
    const char *unbound =
        "Unbounded Deletion of (%i, %i, %i) in structure of %i residues.\n";
    struct coord *nc;
    int i;
    size_t n;
    size_t fs,sc;

    /* Check for a sensible Deletion, or return
       undefined.
    */
    if (oc==NULL) {
        err_printf (this_sub, "---NULL--- coordinates\n");
        return NULL;
    }
    if ((start>=oc->size) || ((start+sl+ql)>oc->size)) {
        err_printf (this_sub, unbound, start, sl, ql, oc->size);
        return NULL;
    }

    /* Make a new coord of appropriate size
       and copy over the deletioned sequence and
       atomset.
    */
    n = oc->size - (ql + sl); /* Final Pairset Size */
    nc = coord_template (oc, n);
    if (! nc) {
        err_printf(this_sub, "---- NO COORDINATE SPACE ----");
        return(NULL);
    }

    nc->units = oc->units;
    for (i = 0; i < ACQ_SIZ; i++)
        nc->pdb_acq[i] = oc->pdb_acq[i];
    nc->chain = oc->chain;
    fs = start;
    sc = ((start + sl + ql) < (oc->size)) ? (n - start) : 0;
    if (sl) {
        size_t slstart = (start + sl);
        if (oc->sec_typ) {
            if (fs)
                memcpy (nc->sec_typ, oc->sec_typ, start*sizeof (*oc->sec_typ));
            if (sc)
                memcpy (nc->sec_typ+start, oc->sec_typ+slstart, sc*sizeof (*oc->sec_typ));
        }
        if (oc->rp_ca) {
            if (fs)
                memcpy (nc->rp_ca, oc->rp_ca, start*sizeof (*oc->rp_ca));
            if (sc)
                memcpy (nc->rp_ca+start, oc->rp_ca+slstart, sc*sizeof (*oc->rp_ca));
        }
        if (oc->rp_cb) {
            if (fs)
                memcpy (nc->rp_cb, oc->rp_cb, start*sizeof (*oc->rp_cb));
            if ((sc))
                memcpy (nc->rp_cb+start, oc->rp_cb+slstart, sc*sizeof (*oc->rp_cb));
        }
        if (oc->rp_n) {
            if (fs)
                memcpy (nc->rp_n, oc->rp_n, start*sizeof (*oc->rp_n));
            if ((sc))
                memcpy (nc->rp_n+start, oc->rp_n+slstart, sc*sizeof (*oc->rp_n));
        }
        if (oc->rp_c) {
            if (fs)
                memcpy (nc->rp_c, oc->rp_c, start*sizeof (*oc->rp_c));
            if ((sc))
                memcpy (nc->rp_c+start, oc->rp_c+slstart, sc*sizeof (*oc->rp_c));
        }
        if (oc->rp_o) {
            if (fs)
                memcpy (nc->rp_o, oc->rp_o, start*sizeof (*oc->rp_o));
            if ((sc))
                memcpy (nc->rp_o+start, oc->rp_o+slstart, sc*sizeof (*oc->rp_o));
        }
        if (oc->icode) {
            if (fs)
                memcpy (nc->icode, oc->icode, start*sizeof (*oc->icode));
            if ((sc))
                memcpy (nc->icode+start, oc->icode+slstart, sc*sizeof (*oc->icode));
        }
        if (oc->psi) {
            nc->psi = E_MALLOC(n * sizeof(*nc->psi));
            if (fs)
                memcpy (nc->psi, oc->psi, start*sizeof (*oc->psi));
            if ((sc))
                memcpy (nc->psi+start, oc->psi+slstart, sc*sizeof (*oc->psi));
        }
        if (oc->orig)
            memcpy (nc->orig, oc->orig, n * sizeof (*oc->orig));
        
    } else { /* ! sl , Simple coordinate set truncation */
        if (oc->sec_typ)
            memcpy (nc->sec_typ, oc->sec_typ, n * sizeof (*oc->sec_typ));
        if (oc->rp_ca)
            memcpy (nc->rp_ca, oc->rp_ca, n * sizeof (*oc->rp_ca));
        if (oc->rp_cb)
            memcpy (nc->rp_cb, oc->rp_cb, n * sizeof (*oc->rp_cb));
        if (oc->rp_n)
            memcpy (nc->rp_n, oc->rp_n, n * sizeof (*oc->rp_n));
        if (oc->rp_c)
            memcpy (nc->rp_c, oc->rp_c, n * sizeof (*oc->rp_c));
        if (oc->rp_o)
            memcpy (nc->rp_o, oc->rp_o, n * sizeof (*oc->rp_o));
        if (oc->orig) {
            if (fs)
                memcpy (nc->orig, oc->orig, start * sizeof (*oc->orig));
            if ((sc))
                memcpy (nc->orig+start, oc->orig+(ql+start), (n-start) * sizeof (*oc->orig));
        }
        if (oc->icode)
            memcpy (nc->icode, oc->icode, n * sizeof(*oc->icode));
        if (oc->psi) {
            nc->psi = E_MALLOC (n * sizeof (*nc->psi));
            memcpy (nc->psi, oc->psi, n * sizeof (*oc->psi));
        }
    }
    /* Write Sequence */
    nc->seq = seq_deletion(oc->seq, start, sl,ql);
    return nc;
}
