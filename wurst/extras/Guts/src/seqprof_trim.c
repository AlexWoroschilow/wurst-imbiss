/* seqprof_trim.c
 * 8-7-2004
 * take a range of values from a sequence profile
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "amino_a.h"
#include "e_malloc.h"
#include "fio.h"
#include "matrix.h"
#include "mprintf.h"
#include "read_blst.h"
#include "read_seq_i.h"
#include "scratch.h"
#include "seq.h"
#include "seqprof.h"

#include "seqprof_trim_i.h"
/* ------------------- seqprof_trim ---------------------------------------
 * start and end are in the range 1...seq_length
 */
struct seqprof *
seqprof_trim ( struct seqprof *prf, size_t start, size_t end) {
    struct seq *nseq;
    struct seqprof *nprf;
    nprf = NULL;
    if (prf) {
        /* range check */
        if (!start) 
            err_printf("seqprof_trim", 
                       "Broken script - start,end take range from 1->seq_length\n");
        
        if (end && (end<start)) 
            err_printf("seqprof_trim", 
                       "Broken script - end is < than start\n");
        
        if (start>prf->seq->length)
            err_printf("seqprof_trim", 
                       "Broken script - start is beyond length\n");
        
        if (end>prf->seq->length)  
            err_printf("seqprof_trim", 
                       "Broken script - end is beyond length\n");
        
        /* make new sequence */
        nseq = E_MALLOC(sizeof(*nseq));
        nseq->comment = NULL; /* we don't keep this */
        nseq->length = end-start+1;
        nseq->seq = E_MALLOC(sizeof(*nseq->seq)*(nseq->length+1));
        nseq->format = prf->seq->format;
        memcpy(nseq->seq, prf->seq->seq+start-1, nseq->length);
        nseq->seq[nseq->length] = '\0';
        
        /* and new profile - following matrix.c:f_matrix v 1.4 */
        nprf = E_MALLOC(sizeof(*nprf));
        nprf->seq = nseq;
        nprf->nres = nseq->length;
        nprf->freq_mat = copy_f_matrix(&prf->freq_mat[start-1], 
                                       end+1-start, MAX_AA);
    }
    return nprf;
}
/* -------------------- seqprof_merge --------------------------------
 * merge two sequence profiles, sequentially.
 */

struct seq *seq_merge ( struct seq *seq1, struct seq *seq2) {
    struct seq *nseq;
    if (!seq1 || !seq2) {
        err_printf("seq_merge","Given null sequences to merge\n");
    }
    nseq = E_MALLOC(sizeof(*nseq));
    nseq->comment = NULL; /* we don't keep this */
    nseq->length = seq1->length+seq2->length;
    nseq->seq = E_MALLOC(sizeof(*nseq->seq)*(nseq->length+1));
    nseq->format = seq1->format;
    memcpy(nseq->seq, seq1->seq,seq1->length);
    if (seq2->format!=nseq->format) {
        /* convert formats */
        struct seq *tseq = seq_copy(seq2);
        if (nseq->format==PRINTABLE) {
            seq_thomas2std(tseq);
        } else {
            if (nseq->format==THOMAS) {
                seq_std2thomas(tseq);
            } else {
                err_printf("seqprof_merge","source inconsistency - unknown sequence format!\n");
            }
        }
        memcpy(nseq->seq+seq1->length, tseq->seq, tseq->length);
        seq_destroy(tseq);

    } else {
        memcpy(nseq->seq+seq1->length, seq2->seq, sizeof(*seq2->seq)*seq2->length);
    }
    nseq->seq[nseq->length] = '\0';
    return(nseq);
}

struct seqprof *
seqprof_merge ( struct seqprof *prf1, struct seqprof *prf2) {
    struct seq *nseq;
    struct seqprof *nprf;
    size_t i,j;
    nprf = NULL;
    if (!prf1) 
        err_printf("seqprof_merge", "prf1 is not defined.\n");
    if (!prf2) 
        err_printf("seqprof_merge", "prf2 is not defined.\n");

    /* make new sequence */
    nseq = seq_merge(prf1->seq, prf2->seq);    
    /* and new profile - following matrix.c:f_matrix v 1.4 */
    nprf = E_MALLOC(sizeof(*nprf));
    nprf->seq = nseq;
    nprf->nres = nseq->length;
    nprf->freq_mat = f_matrix(prf1->nres+prf2->nres, MAX_AA);
    memcpy(nprf->freq_mat[0], prf1->freq_mat[0], sizeof(prf1->freq_mat[0][0])*prf1->nres*MAX_AA);
    memcpy(nprf->freq_mat[prf1->nres], prf2->freq_mat[0], sizeof(prf2->freq_mat[0][0])*prf2->nres*MAX_AA); 

    return nprf;
}
