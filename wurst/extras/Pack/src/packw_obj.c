/* 21st August 2002
 * Routines for packing sequences, coordinates,
 * secondary structure predictions and sequence
 * profiles.

 * Currently untyped (ie we don't really check if we are
 * given garbage to unpack).  Typically, this results
 * in an out of memory error, but there is a simple
 * checksum.
 *
 * Useful for perl based storage, retrieval
 * and transfer of binary coordinate/sequence data in
 * a homogenous environment.
 * 
 * Depends on units.h, coord.h, seq.h, sec_s.h, and
 * seqprof.h
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "e_malloc.h"
#include "mprintf.h"
#include "coord.h"
#include "amino_a.h"
#include "seq.h"
#include "coord_i.h"
#include "sec_s.h"
#include "seqprof.h"
#include "read_blst.h"
#include "packw_obj_i.h"

#if !defined (lint) && !defined (DONT_SEE_RCS)
static const char *rcsid =
    "$Id: packw_obj.c,v 1.1 2007/09/28 16:57:20 mmundry Exp $";
#endif
/* !defined (lint) && !defined (DONT_SEE_RCS) */



/* ---------------- pack_seq -------------------------------------------
 * Pack a sequence structure into an allocated
 * block of memory
 */

static void *pack_seq(struct seq *seq)
{
    struct seq *seqp;
    char *cp;
    size_t tot_size, c_size, *tp;
    char *block;

    block = NULL;
    if (seq != NULL) {
        tot_size = sizeof(size_t) + sizeof(struct seq)
            + ((seq->comment != NULL) ?
               (c_size = (sizeof(char) * (1 + strlen(seq->comment)))) : 0)
            +
            ((seq->seq != NULL) ? (sizeof(char) * (1 + seq->length)) : 0);
        block = E_MALLOC(tot_size);
        tp = (size_t *) block;
        seqp = (struct seq *) (tp + 1);
        cp = (char *) (seqp + 1);
        *tp = tot_size;
        seqp->length = seq->length;
        seqp->format = seq->format;
        if (seq->seq) {
            seqp->seq = (char *) (((size_t) cp) - (size_t) block);
            memcpy(cp, seq->seq, 1 + seq->length);
            cp += 1 + seq->length;
        } else {
            seqp->seq = NULL;
        }

        if (seq->comment) {
            seqp->comment = (char *) (((size_t) cp) - (size_t) block);
            memcpy(cp, seq->comment, c_size);
        } else {
            seqp->comment = NULL;
        }

    }
    return (void *) block;
}


/* -------- reloc_pseq -------------------------------------------------
 * Relocate the pointers in the void block, returning the 
 * struct_seq pointer 
 */

static struct seq *reloc_pseq(void *sq)
{
    const char *tfn = "reloc_pseq";
    struct seq *sp;
    char *cpe;
    size_t base, bend;
    sp = NULL;
    if (sq != NULL) {
        base = (size_t) sq;
        bend = base + *(size_t *) sq;
        cpe = (char *) bend;
        sp = (struct seq *) (((size_t *) sq) + 1);

        /* Assertions */
        if (sp->seq != NULL) {
            if (((base + (size_t) sp->seq) + sp->length + 1) <= bend) {
                sp->seq = (char *) (base + (size_t) sp->seq);
            } else {
                err_printf(tfn, "seq pointer is out of range of block!\n");
                return (NULL);
            }
        } else {
            if (sp->length != 0) {
                err_printf(tfn,
                           "length is non-zero for null sequence!\n|");
                err_printf(tfn, "(probable corruption)\n");
                return (NULL);
            }
        }

        if (sp->comment != NULL) {
            if ((base + (size_t) sp->comment) < (bend - 1)) {
                if (*(cpe - 1) == '\0') {
                    sp->comment = (char *) (base + (size_t) sp->comment);
                } else {
                    err_printf(tfn, "*cpe is non-zero character\n");
                    return (NULL);
                }
            } else {
                err_printf(tfn,
                           "Comment pointer is out of range of block!\n");
                return (NULL);
            }
        }

    } else {
        err_printf(tfn, "Received a null sequence block pointer\n");
    }
    return sp;
}

/* ------coord_pack ----------------------------------------------------
 * Pack a coord structure into an allocated
 * block of memory
 */
void *coord_pack(struct coord *coord)
{
    char *block, *seq_block;
    size_t tot_size, i;
    struct coord *cb;
    struct RPoint *atpos;
    char *cp;
    float *flp;
    short *sh_p;
    int *sec_p;
    block = NULL;
    seq_block = NULL;
    if (coord != NULL) {

        if (coord->seq != NULL)
            seq_block = pack_seq(coord->seq);

        tot_size = 2 * sizeof(size_t)
            + ((seq_block != NULL) ? *(size_t *) seq_block : 0)
            + sizeof(*coord)
            +
            ((coord->rp_n !=
              NULL) ? sizeof(*coord->rp_n) * coord->size : 0) +
            ((coord->rp_c !=
              NULL) ? sizeof(*coord->rp_c) * coord->size : 0) +
            ((coord->rp_o !=
              NULL) ? sizeof(*coord->rp_o) * coord->size : 0) +
            ((coord->rp_ca !=
              NULL) ? sizeof(*coord->rp_ca) * coord->size : 0) +
            ((coord->rp_cb !=
              NULL) ? sizeof(*coord->rp_cb) * coord->size : 0) +
            ((coord->psi !=
              NULL) ? sizeof(*coord->psi) * coord->size : 0) +
            ((coord->sec_typ !=
              NULL) ? sizeof(*coord->sec_typ) * coord->size : 0) +
            ((coord->icode !=
              NULL) ? sizeof(*coord->icode) * coord->size : 0) +
            ((coord->orig !=
              NULL) ? sizeof(*coord->orig) * coord->size : 0);

        block = E_MALLOC(tot_size);

        *(size_t *) block = tot_size;   /* Next, the 'checksum' */
        *(1 + (size_t *) block) = tot_size - sizeof(struct coord);
        cb = (struct coord *) (((size_t *) block) + 2);
        atpos = (struct RPoint *) (&cb[1]);
        cb->size = coord->size;
        cb->units = coord->units;
        cb->chain = coord->chain;
        for (i = 0; i < ACQ_SIZ; i++)
            cb->pdb_acq[i] = coord->pdb_acq[i];

        if (coord->rp_n) {
            cb->rp_n =
                (struct RPoint *) (((size_t) atpos) - (size_t) block);
            memcpy(atpos, coord->rp_n, sizeof(*coord->rp_n) * coord->size);
            atpos += coord->size;
        } else {
            cb->rp_n = NULL;
        }

        if (coord->rp_ca) {
            cb->rp_ca =
                (struct RPoint *) (((size_t) atpos) - (size_t) block);
            memcpy(atpos, coord->rp_ca,
                   sizeof(*coord->rp_ca) * coord->size);
            atpos += coord->size;
        } else {
            cb->rp_ca = NULL;
        }

        if (coord->rp_cb) {
            cb->rp_cb =
                (struct RPoint *) (((size_t) atpos) - (size_t) block);
            memcpy(atpos, coord->rp_cb,
                   sizeof(*coord->rp_cb) * coord->size);
            atpos += coord->size;
        } else {
            cb->rp_cb = NULL;
        }

        if (coord->rp_c) {
            cb->rp_c =
                (struct RPoint *) (((size_t) atpos) - (size_t) block);
            memcpy(atpos, coord->rp_c, sizeof(*coord->rp_c) * coord->size);
            atpos += coord->size;
        } else {
            cb->rp_c = NULL;
        }

        if (coord->rp_o) {
            cb->rp_o =
                (struct RPoint *) (((size_t) atpos) - (size_t) block);
            memcpy(atpos, coord->rp_o, sizeof(*coord->rp_o) * coord->size);
            atpos += coord->size;
        } else {
            cb->rp_o = NULL;
        }

        flp = (float *) atpos;
        if (coord->psi) {
            cb->psi = (float *) (((size_t) flp) - (size_t) block);
            memcpy(flp, coord->psi, sizeof(float) * coord->size);
            flp += coord->size;
        } else {
            cb->psi = NULL;
        }

        sec_p = (int *) flp;
        if (coord->sec_typ) {
            cb->sec_typ = (int *) (((size_t) sec_p) - (size_t) block);
            memcpy(sec_p, coord->sec_typ, sizeof(int) * coord->size);
            sec_p += coord->size;
        } else {
            cb->sec_typ = NULL;
        }

        cp = (char *) sec_p;
        if (coord->seq) {
            /* this points to the relative position of the
             * packed sequence structure */
            cb->seq =
                (struct seq *) (sizeof(size_t) + ((size_t) cp) -
                                ((size_t) block));
            memcpy(cp, seq_block, *(size_t *) seq_block);
            cp = (char *) (((size_t) cp) + *(size_t *) seq_block);
            free(seq_block);
        } else {
            cb->seq = NULL;
        }

        sh_p = (short *) cp;
        if (coord->orig) {
            cb->orig = (short *) (((size_t) sh_p) - (size_t) block);
            memcpy(sh_p, coord->orig, sizeof(short) * coord->size);
            sh_p += coord->size;
        } else {
            cb->orig = NULL;
        }

        cp = (char *) sh_p;
        if (coord->icode) {
            cb->icode = (char *) (((size_t) cp) - (size_t) block);
            memcpy(cp, coord->icode, sizeof(char) * coord->size);
            cp += coord->size;
        } else {
            cb->icode = NULL;
        }
    }
    return (void *) block;
}

/* ---------------reloc_pcoord------------------------------------------
 * relocate the internal pointers of struct coord
 */

static struct coord *reloc_pcoord(void *cblk)
{
    const char *tfn = "reloc_pcoord";
    struct coord *c;
    size_t base, bend;
    c = NULL;
    if (cblk != NULL) {
        base = (size_t) cblk;
        bend = base + *(size_t *) cblk;
        c = (struct coord *) (((size_t *) cblk) + 2);
        if ((base + 1 + sizeof(*c)) >= bend) {
            err_printf(tfn, "Block too small for coord structure\n");
        } else {

            if (c->seq) {
                if ((c->seq = reloc_pseq((char *)
                                         (base - sizeof(size_t)
                                          + (size_t) c->seq))) == NULL)
                    err_printf(tfn, "reloc_pseq failed.\n");
            }

            if (c->rp_n) {
                if ((base + c->size * sizeof(struct RPoint) +
                     (size_t) c->rp_n) < bend)
                    c->rp_n = (struct RPoint *) (base + (size_t) c->rp_n);
                else
                    err_printf(tfn, "rp_n cannot be relocated.\n");
            }

            if (c->rp_ca) {
                if ((base + c->size * sizeof(struct RPoint) +
                     (size_t) c->rp_ca) < bend)
                    c->rp_ca =
                        (struct RPoint *) (base + (size_t) c->rp_ca);
                else
                    err_printf(tfn, "rp_ca cannot be relocated.\n");
            }

            if (c->rp_cb) {
                if ((base + c->size * sizeof(struct RPoint) +
                     (size_t) c->rp_cb) < bend)
                    c->rp_cb =
                        (struct RPoint *) (base + (size_t) c->rp_cb);
                else
                    err_printf(tfn, "rp_cb cannot be relocated.\n");
            }

            if (c->rp_c) {
                if ((base + c->size * sizeof(struct RPoint) +
                     (size_t) c->rp_c) < bend)
                    c->rp_c = (struct RPoint *) (base + (size_t) c->rp_c);
                else
                    err_printf(tfn, "rp_c cannot be relocated.\n");
            }

            if (c->rp_o) {
                if ((base + c->size * sizeof(struct RPoint) +
                     (size_t) c->rp_o) < bend)
                    c->rp_o = (struct RPoint *) (base + (size_t) c->rp_o);
                else
                    err_printf(tfn, "rp_o cannot be relocated.\n");
            }

            if (c->orig) {
                if ((base + c->size * sizeof(*c->orig) +
                     (size_t) c->orig) < bend)
                    c->orig = (short *) (base + (size_t) c->orig);
                else
                    err_printf(tfn, "orig cannot be relocated.\n");
            }

            if (c->psi) {
                if ((base + c->size * sizeof(*c->psi) +
                     (size_t) c->psi) < bend)
                    c->psi = (float *) (base + (size_t) c->psi);
                else
                    err_printf(tfn, "psi cannot be relocated.\n");
            }

            if (c->sec_typ) {
                if ((base + c->size * sizeof(*c->sec_typ) +
                     (size_t) c->sec_typ) < bend)
                    c->sec_typ = (int *) (base + (size_t) c->sec_typ);
                else
                    err_printf(tfn, "sec_typ cannot be relocated.\n");
            }

            if (c->icode) {
                if ((base + c->size * sizeof(*c->icode) +
                     (size_t) c->icode) <= bend)
                    c->icode = (char *) (base + (size_t) c->icode);
                else
                    err_printf(tfn, "icode cannot be relocated.\n");
            }
        }
    } else {
        err_printf(tfn, "Null block!\n");
    }
    return c;
}

/* --------coord_unpack ------------------------------------------------
 * 
 * We can use the 'copy constructors' for
 * coord on the packed objects, but only after we have
 * shifted the pointers to their absolute addresses
 * in the packed object
 */

extern struct coord *coord_template(struct coord *, size_t);

struct coord *coord_unpack(void *pak_coord)
{
    const char *tfn = "coord_unpack";
    struct coord *mycoord, *new_coord;
    char *myb;
    size_t i;

    if (pak_coord != NULL) {
        size_t *v;
        v = (size_t *) pak_coord;
        /* 'Checksum' */
        if ((v[0] - v[1]) == sizeof(struct coord)) {

            myb = E_MALLOC(*(size_t *) pak_coord);

            memcpy(myb, pak_coord, *(size_t *) pak_coord);
            mycoord = NULL;
            new_coord = NULL;

            if ((mycoord = reloc_pcoord(myb))
                && (new_coord = coord_template(mycoord, mycoord->size))) {
                new_coord->seq = NULL;

                if (mycoord->seq)
                    new_coord->seq = coord_get_seq(mycoord);

                for (i = 0; i < ACQ_SIZ; i++)
                    new_coord->pdb_acq[i] = mycoord->pdb_acq[i];

                new_coord->chain = mycoord->chain;
            }
            free(myb);
            return new_coord;
        } else {
            /* Too lazy to try to flip endianism before giving up */
            err_printf(tfn, "Silly checksum failed.\n");
            err_printf(tfn, "(Not a packed coord structure?)\n");
        }
    } else {
        err_printf(tfn, "Given a null pointer to unpack\n");
    }
    return (NULL);
}

/* -------------------- pack_sec_s   -------------------------- 
 * Does the same as for the other opaque datatypes above
 * There are no checks, however, beyond simple buffer overruns.
 */

void *sec_s_pack(struct sec_s_data *msdat)
{
    char *pb;
    size_t i;

    pb = E_MALLOC(i = (2 * sizeof(i) + (sizeof(*msdat->data) * msdat->n)));
    *(size_t *) pb = i;
    *((size_t *) pb + 1) = msdat->n;
    memcpy((char *) ((size_t *) pb + 2), msdat->data,
           msdat->n * (sizeof(*msdat->data)));
    return ((void *) pb);
}

/* -------------------- unpack_sec_s -------------------------- 
 *
 */
struct sec_s_data *sec_s_unpack(void *pb)
{
    const char *tfn = "sec_s_unpack";
    struct sec_s_data *pdat, t;
    size_t *p;
    p = (size_t *) pb;
    pdat = NULL;
    if (pb) {
        if (p[0] != (2 * sizeof(*p) + p[1] * sizeof(*t.data))) {
            err_printf(tfn, "Silly simple checksum failed. Giving up!\n");
            err_printf(tfn, "(Not a packed sec_s_data structure ?)\n");
        } else {
            pdat = E_MALLOC(sizeof(struct sec_s_data));
            pdat->data = E_MALLOC(sizeof(*pdat->data) * p[1]);
            memcpy(pdat->data, p + 2, sizeof(*pdat->data) * p[1]);
            pdat->n = p[1];
        }
    } else
        err_printf(tfn, "Given a void pointer to unpack\n");

    return (pdat);
}

/* ---------------------- seqprof_pack ------------------------
 */

void *seqprof_pack(struct seqprof *prof)
{
    char *pp, *p, *pseq;
    size_t j, *v;
    if (prof != NULL) {
        if (prof->seq != NULL)
            pseq = pack_seq(prof->seq);
        pp = E_MALLOC(j = ((3 * sizeof(size_t) +
                            ((pseq != NULL) ? *((size_t *) pseq) : 1)
                            +
                            sizeof(**prof->freq_mat) * prof->nres *
                            MAX_AA)));
        v = (size_t *) pp;
        v[0] = j;
        v[1] = j - sizeof(*prof);
        v[2] = prof->nres;
        if (pseq == NULL) {
            v[3] = 1;           /* length of the block for null sequence */
            p = (char *) v + 4;
        } else {
            memcpy(v + 3, pseq, *(size_t *) pseq);
            p = (char *) ((size_t) (v + 3)) + *(size_t *) pseq;
        }
        /* block copy the f_matrix - following matrix.c */
        memcpy(p, prof->freq_mat[0],
               sizeof(*prof->freq_mat[0]) * prof->nres * MAX_AA);
    } else
        pp = NULL;
    return (void *) pp;
}

/* ----------------------seqprof_unpack---------------------------------
 */
struct seqprof *seqprof_unpack(void *pp)
{
    const char *tfn = "seqprf_unpack";
    struct seqprof *prof;
    size_t *v, i;

    if (pp != NULL) {

        v = (size_t *) pp;

        if (sizeof(struct seqprof) == (v[0] - v[1])) {

            prof = E_MALLOC(sizeof(struct seqprof));
            prof->nres = v[2];

            if (v[3] == 1) {
                prof->seq = NULL;
            } else {
                /* nasty self-reference copy here 
                   Expect it to fail if seqprof_get_seq 
                   has any assertions about the validity of
                   *prof
                 */
                if ((prof->seq = reloc_pseq((char *) (v + 3))) == NULL)
                    err_printf(tfn,
                               "reloc_pseq failed(packed seqprof?)\n");

                prof->seq = seqprof_get_seq(prof);
            }

            if ((v[2] * MAX_AA * sizeof(**prof->freq_mat)
                 + v[3] + 3 * sizeof(size_t)) != v[0]) {
                err_printf(tfn,
                           "Size of freq_mat and Sequence don't match packed block\n");
            } else {
                prof->freq_mat = E_MALLOC(sizeof(*prof->freq_mat) * v[2]);
                prof->freq_mat[0] =
                    E_MALLOC(sizeof(**prof->freq_mat) * v[2] * MAX_AA);

                memcpy(prof->freq_mat[0],
                       (char *) (((size_t) (v + 3)) + v[3]),
                       sizeof(**prof->freq_mat) * v[2] * MAX_AA);

                for (i = 1; i < v[2]; i++)
                    prof->freq_mat[i] = prof->freq_mat[i - 1] + MAX_AA;
            }
        } else {
            err_printf(tfn, "Simple checksum failed.\n");
        }
    } else {
        err_printf(tfn, "Given a void pointer to unpack\n");
    }
    return (prof);
}
