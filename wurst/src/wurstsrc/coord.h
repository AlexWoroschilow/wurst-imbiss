/* 23 Oct 2001
 * Structures for coordinates.
 * This does *not* define the interface to the routines.
 * It may only be included by routines which need to manipulate
 * the innards of a coordinate structure.
 *
 * This may only be included after "units.h".
 * rcsid = $Id: coord.h,v 1.3 2014/07/04 10:54:46 torda Exp $
 */
#ifndef COORD_H
#define COORD_H

#include "stdlib.h"

enum { ACQ_SIZ = 5 };  /* Length of PDB acquisition code */

enum units {         /* Do we have coordinates in Angstrom or nanometers ? */
    nm,
    angstrom,
    broken
};

struct RPoint {
    float x;
    float y;
    float z;
};

struct seq;
struct coord  {
    struct RPoint *rp_ca;     /* Alpha carbon coordinates */
    struct RPoint *rp_cb;     /* Beta  carbon coordinates */
    struct RPoint *rp_n;      /* Nitrogen coordinates */
    struct RPoint *rp_c;      /* Carbonyl coordinates */
    struct RPoint *rp_o;      /* Carbonyl coordinates */
    short *orig;              /* Original sequence numbering in PDB file */
    char *icode;              /* Insertion codes */
    float *omega;             /* For storing calculated dihedral angles */
    float *psi;               /* For storing calculated dihedral angles */
    float *phi;               /* For storing calculated dihedral angles */
    float *theta;             /* For storing calculated pseudo bond angles */
    float *tau;               /* For storing calculated pseudo dih. angles */
    int *sec_typ;             /* Secondary struct info */
    char *compnd;             /* Compound information  */
    struct seq *seq;          /* Sequence structure */
    size_t size;
    size_t compnd_len;
    enum units units;
    char pdb_acq[ACQ_SIZ];    /* PDB acquisition code */
    char chain;               /* Chain identifier or _ for no specific chain */
};



#endif  /* COORD_H */
