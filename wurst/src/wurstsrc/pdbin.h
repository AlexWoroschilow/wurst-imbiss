/*
 * 10 Juni 2014
 * Interface to pdb internal structeres and routines.
 * Needed for the integration into ost project 
 */
#ifndef PDBIN_H
#define PDBIN_H

/* ---------------- Constants   -------------------------------------
 * We have a bunch of switches in the code. These could /
 * should be passed in. For the moment, we set them here.
 */

static const float CC_BOND = 1.54;

/* 28 May 96
 * Format for ATOM records in a PDB file.
 */

enum {
    SERIAL =       7,      /* int, atom serial number */
    SERIAL_LEN =   5,
    AT_NAM =      12,      /* char *, atom name */
    AT_NAM_LEN =   4,
    ALT_LOC =     16,      /* alternate location indicator */
    ALT_LOC_LEN =  1,
    RES_NAM =     17,
    RES_NAM_LEN =  3,
    CHAIN_ID =    21,
    CHAIN_ID_LEN = 1,
    RES_SEQ =     22,
    RES_SEQ_LEN =  4,
    ICODE =       26,      /* Code for insertion of residues */
    ICODE_LEN =    1,
    X_START =     30,
    X_LEN =        8,
    Y_START =     38,
    Z_START =     46,
    OCC =         54,
    OCC_LEN =      6,
    BFAC =        60,
    BFAC_LEN =     6,
    SEG_ID =      72,      /* Segment identifier, left-justified */
    SEG_ID_LEN =   4,
    ELEM =        76,     /* Element symbol, right-justified */
    ELEM_LEN =     2,
    CHARGE =      78,
    CHARGE_LEN =   2
};



/* ---------------- Structures  -------------------------------------
 */
/* This is really a copy of an ATOM record from a PDB file.
 */
struct atompdb {
    int at_num;
    char at_name[5];
    char alt_loc;
    char res_name[4];
    char chain_id;
    int res_num;
    char icode;
    float x, y, z;
    enum yes_no wanted;
    char error;
};

struct allatoms {
    struct atompdb *atoms;
    size_t n;
};

struct bin_res { /* Bare number of atoms to make a residue worthwhile */
    struct RPoint rp_ca,
                  rp_cb,
                  rp_n,
                  rp_c,
                  rp_o;
};

enum at_type {
    BORING,
    CA,          /* A list of atom types we save and write out */
    CB,
    N,
    C,
    O
};

struct atom_ok {
    unsigned int ca : 1;
    unsigned int cb : 1;
    unsigned int n  : 1;
    unsigned int c  : 1;
    unsigned int o  : 1;
};

struct atompdb 
nothing_atom ( void );

struct atompdb
line2atom (char *inbuf);

struct coord *
atoms2mdl (struct allatoms allatoms);
#endif /* PDBIN_H */