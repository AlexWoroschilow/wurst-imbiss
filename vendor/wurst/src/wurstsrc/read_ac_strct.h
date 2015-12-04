/*
 * 12 Dec 2005
 * rcsid = $Id: read_ac_strct.h,v 1.2 2014/06/18 11:32:00 torda Exp $
 */
#ifndef READ_AC_STRCT_H
#define READ_AC_STRCT_H

struct clssfcn;

struct aa_strct_clssfcn {
    struct clssfcn *strct;
    size_t n_att;
    float ***log_pp; /* log of attribute probability for residue */
};

#endif
