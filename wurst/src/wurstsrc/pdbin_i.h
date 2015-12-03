/*
 * 23 March 2002
 * Interface to pdb reading routines.
 * rcsid = $Id: pdbin_i.h,v 1.2 2012/02/14 15:26:58 ploeffler Exp $
 */
#ifndef PDBIN_I_H
#define PDBIN_I_H

struct coord;
struct coord *
pdb_read ( const char *fname, const char *acq_c, const char chain);
void check_cmpnd ( struct coord *c);

#endif /* PDBIN_I_H */
