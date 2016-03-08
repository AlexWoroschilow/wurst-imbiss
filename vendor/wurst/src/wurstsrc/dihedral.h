/*
 * 26 April 2001
 * This can only be #included after struct RPoint is defined, so
 * one has to #include "coord.h".
 * $Id: dihedral.h,v 1.2 2008/06/27 08:58:07 schenk Exp $
 */

#ifndef DIHEDRAL_H
#define DIHEDRAL_H

float
dihedral (const struct RPoint i, const struct RPoint j,
          const struct RPoint k, const struct RPoint l);

float
bondangle(const struct RPoint i, const struct RPoint j,
          const struct RPoint k);

#endif /*  DIHEDRAL_H */
