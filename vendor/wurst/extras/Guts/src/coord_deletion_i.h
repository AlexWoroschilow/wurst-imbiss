/* coord_deletion_i.h
 * 17-06-2004
   Simple deletions to structures
   J.B. Procter 2002 April (original main src branch)
*/

#ifndef COORD_DELETION_I_H
#define COORD_DELETION_I_H


struct seq *seq_deletion(struct seq *, size_t, size_t,size_t);
struct coord *coord_deletion(struct coord *, size_t, size_t,size_t);
struct coord *coord_segment( struct coord *c, size_t start, size_t end);
struct coord *coord_merge( struct coord *c1, struct coord *c2);

#endif
