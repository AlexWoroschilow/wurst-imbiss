/*
 * 9 January 2002
 * $Id: scratch.h,v 1.1 2002/02/07 06:41:03 torda Exp $
 */

#ifndef SCRATCH_H
#define SCRATCH_H

void  scr_reset (void);
char *scr_printf (const char *fmt, ...)
#ifdef __GNUC__
__attribute__ ((format (printf, 1, 2)))
#endif /* __GNUC__ */
;
void  free_scratch (void);

#endif /* SCRATCH_H */
