/*
 * 22 March 2001
 * rcsid = $Id: fio.h,v 1.3 2006/03/21 21:28:29 torda Exp $
 */
#ifndef FIO_H
#define FIO_H

FILE *mfopen (const char *fname, const char *mode, const char *s);
int   file_no_cache (FILE *fp);
int   file_clear_cache (FILE *fp);
#endif /* FIO_H */
