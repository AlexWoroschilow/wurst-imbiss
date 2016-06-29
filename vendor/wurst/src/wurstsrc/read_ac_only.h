/*
 * 16 Feb 2012
 * rcsid = $Id: read_ac_only.h,v 1.2 2012/04/04 13:18:00 ibondarenko Exp $
 * There are functions which are used in two places for parsing
 * autoclass output. Define only those functions here.
 * This file can only be included after regex.h has been read.
 */
#ifndef READ_AC_ONLY_H
#define READ_AC_ONLY_H

int m_regcomp (regex_t *r_compiled, const char *regex);
char *find_line(char *buf, const int size, FILE *fp, const char *s);
char *
find_regex (char *buf, const int size, FILE *fp, const char *regex, int max);
size_t get_n_class (FILE *fp, char *buf, const int bufsiz);
float get_class_wt (const char *buf);


#endif /* READ_AC_ONLY_H */

