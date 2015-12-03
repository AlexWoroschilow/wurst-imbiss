/*
 * 21 March 2001
 * Printf wrapper to allow us to control I/O
 */

#include <errno.h>
#include <stdarg.h>
#include <stdio.h>
#include <string.h>

#include "mprintf.h"

#if !defined (lint) && !defined (DONT_SEE_RCS)
    static const char *rcsid =
    "$Id: mprintf.c,v 1.1.1.1 2001/09/21 05:19:50 torda Exp $";
#endif /* !defined (lint) && !defined (DONT_SEE_RCS) */

/* ---------------- err_printf --------------------------------
 * This prints out the name of the programme, protein file name
 * and then the rest of the arguments to printf.
 * Interesting quirk.. If we are being called (under linux anyway),
 * the calls to Tcl bits, make a mess of "errno". This means that
 * a subsequent call to perror/mperror will produce rubbish. For
 * that reason, we store errno on function entry and take the
 * liberty of setting it before we return.
 */
int
err_printf (const char *fnc_name, const char *fmt, ...)
{
    int r;
    va_list ap;
    va_start (ap, fmt);
    mfprintf (stderr, "Function %s: ", fnc_name);
    r = vfprintf (stderr, fmt, ap);
    va_end (ap);
    return r;
}

/* ---------------- mprintf -----------------------------------
 * printf() wrapper which allows us to send output through Tcl's
 * Tcl_write routine.  Unfortunately, it is not a good idea to
 * store "outchannel" below.  If the user redirects output, we
 * want this to be picked up.
 */
int
mprintf (const char *fmt, ...)
{
    int ret;
    va_list ap;

    va_start (ap, fmt);
    ret = vfprintf (stdout, fmt, ap);
    va_end (ap);
    return ret;
}
/* ---------------- mputchar ----------------------------------
 * The tcl routine expects a char buffer, but putchar() type
 * routines operate on int's.  Hence the casting below.
 */
int
mputchar ( int c)
{
        return (putchar (c));
}

/* ---------------- mfprintf  ---------------------------------
 * Wrapper for fprintf.
 * If we are writing to stdout, use mprintf so our application
 * works under Tcl.  If we are writing to some other file
 * pointer, just call vfprintf to do the work.
 */
int
mfprintf ( FILE *fp, const char *fmt, ...)
{
    int r;
    va_list ap;
    va_start (ap, fmt);

    r = vfprintf (fp, fmt, ap);
    va_end (ap);
    return r;
}

/* ---------------- mfputc ------------------------------------
 * Wrapper for putc ().
 * If we are not writing to stdout, then don't interfere. Just
 * call standard putc().
 * If we are writing to stdout, then steal force use of Tcl
 * routines via mputchar.
 */
int
mfputc (int c, FILE *fp)
{
    if (fp != stdout)
        return (putc (c, fp));
    else
        return (mputchar (c));
}

/* ---------------- mputc  ------------------------------------
 * In the spirit of putc(), mputc() is defined as a macro.
 * See the mprintf.h header file.
 */


/* ---------------- mputs  ------------------------------------
 * Wrapper for puts.
 */
int mputs (const char *s)
{
    return (mprintf ("%s\n", s));
}
/* ---------------- mfputs ------------------------------------
 *
 */
int mfputs (const char *s, FILE *fp)
{
    if (fp == stdout)
        return ( mprintf ("%s", s));
    else
        return (fputs (s, fp));
}
/* ---------------- mperror -----------------------------------
 * This has the same calling interface as perror, but calls
 * err_printf() so one can treat any writes to stderr specially
 */
void
mperror (const char *s)
{
    err_printf (s, "%s\n", strerror (errno));
}
