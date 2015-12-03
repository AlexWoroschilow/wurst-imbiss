/*
 * 20 nov 2003
 *
 * Put commonly used strings here, rather than duplicate them in
 * lots of other files.
 * To use a string from here, #include common.h (maybe).
 * Really to use a string from here, just say
      extern const char *prog_bug
 * for example.
 *
 * $Id: common.c,v 1.3 2013/07/23 11:21:52 mosisch Exp $
 */

const char *prog_bug   = "Programming bug %s %d\n";
const char *null_point = "Called with null pointer\n";
const char *mismatch  ="Mismatch of score matrix with elements for scoring.\n\
Matrix set up for %d x %d,\n\
but items for scoring are %d x %d.\n\
Giving up.\n";
const float BAD_ANGLE          = -999.9;     /* some unlikely angle */
const float BAD_ANGLE_FLOAT    = -999.9;     /* some unlikely angle */
const float BAD_DISTANCE       = 0.0;     /* some unlikely distance */
const float BAD_DISTANCE_FLOAT = 0.0;     /* some unlikely distance */
