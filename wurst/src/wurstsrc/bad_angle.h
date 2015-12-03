/*
 * 25 Nov 2005
 * There is a protocol for saying that a dihedral angle cannot be
 * calculated. Flag it with this value.
 * rcsid = $Id: bad_angle.h,v 1.3 2013/07/23 09:55:28 mosisch Exp $
 */
#ifndef BAD_ANGLE_H
#define BAD_ANGLE_H
extern const float BAD_ANGLE_FLOAT;            /* some unlikely angle, initialization in common.c */
extern const float BAD_DISTANCE_FLOAT;         /* some unlikely distance, initialization in common.c */
extern const float BAD_ANGLE;                  /* some unlikely angle, initialization in common.c */
extern const float BAD_DISTANCE;               /* some unlikely distance, initialization in common.c */
#endif /* BAD_ANGLE_H */
