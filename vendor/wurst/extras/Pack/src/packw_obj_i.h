/* 17th June 2004
 *  Public interface to routines in packw_obj.c
 * rcsid = $Id: packw_obj_i.h,v 1.1 2007/09/28 16:57:20 mmundry Exp $;
 */

#ifndef PACKW_OBJ_I_H
#define PACKW_OBJ_I_H

struct coord;
struct sec_s_data;
struct seqprof;

void *coord_pack(struct coord *coord);
struct coord *coord_unpack(void *pak_coord);

void *sec_s_pack(struct sec_s_data *msdat);
struct sec_s_data *sec_s_unpack(void *pb);

void *seqprof_pack(struct seqprof *prof);
struct seqprof *seqprof_unpack(void *pp);

#endif				/* PACKW_OBJ_I_H */
