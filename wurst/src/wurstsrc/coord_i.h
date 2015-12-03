/*
 * 23 October 2001
 * This defines the interface to the coordinate routines.
 * It does *not* define the internal structures.
 * rcsid = $Id: coord_i.h,v 1.3 2012/04/18 12:02:04 ibondarenko Exp $
 */

#ifndef COORD_I_H
#define COORD_I_H

struct coord;
int            coord_2_bin (struct coord *c, const char *fname);
char *         coord_name (struct coord *c);
char *         coord_get_numbering (struct coord *c);
char *         coord_name_thr (struct coord *c);
size_t         coord_size (const struct coord *c);
struct coord * coord_read (const char *fname);
struct seq   * coord_get_seq (const struct coord *c);
int            coord_copy_seq (const struct coord *src, struct coord *dst);
int            coord_move_seq (struct coord *src, struct coord *dst);
void           coord_reset_phi_psi (struct coord *c);
void           coord_calc_psi (struct coord *c);
void           coord_calc_phi (struct coord *c);
void           coord_calc_omega (struct coord *c);
void           coord_calc_theta (struct coord *c);
void           coord_calc_tau (struct coord *c);
float          coord_omega (struct coord *c, const size_t i, 
                            const float shift_min);
float          coord_psi (struct coord *c, const size_t i, 
                          const float shift_min);
float          coord_phi (struct coord *c, const size_t i,
                          const float shift_min);
float          coord_theta (struct coord *c, const size_t i,
                            const float shift_min);
float          coord_tau (struct coord *c, const size_t i,
                          const float shift_min);
void           coord_nm_2_a (struct coord *c);
void           coord_a_2_nm (struct coord *c);
struct coord * coord_template (const struct coord *c, size_t i );
struct coord * coord_copy (const struct coord *c);
struct coord * coord_trim (struct coord *c, const size_t size);
struct coord * coord_construct_frgmt (size_t n);
void           coord_add_ca (struct coord *frgmt, const size_t pos,
                             float x, float y, float z);
void           coord_destroy (struct coord *c);
int            coord_has_sec_s (const struct coord *c);
float          coord_c_n_dist (const struct coord *c,
                               const unsigned int i, const unsigned int j,
                               const unsigned int sqrt_flag);
float          coord_ca_ca_dist (const struct coord *c,
                               const unsigned int i, const unsigned int j,
                               const unsigned int sqrt_flag);
size_t         coord_neighbours_cb(const struct coord *c,
                                   const size_t i, const float radius);
int            coord_has_steric_clash(const struct coord *c);
float          coord_gyration(const struct coord *c);
int            coord_ok (const struct coord *c);
#endif  /* COORD_I_H */
