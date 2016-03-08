/* mcigraph_io_i.h
 * 18-06-2004
 * Routines for writing out score matrices as weighted graphs 
 * in the mci matrix format
 */
#ifndef MCIGRAPH_IO_I_H
#define MCIGRAPH_IO_I_H
int write_mci_mat(const char *fname, size_t sz, float **cm);
#endif
