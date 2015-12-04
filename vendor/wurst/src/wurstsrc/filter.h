
/*
 *  filter.h
 *  
 *
 *  @author Laura Weidmann
 *  @date 06/03/2015 02:17:15 PM CEST
 *  
 */

float get_float_value(float* values, int index);
float* get_focf(char* struct_list, int size, struct prob_vec* query_vec, char* bin_fpc);
float* get_threshs(float* focf, int num_focf, int first_perc, int num);
int write_fpc_bin(char* vec_lib, char* target);


