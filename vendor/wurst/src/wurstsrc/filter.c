/*
 *  filter.c
 *  new measurement of protein similarity: focf
 *  
 *  @author Laura Weidmann
 *  @date 06/03/2015 02:17:15 PM CEST
 *  
 */

#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <dirent.h>
#include <libgen.h>

#include "e_malloc.h"
#include "fio.h"
#include "matrix.h"
#include "mgc_num.h"
#include "mprintf.h"
#include "prob_vec.h"
#include "prob_vec_i.h"
#include "scratch.h"
#include "seq.h"
#include "yesno.h"

float calc_focf_of(float* fpc1, int len1, float* fpc2, int len2, int classes);
void prob_vec_2_frag_per_class(struct prob_vec* pv, float* fpc);
float num_of_frags(int class, struct prob_vec* pv);
void generate_file_name(char* compare_file, char* id, char* dir);

/* ------------------ get_focf ----------------------------
 * returns list of focf for a structure list and a query structure
 * called by perl
 */

float * 
get_focf(char* struct_list, int size, struct prob_vec* query_vec, char* bin_fpc){
   float* focf_list = E_MALLOC(sizeof(float) * size);
   int classes = query_vec->n_class;
   float* query_fpc = E_MALLOC(sizeof(float) * classes);
   int i, j, compare_len;
   for (i=0; i<classes; i++){
      query_fpc[i] = num_of_frags(i, query_vec);
   }
   prob_vec_2_frag_per_class(query_vec, query_fpc);
   char* compare_id = strtok(struct_list, ",");
   FILE* fpc_file = fopen(bin_fpc, "r");
   char bin_id[8];
   float* compare_fpc = E_MALLOC(sizeof(float) * classes);
   int strcmp_value;
   for (i=0; i<size; i++){
      focf_list[i] = -1.0; // if not found 
      for (j=0;j<263207;j++){
         fread(bin_id, sizeof(char), 8, fpc_file);
         fread(&compare_len, sizeof(int), 1, fpc_file);
         fread(compare_fpc, sizeof(float), classes , fpc_file);
         strcmp_value = strcmp(compare_id, bin_id);
         if (strcmp_value > 0){ continue; } 
         else if (strcmp_value == 0){
           focf_list[i] = calc_focf_of(query_fpc, query_vec->prot_len, compare_fpc, compare_len, classes);
           break;
         } else { break; }
      }
      compare_id = strtok(NULL, ",");
   }
   free(query_fpc);
   free(compare_fpc);
   return focf_list;
}

/* ------------------ get_threshs----------------------------
 * returns list of focf-cutoffs
 * modifiable by first_perc and num of cutoffs
 */
float* get_threshs(float* focf, int num_focf, int first_perc, int num){
  int* bins = E_MALLOC(sizeof(int) * 1001);
  float* threshs = E_MALLOC(sizeof(float) * num);
  int* percs = E_MALLOC(sizeof(int) * num);
  percs[0] = first_perc*num_focf / 100;
  threshs[0] = num_focf;
  int rest_perc = 100-first_perc;
  int i;
  for (i=1; i<num; ++i)  {
    percs[i] = percs[i-1] + (int) (num_focf*rest_perc/(100*num));
    threshs[i] = 0;
  }
  for (i=0; i<1001; ++i) bins[i] = 0;
  for (i=0; i<num_focf; ++i){
     bins[(int)(1000*focf[i])] += 1; 
  }
  int sum = 0; 
  int current_perc = 0;
  for (i=1000; i>=0; --i){
     sum += bins[i];
     if (sum >= percs[current_perc]) {
        threshs[current_perc] = (float) i/ 1000.0;
        current_perc ++;
        if (current_perc == num) break;
     }
  }
  free_if_not_null(bins);
  free_if_not_null(percs);
  return threshs;
}

/* ------------------ get_float_value----------------------------
 * returns a float of a float list, needed by Perl interface
 */
float
get_float_value(float* values, int index){
   return values[index];
}

/* ------------------ prob_vec_2_frag_per_class ---------------------------
 * writes to frag_per_class sum of all prob/fragment for each class
 */
void
prob_vec_2_frag_per_class(struct prob_vec* pv, float* frag_per_class){
   int c; int i;
   if (pv->mship == NULL) prob_vec_expand(pv);
   for (c=0; c<pv->n_class; c++){
      frag_per_class[c] = 0.0;
      for (i=0; i<(pv->prot_len-6) ; i++){
         frag_per_class[c] += pv->mship[i][c];
      }
   }
}

/* ------------------ num_of_frags ---------------------------
 * return a float that contains number of fragments of a certain class
 */
float
num_of_frags(int class, struct prob_vec* pv){
   if (pv->mship == NULL) prob_vec_expand(pv);
   int i;
   float frags = 0.0;
   for (i=0; i<pv->prot_len-6; i++){
      frags += pv->mship[i][class];
   }
   return frags;
}

/* ----------------- calc_focf_of ----------------------------------
 * for two proteins calculate focf
 */
float
calc_focf_of(float* fpc1, int len1, float* fpc2, int len2, int classes){
   int i=0;
   int c=0;
   float focf = 0.0;
   for (c=0; c<classes; c++){
      focf += fminf(fpc1[c], fpc2[c]);      
   }
   // normalize by size of smallest protein
   if (len1 < len2) { focf = focf/(float) (len1);
   } else { focf = focf/(float) (len2); }

   return focf;
}


/* ------------------ generate_file_name ---------------------------
 */
void
generate_file_name(char* compare_file, char* id, char* dir){
  compare_file[0] = 0;
  strcat(compare_file, dir);
  compare_file[strlen(dir)] = '/';
  compare_file[strlen(dir)+1] = 0;
  strcat(compare_file, id);
  strcat(compare_file, ".vec");
  compare_file[strlen(dir)+10] = 0;
}


/* ------------------ write_fpc_bin ----------------------------------
 * for a directory with .vec files write a fpc file
 */
int
write_fpc_bin(char* vec_lib, char* target){
  struct dirent **entries;
  int n = scandir(vec_lib, &entries, NULL, alphasort);
  struct prob_vec** lib_vec = E_MALLOC(sizeof(struct prob_vec*) * n);
  int i, j;
  char id[16];
  for (i=0; i<16; ++i) {
      id[i] = 0;
  }
  int length = 0;
  char compare_file[strlen(vec_lib)+10];
  strncpy(id, entries[3]->d_name, 5); id[5] = 0; 
  generate_file_name(compare_file, id, vec_lib);
  lib_vec[0] = prob_vec_read(compare_file);
  int classes = lib_vec[0]->n_class;
  float* fpc = E_MALLOC(sizeof(float) * classes);
  for (i=0; i<classes; ++i) {
     fpc[i] = 0.0;
  }
  FILE* bin_file = fopen(target, "w");

  for (i=0; i<n; ++i) {
     if (!strstr(entries[i]->d_name, ".vec")) continue;
     // 1) ID from filename
     id[0] = 0; strncpy(id, entries[i]->d_name, 5); id[5] = 0;
     // 2) fpc -> read .vec to a prob_vec
     generate_file_name(compare_file, id, vec_lib);
     lib_vec[i] = prob_vec_read(compare_file);
     // 3) get length form prob_vec
     length = lib_vec[i]->prot_len;
     prob_vec_2_frag_per_class(lib_vec[i], fpc);
     fwrite(id, sizeof(char) , 8, bin_file);
     fwrite(&length, sizeof(int), 1, bin_file);
     fwrite(fpc, sizeof(float) , classes, bin_file);
     free(entries[i]);
     prob_vec_destroy(lib_vec[i]);
  }
  free(entries);
  free(lib_vec);
  id[0] = 0; strcpy(id, "00EOF"); id[5] = 0; 
  fwrite(id, sizeof(char), 8, bin_file);
  fclose(bin_file);  
  return n;
}








