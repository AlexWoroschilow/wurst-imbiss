/*
 * pair_set_sel.h
 *
 *  Created on: Sep 8, 2011
 *      Author: ibondarenko
 *
 * This defines the structures for passing alignments selected pairs or pair_set_sel around.
 * The functions are declared in pair_set_sel_i.h
 */

#ifndef PAIR_SET_SEL_H_
#define PAIR_SET_SEL_H_

/*
 *
 **/
char*
pair_set_sel_geti(struct pair_set *align);

/*
*
**/
char*
pair_set_sel_print(char* pss, size_t bs_len);

/*
 * converts all overposited letters to its uppercase equivalent
 **/
void
pair_set_print_prepare(struct seq *s1, struct seq *s2,
                       struct pair_set *align, char* pss);

/*
 * returns a pair_set, that contains only the subset of residues
**/
struct pair_set *
selected_pair_set_get(struct pair_set *source, char* pss);

/*
 * takes an alignment and two structures and
 * deletes the flag of a pair with max distance from the source.
 **/
char*
pair_set_sel_delmaxdistance(struct coord *c1, struct coord *c2,
                            struct pair_set *align, char* source);

size_t 
pair_set_sel_len(struct pair_set *align);


void 
printchar (int t);


#endif /* PAIR_SET_SEL_H_ */
