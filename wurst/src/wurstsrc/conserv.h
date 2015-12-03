/*
 * conserv.h
 *
 *  Created on: Jan 4, 2012
 *      Author: ibondarenko
 */

#ifndef CONSERV_H_
#define CONSERV_H_

struct prob_vecs {
    size_t maxnum;
    size_t curnum;
    struct prob_vec **pv;
};

struct sequences_list {
    size_t maxnum;
    size_t curnum;
    struct seq **ps;
};

struct prob_vecs *
initveclist(size_t num);

void
cleanveclist(struct prob_vecs *ps);

void
deleteveclist(struct prob_vecs *ps);

struct sequences_list *
initseqlist(size_t num);

void
cleanseqlist(struct sequences_list *ps);

void
deleteseqlist (struct sequences_list *ps);

void
addvec(struct prob_vecs* pvs, struct prob_vec *vec);

void
addseq(struct sequences_list* pvs, struct seq *seq);

float*
get_seq_conserv(struct pair_set *pair_set, struct sequences_list* ps);

float *
getconserv(struct prob_vec* pv);

float *
getconservvec(struct pair_set *pair_set, struct prob_vecs* pvs);

char *
printconserv(float * consv, const size_t len, struct pair_set *pair_set,
             const size_t startpos, const size_t maxlen,
             const size_t strNum1);

char *
conserv_str(float * consv, size_t len);

#endif /* CONSERV_H_ */
