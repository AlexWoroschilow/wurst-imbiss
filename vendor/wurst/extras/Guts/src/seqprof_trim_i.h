/* seqprof_trim_i.h
 * 8-7-2004
 * extract range of values from a sequence profile
 */
#ifndef SEQPROF_TRIM_I_H
#define SEQPROF_TRIM_I_H

struct seqprof *
seqprof_trim ( struct seqprof *prf, size_t start, size_t end);
struct seqprof *
seqprof_merge ( struct seqprof *prf1, struct seqprof *prf2);
struct seq *
seq_merge ( struct seq *seq1, struct seq *seq2);

#endif
