/*
 * 29 March 2005
 * Gundolf Schenk
 * $Id: altscores.c,v 1.2 2011/11/04 16:06:33 torda Exp $
 */

#ifndef _XOPEN_SOURCE
#    define _XOPEN_SOURCE 500
#endif

#include <stdio.h>  /* For debugging, maybe not necessary at end */
#include <stdlib.h>

#include "altscores.h"
#include "e_malloc.h"
#include "mprintf.h"  /* For debugging, maybe not necessary at end */
#include "pair_set.h"
#include "score_mat.h"

#define USE_THREADS

#ifdef USE_THREADS
/*
 * There is a buffer for random numbers. The rand function just hands back
 * the number at the top. If it notices that there are less than R_MIN
 * numbers left, he starts up the refill routine in the background by
 * posting a semaphore.
 * The refill routine can start, each time calculating N_FILL
 * numbers. He then keeps trying to put numbers on to the buffer.
 * When there is no room left, he goes into a cond_wait.
 * "refilling" is locked by "refill_lock". It is set on when refilling
 * takes place.  The signal to refill will not be sent if it is already
 * refilling.
 * Taking a package of random numbers out of the buffer and putting them
 * into the consumer buffer also means changing the refill buffer, so
 * this also takes "refill_lock" when necessary.
 */
#include <pthread.h>
#include <sched.h>
#include <semaphore.h>
#include <string.h>

enum {R_BUF_SIZ = 32768 };
enum {R_MIN = 20000 };

struct t_arg {
    long int r_buf [R_BUF_SIZ]; /* Where we store random numbers */
    unsigned *pn_rbuf;
    pthread_rwlock_t *r_buf_lock;
    sem_t *refill_sem;
};

/* ---------------- r_refill   --------------------------------------------
 * Get random numbers and put in a buffer.
 */
static void*
r_refill ( void *pa)
{
    enum { N_FILL = 1024 };
    long int fill_buf [N_FILL];
    const long int *plast = fill_buf + N_FILL;
    struct t_arg *a = pa;
    unsigned *pn_rbuf = a->pn_rbuf;
    long int *r_buf = a->r_buf;
    pthread_rwlock_t *r_buf_lock   = a->r_buf_lock;
    sem_t  *refill_sem             = a->refill_sem;

    while (1) {
        while ( (R_BUF_SIZ - *pn_rbuf) >= N_FILL) {
            long int *p;
            for (p = fill_buf; p < plast; p++)
                *p = lrand48();
            pthread_rwlock_wrlock (r_buf_lock);
            memcpy (r_buf + *pn_rbuf, fill_buf, N_FILL * sizeof (r_buf[0]));
            *pn_rbuf += N_FILL;
            pthread_rwlock_unlock (r_buf_lock);
        }
        sem_wait (refill_sem);
    }
    return NULL;
}

/* ---------------- get_n_rbuf --------------------------------------------
 */
static unsigned
get_n_rbuf ( unsigned *pn_rbuf, pthread_rwlock_t *r_buf_lock )
{
    unsigned r;
    pthread_rwlock_rdlock(r_buf_lock);
    r = *pn_rbuf;
    pthread_rwlock_unlock(r_buf_lock);
    return r;
}

/* ---------------- t_lrand48  --------------------------------------------
 * Threaded lrand48().
 * First version used to lock the buffer get a number, unlock and return.
 * It turns out, it was spending 20% of its time locking and unlocking.
 * Now, we have a local consum_buf and we refill this when necessary.
 */
static long int
t_lrand48 ( void )
{
    enum { N_CNSM = 8192 };
    static long int cnsm_buf [N_CNSM];
    static pthread_t tid;
    static unsigned n_rbuf;
    static struct t_arg t_arg;
    static pthread_rwlock_t r_buf_lock;
    static sem_t refill_sem;
    static unsigned short n_cnsm_buf = 0;
    static const char *this_sub = "t_lrand48";
    static char first_call = 1;

    if (first_call) {
        srand48 (55919929133);
        first_call = 0;
        t_arg.pn_rbuf     = &n_rbuf;
        t_arg.r_buf_lock  = &r_buf_lock;
        t_arg.refill_sem  = &refill_sem;
        pthread_rwlock_init (&r_buf_lock, NULL);
        sem_init (&refill_sem, 0, 0);
        if (pthread_create (&tid, NULL, r_refill, &t_arg)) {
            mperror (this_sub);
            err_printf (this_sub, "create random number thread failed\n");
            exit (EXIT_FAILURE);
        }
        while (get_n_rbuf(&n_rbuf, &r_buf_lock) < N_CNSM)
            sched_yield();
    }

    if ( n_cnsm_buf )                 /* Check easiest and most common case */
        return (cnsm_buf[--n_cnsm_buf]);  /* This order may avoid all locks */

    /*
     * The local consumers buffer is empty, but we can probably get
     * another load of numbers.
     */
    if (n_rbuf >= N_CNSM) {
        long int *src;
        pthread_rwlock_wrlock (&r_buf_lock);
        src = t_arg.r_buf + n_rbuf - N_CNSM;
        memcpy (cnsm_buf, src, N_CNSM * sizeof (t_arg.r_buf[0]));
        n_rbuf -= N_CNSM;
        pthread_rwlock_unlock (&r_buf_lock);
        n_cnsm_buf = N_CNSM;
    }

    if (n_rbuf < R_MIN) {                    /* Check if refill thread */
        int tmp;                             /* should be started */
        sem_getvalue (&refill_sem, &tmp);
        if (tmp == 0)
            sem_post (&refill_sem);
    }

    if (n_cnsm_buf)
        return (cnsm_buf[--n_cnsm_buf]);

    /*
     * Now we are in a bad place.
     * We did not even manage to get one package, so we have to
     * wait until it refills. In my few tests, this
     * part of the code was never reached - the random number
     * generator was always fast enough.
     */

#   ifdef want_verbose_checks
    {
        unsigned starve = 0;
        while (get_n_rbuf(&n_rbuf, &r_buf_lock) < N_CNSM) {
            sched_yield();
            starve++;
        }
        if (starve)
            err_printf ("l_rand48", "starved %u\n", starve);
    }
#   else              /* We do not want warnings about thread starvation */
        while (get_n_rbuf(&n_rbuf, &r_buf_lock) < N_CNSM)
            sched_yield();
#   endif             /* want_verbose_checks */
    {           /* Because we exited the loop above, we know the big buffer */
        long int *src;                          /* has enough numbers. Just */
        pthread_rwlock_wrlock (&r_buf_lock);    /* fill up our local buffer */
        src = t_arg.r_buf + n_rbuf - N_CNSM;
        memcpy (cnsm_buf, src, N_CNSM * sizeof (t_arg.r_buf[0]));
        n_rbuf -= N_CNSM;
        pthread_rwlock_unlock (&r_buf_lock);
    }
    n_cnsm_buf = N_CNSM;
    return (cnsm_buf[--n_cnsm_buf]);
}

#else /* do not use threads to generate random numbers */

#   define t_lrand48 lrand48

#endif /* USE_THREADS */

/* ---------------- pertubate_random_row_index_vec ------------------------
 * pertubates a random vector with numbers from 1 to (n_rows-1)
 * and indices from 1 to (n_rows-2)
 */
static size_t *
perturbate_vec (size_t * vec, size_t n)
{
    size_t i;
    size_t *p;
    size_t *last = vec + n;
    size_t *ndx = E_CALLOC ((n), sizeof (ndx[0]));
    size_t *ndx_beg = ndx;

    for (i = 0 ; i < n; i++)
        ndx[i] = t_lrand48() % n;

    for (p = vec; p < last; p++, ndx++) {
        size_t *dst = vec + *ndx;
        size_t tmp = *dst;
        *dst = *p;
        *p = tmp;
    }

    free (ndx_beg);

    return vec;
}

/* ---------------- create_random_row_index_vec ---------------------------
 * creates a random vector
 */
static size_t *
create_random_row_index_vec (size_t ** random_row_index_vec,
                             size_t * vec_length,
                             const struct pair_set * pair_set)
{
    size_t i, k;
    *random_row_index_vec = E_CALLOC ((pair_set->n), sizeof (size_t));

    for (i = 0, k = 0; i < pair_set->n; i++) {
        if (pair_set->indices[i][0] !=GAP_INDEX
            && pair_set->indices[i][1] !=GAP_INDEX) {
            (*random_row_index_vec)[k]=pair_set->indices[i][0];
            k++;
        }
    }
    *vec_length=k;

    *random_row_index_vec =
        perturbate_vec (*random_row_index_vec, *vec_length);

    return *random_row_index_vec;
}

/* ---------------- find_alt_path_score -----------------------------
 * Computes a suboptimal score given by the path ordered by
 * random_row_index_vec from left to right in score_mat of
 * size n_rows x n_cols ignoring borders.
 */
float
find_alt_path_score (const struct score_mat * score_mat,
                     const size_t * random_row_index_vec,
                     const size_t vec_length,
                     const struct pair_set * pair_set)
{
    float score = 0;
    size_t j, k;
    const size_t n = pair_set->n;
    for (j = 0, k = 0; j < n && k < vec_length; j++) {
        if (pair_set->indices[j][0] != GAP_INDEX
            && pair_set->indices[j][1] != GAP_INDEX) {
            score += score_mat->mat
                [random_row_index_vec[k]][pair_set->indices[j][1]];
            k++;
        }
    }

    return score;
}

/* ---------------- find_alt_path_score_simple ---------------------------
 * Computes a suboptimal score given by a random path
 * from left to right in score_mat
 * ignoring borders.
 *
 * This is with Gaps!
 */
float
find_alt_path_score_simple (const struct score_mat * score_mat,
                            const struct pair_set * pair_set)
{
    size_t * random_row_index_vec = NULL;
    size_t vec_length=0;
    float score;

    create_random_row_index_vec (&random_row_index_vec, &vec_length, pair_set);
    score = find_alt_path_score (score_mat, random_row_index_vec, vec_length,
                                 pair_set);

    free(random_row_index_vec);
    return score;
}
