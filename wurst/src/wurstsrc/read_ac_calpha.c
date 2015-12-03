/*
 * 5 sep 2011
 * Read the influence file from an autoclass run for a calpha structure
 * classification and put the numbers into a table.
 *
 * This assumes our input is the text results from an autoclass
 * calpha-classification and it assumes a corresponding layout of the
 * data. This means we look for characteristic strings and
 * patterns. We do this with a mixture of the posix regex
 * functions and calls to strstr(). In principle, we do not need
 * both, but string matching is so much simpler, it is used to
 * quickly hop over large parts of the input file.
 *
 */

#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <regex.h>
#include <stdlib.h>
#include <string.h>

#include "bad_angle.h"
#include "coord.h"
#include "coord_i.h"
#include "e_malloc.h"
#include "fio.h"
#include "mprintf.h"
#include "gsldir/gsl_sf_erf.h"
#include "gsldir/gsl_linalg.h"
#include "gsldir/gsl_permutation.h"
#include "gsldir/gsl_matrix_double.h"
#include "gsldir/gsl_blas.h"
#include "prob_vec.h"
#include "prob_vec_i.h"
#include "read_ac_calpha_i.h"
#include "read_ac_only.h"
#include "yesno.h"

/* ---------------- Structures -------------------------------- */
/* number of class-properties for each descriptor */
/* and the class properties itself */
struct calpha_clssfcn {
    size_t n_class;                /* number of classes */
    size_t n_att;                  /* number of attributes per class */
    size_t cor_n_att;              /* number of corr attribs per class */
    float ***properties;           /* attribute properties */
    float *class_wt;               /* Normalised weight of each class */
    float *rcs;                    /* (heuristic) Relative class strength
                                      for each class */
    float **cov_matrix;            /* Covariance matrix */
    float ca_dist_error;
    float tau_error;
};

/* ---------------- Constants  -------------------------------- */
enum {CPROPERTIES = 20};
enum {DISTRIBUTION = 0}; /* 0 is SNcn; 1 is MNcn */
enum {DESCRIPTOR = 1};
enum {DESCRIPTORNUM = 2};
enum {INFLUENCE = 3};
enum {MEAN = 4};
enum {STDEV = 5};
enum {MEAN_DIV_STDEV = 6};
enum {MEAN_ = 7};
enum {STDEV_ = 8};
enum {DEFINED = 9};
enum {CORR_PARTNER = 10}; /* AttNum of correlation partner */
enum {CORR_VAL = 11};   /* The correlation value */
enum {CORR_PARTNER1 = 12};
enum {CORR_VAL1 = 13};
enum {CORR_PARTNER2 = 14};
enum {CORR_VAL2 = 15};
enum {CORR_PARTNER3 = 16};
enum {CORR_VAL3 = 17};
enum {CORR_PARTNER4 = 18};
enum {CORR_VAL4 = 19};



/* ---------------- calpha_get_n_att ---------------------------------
 * Keep reading the file and find the number of attributes in the
 * classification.
 */
static size_t
calpha_get_n_att (FILE *fp, char *buf, const int bufsiz)
{
    size_t n_att = 1;
    const char *heading = "num                        description ";
    const char *log = " Log Calpha";
    const char *magic_line_angles = "[0-9]+  [sinco()tau0-9]+ +[0-9].[0-9]+";
    const char *magic_line_log =    "[0-9]+  [Log]+ [Calpha()0-9distance+-]+ +"
                                    "[0-9].[0-9]+";
    if ( ! find_line (buf, bufsiz, fp, heading)) {
        return 0;
    }
    if ( ! find_regex (buf, bufsiz, fp, magic_line_angles, 0)) {
        return 0;
    }
    while ( find_regex (buf, bufsiz, fp, magic_line_angles, 0)) {
        n_att++;
    }
    rewind(fp);
    if ( find_line (buf, bufsiz, fp, log)) {
        n_att++;
        while ( find_regex (buf, bufsiz, fp, magic_line_log, 0))
            n_att++;
    }
    rewind(fp);
    return n_att;
}

/* ---------------- read_classes_calpha -----------------------
 * We have the overview of the classification, now read up the
 * classes and the values for the descriptors Ca-(i)-Ca-(i+2)-
 * distance and sine/cosine Tau
 */
static int
read_classes_calpha (FILE *fp, char *buf, const int bufsiz,
            struct calpha_clssfcn *ca_clssfcn, int corr_num)
{
    int i_class = 0;
    int i_atts = 0;
    float wt;
    int i,i_,j,r,l,l_,l__,k,k_,k__;
    float corval;
    int s = 0;
    int cur_ids[4];
    int this_num = 0;
    char distribution[10];
    char descriptor1[20];
    char descriptor2[30];
    float f_dummy;
    size_t ret = EXIT_SUCCESS; /* On error, set this and go to escape */
    regex_t get_prop, get_att_num, get_clss_wt, get_corr_matrix;
    const char *this_sub = "read_classes_calpha";
    const char *heading1  = " numb  t mtt  description           I-jk";
    const char *s_class_wt    = "normalized weight ";
    const char *s_get_att_num = "[0-9]+ [DR] [MNcnS]+ [sincoLog]+";
    const char *s_get_prop    = " +[0-9]+.[0-9]+ [(][? ?-][0-9]+"
                                "[.][0-9]+e[+-][0-9]+";
    const char *s_corr_matrix    = " Correlation matrix";
    const char *broke_ijk     = "Broke looking for I-jk value in \"%s\"\n";
    if (corr_num != 2 && corr_num != 4) {
        err_printf(this_sub,"At the moment, only 2 or 4 correlated"
                   " values are implemented for sine/cosine! You wanted: %d"
                   "Exiting..\n",corr_num);
        return EXIT_FAILURE;
    }
    if ( ! find_line (buf, bufsiz, fp, s_class_wt)) {
        return EXIT_FAILURE;
    }
    if ((wt = get_class_wt (buf)) < 0) {
        err_printf (this_sub, "Failed finding class weight on %s\n", buf);
        return EXIT_FAILURE;
    }
    if ( ! find_line (buf, bufsiz, fp, heading1)) {
        err_printf(this_sub,"Could not find heading, is this really an "
                            "AutoClass file?\n");
        return EXIT_FAILURE;
    }
    if (m_regcomp (&get_att_num, s_get_att_num) == EXIT_FAILURE) {
        err_printf(this_sub,"Could not get the Attributes, is this really an "
                            "AutoClass file?\n");
        return EXIT_FAILURE;
    }
    if (m_regcomp (&get_prop, s_get_prop) == EXIT_FAILURE) {
        err_printf(this_sub,"Could not get the Properties, is this really an "
                            "AutoClass file?\n");
        return EXIT_FAILURE;
    }
    if (m_regcomp (&get_clss_wt, s_class_wt) == EXIT_FAILURE) {
        err_printf(this_sub,"RegComp of Class weight wrong!\n");
        return EXIT_FAILURE;
    }
    if (m_regcomp (&get_corr_matrix, s_corr_matrix) == EXIT_FAILURE) {
        err_printf(this_sub,"RegComp of Correlation Matrix wrong!\n");
        return EXIT_FAILURE;
    }
    fseek(fp,0,SEEK_SET);

    while ( fgets (buf, bufsiz, fp) != NULL) {
        char *p = buf;
        regmatch_t pmatch[1];
        const size_t nmatch = 1;
        const int eflags = 0;
        r = regexec (&get_clss_wt, p, nmatch, pmatch, eflags);
        if (r != 0) {
            r = regexec (&get_att_num, p, nmatch, pmatch, eflags);
            if (r != 0) {
                r = regexec (&get_prop, p, nmatch, pmatch, eflags);
                if (r != 0) {
                    r = regexec (&get_corr_matrix, p, nmatch, pmatch, eflags);
                    if (r != 0) {
                        continue;
                    }
                    /* I got the lines with the Correlation matrices */
                    /* which looks like */
                    /* DATA_CORR_MATRIX
                       Correlation matrix (row & column indices are
                                           attribute numbers) */
                    /*         24     25           */
                    /* 24  1.000 -0.329            */
                    /* 25 -0.329  1.000            */

                    if (corr_num == 2) {
                      /* hard coded for sine/cosine correlations */
                      if (fscanf(fp, "%d", &this_num) != 1) {
                          err_printf(this_sub,"fscanf failed in reading"
                                     "correlations 1\n");
                          ret = EXIT_FAILURE;
                      }
                      if (fscanf(fp, "%*d %*d %*f %f %*d %*f %*f",
                          &ca_clssfcn->properties[i_class-1]
                                                 [this_num][CORR_VAL]) != 1) {
                          err_printf(this_sub,"fscanf failed in reading"
                                     "correlations 2\n");
                          ret = EXIT_FAILURE;
                      }
                      ca_clssfcn->properties[i_class-1]
                                [this_num][CORR_PARTNER] = (int) (this_num + 1);
                      ca_clssfcn->properties[i_class-1]
                                [this_num+1][CORR_PARTNER] = (int) (this_num);
                      ca_clssfcn->properties[i_class-1]
                                [this_num+1] [CORR_VAL] =
                                    ca_clssfcn->properties[i_class-1]
                                    [this_num][CORR_VAL];
                      /* count the correlated values for one class */
                      if (i_class == 1) {
                          ca_clssfcn->cor_n_att += 2;
                      }

                    /*        21     22     24     25 */
                    /* 21  1.000 -0.197 -0.024  0.033 */
                    /* 22 -0.197  1.000  0.033 -0.039 */
                    /* 24 -0.024  0.033  1.000 -0.226 */
                    /* 25  0.033 -0.039 -0.226  1.000 */
                    } else {
                      /* hard coded for neighbour-tau correlations */
                      if (fscanf(fp,"%d %*u %*u %*u", &this_num) != 1) {
                          err_printf(this_sub,"fscanf failed in reading"
                                     "correlations 1\n");
                          ret = EXIT_FAILURE;
                      }
                      for (j = 0 ; j < 5 ; j++) {
                        if (fscanf(fp, "%*u %f %f %f %f",
                          &ca_clssfcn->properties[i_class-1]
                                               [this_num+j][CORR_VAL1],
                          &ca_clssfcn->properties[i_class-1]
                                               [this_num+j][CORR_VAL2],
                          &ca_clssfcn->properties[i_class-1]
                                               [this_num+j][CORR_VAL3],
                          &ca_clssfcn->properties[i_class-1]
                                               [this_num+j][CORR_VAL4]) != 4) {
                            err_printf(this_sub,"fscanf failed in reading"
                                       "correlations 2\n");
                            ret = EXIT_FAILURE;
                        }
                        ca_clssfcn->properties[i_class-1][this_num+j]
                                         [CORR_PARTNER1] = (int) (this_num);
                        ca_clssfcn->properties[i_class-1][this_num+j]
                                         [CORR_PARTNER2] = (int) (this_num + 1);
                        ca_clssfcn->properties[i_class-1][this_num+j]
                                         [CORR_PARTNER3] = (int) (this_num + 3);
                        ca_clssfcn->properties[i_class-1][this_num+j]
                                         [CORR_PARTNER4] = (int) (this_num + 4);
                        if (j == 1) {
                          j++;
                        }
                      }
                      /* count the correlated values for a class
                         to use them later */
                      if (i_class == 1) {
                         ca_clssfcn->cor_n_att += 4;
                      }
                    }
                } else {
                    i_atts += 1;
                    ca_clssfcn->properties [i_class] [this_num] [DEFINED] = 1.0;
                    /* I got the line with the numbers */
                    /* which looks like "1.037 ( 6.25e-01  2.00e-01)  3.21e+00"
                                        "  (-1.70e-02  7.04e-01)" */
                    /* save all values in properties I may need them later */
                    p += pmatch->rm_so;
                    sscanf(buf, "%f ( %e  %e )  %e  ( %e  %e )",
                      &ca_clssfcn->properties[i_class][this_num][INFLUENCE],
                      &ca_clssfcn->properties[i_class][this_num][MEAN],
                      &ca_clssfcn->properties[i_class][this_num][STDEV],
                      &ca_clssfcn->properties[i_class][this_num]
                                             [MEAN_DIV_STDEV],
                      &ca_clssfcn->properties[i_class][this_num][MEAN_],
                      &ca_clssfcn->properties[i_class][this_num][STDEV_]);
                }
            } else {
                i_atts += 1;
                /* I got the line with the descriptors */
                /* which looks like */
                /* "002 004 R MNcn cos(tau(2))" or */
                /* "000 005 R SNcn Log Calpha-Calpha(+2)-distance(0)" */
                sscanf(buf, "%*u %d %*s %s %s %s",
                  &this_num,distribution,descriptor1,descriptor2);
                /* this attribute is defined */
                ca_clssfcn->properties[i_class][this_num][DEFINED] = 1.0;
                /* save the descriptor */
                if (descriptor1[0] == 'L') {
                    /* got Log of Calphas */
                    sscanf(descriptor2,"%*27c%f",&f_dummy);
                    ca_clssfcn->properties[i_class][this_num][DESCRIPTOR] = 2.0;
                    ca_clssfcn->properties[i_class][this_num][DESCRIPTORNUM] =
                        f_dummy;
                } else if (descriptor1[0] == 's') {
                    /* got sinus */
                    sscanf(descriptor1,"%*8c%f",&f_dummy);
                    ca_clssfcn->properties[i_class][this_num][DESCRIPTOR] = 0.0;
                    ca_clssfcn->properties[i_class][this_num][DESCRIPTORNUM] =
                        f_dummy;
                } else if (descriptor1[0] == 'c') {
                    /* got cosinus */
                    sscanf(descriptor1,"%*8c%f",&f_dummy);
                    ca_clssfcn->properties[i_class][this_num][DESCRIPTOR] = 1.0;
                    ca_clssfcn->properties[i_class][this_num][DESCRIPTORNUM] =
                        f_dummy;
                } else {
                    err_printf(this_sub,"Sorry, I cant read the descriptor, I"
                           "only know s, c, L and I got %d\n",descriptor1[0]);
                    ret = EXIT_FAILURE;
                }
                /* save the distribution */
                if (strcmp(distribution,"SNcn") == 0) {
                    ca_clssfcn->properties[i_class][this_num][DISTRIBUTION]=0.0;
                } else if (strcmp(distribution,"MNcn") == 0) {
                    ca_clssfcn->properties[i_class][this_num][DISTRIBUTION]=1.0;
                } else {
                    err_printf(this_sub,"Sorry, I cant read the distribution, I"
                                        "only know MNcn, SNcn\n");
                    ret = EXIT_FAILURE;
                }
            }
            if (r != 0) {
                err_printf (this_sub, broke_ijk, buf);
                ret = EXIT_FAILURE;
                goto escape;
            }
            if (i_atts == 2*ca_clssfcn->n_att) {
                i_atts = 0;
                i_class += 1;
            }
        } else {
            /* I got the line with the class weight, save the normalized one */
            /* CLASS  0 - weight 28534   normalized weight  1.621e-04 */
            /* .. relative strength  7.83e-01 ******* */
            sscanf(buf, "%*s %*u - weight %*u   %*s %*s %e %*s %*s %e",
                   &ca_clssfcn->class_wt[i_class], &ca_clssfcn->rcs[i_class]);
        }
    }
    /* ALLOC 2 dimensional array cov matrix */
    /* err_printf(this_sub,"CORNATT is %lu\n",ca_clssfcn->cor_n_att); */
    /* fflush(stderr); */
    ca_clssfcn->cov_matrix = (float **) E_CALLOC(ca_clssfcn->n_class
                                                 ,sizeof(float*));
    for (k = 0 ; k < ca_clssfcn->n_class ; k++) {
        ca_clssfcn->cov_matrix[k] = (float*) E_CALLOC((ca_clssfcn->cor_n_att *
                                                       ca_clssfcn->cor_n_att),
                                                       sizeof(float));
    }
    /* calculate the cov matrix for 2 correlated values (sin/cos) */
    if (corr_num == 2) {
      for (i = 0 ; i < ca_clssfcn->n_class ; i++) {
        l_ = 0;
        k_ = 0;
        for (k = 0 ; k < (ca_clssfcn->n_att + (ca_clssfcn->n_att
                          - ca_clssfcn->cor_n_att)) ; k++) {
            if (ca_clssfcn->properties[i][k][CORR_PARTNER] >= 1
                    && ca_clssfcn->properties[i][k][DEFINED] == 1) {
                ca_clssfcn->cov_matrix[i][k_+(l_*ca_clssfcn->cor_n_att)] =
                     ca_clssfcn->properties[i][k][STDEV]
                     *  ca_clssfcn->properties[i][k][STDEV];
                ca_clssfcn->cov_matrix[i][k_+1+(l_*ca_clssfcn->cor_n_att)] =
                    (ca_clssfcn->properties[i][k][CORR_VAL]
                     * ca_clssfcn->properties[i][k][STDEV]
                     * ca_clssfcn->properties[i][k+1][STDEV]);
                l_++;
                ca_clssfcn->cov_matrix[i][k_+1+(l_*ca_clssfcn->cor_n_att)] =
                        ca_clssfcn->properties[i][k+1][STDEV] * 2;
                ca_clssfcn->cov_matrix[i][k_+(l_*ca_clssfcn->cor_n_att)] =
                    (ca_clssfcn->properties[i][k+1][CORR_VAL]
                     * ca_clssfcn->properties[i][k][STDEV]
                     * ca_clssfcn->properties[i][k+1][STDEV]);
                l_++;
                k_+=2;
                k++;
            }
        }
      }
    } else { /* do it for 4 correlated values (sin1/cos1) & (sin2/cos2)*/
      for (i = 0 ; i < ca_clssfcn->n_class ; i++) {
        l = 0;
        l_ = 0;
        k__ = 0;
        l__ = 0;
        for (k = 0 ; k < (ca_clssfcn->n_att + (ca_clssfcn->n_att
                        - ca_clssfcn->cor_n_att)+1) ; k++) {
            if (l__ == 4 || l__ == 0) {
                l__ = 0;
                for (i_ = 0 ; i_ < 4 ; i_++) {
                    cur_ids[i_] = k + i_ + s;
                    if (i_ == 1) {
                        s++;
                    }
                }
                s = 0;
            }
            if (ca_clssfcn->properties[i][k][CORR_PARTNER1] >= 1
                  && ca_clssfcn->properties[i][k][DEFINED] == 1) {
                for (k_ = 0 ; k_ < 4 ; k_++) {
                    if (k_ == 0) {
                        corval = ca_clssfcn->properties[i][k][CORR_VAL1];
                    } else if (k_ == 1) {
                        corval = ca_clssfcn->properties[i][k][CORR_VAL2];
                    } else if (k_ == 2) {
                        corval = ca_clssfcn->properties[i][k][CORR_VAL3];
                    } else if (k_ == 3) {
                        corval = ca_clssfcn->properties[i][k][CORR_VAL4];
                    }
                    if (k_ == l__) {
                        ca_clssfcn->cov_matrix[i]
                                    [l+l_+(k__ * ca_clssfcn->cor_n_att)] =
                                      ca_clssfcn->properties[i][k][STDEV]
                                      * ca_clssfcn->properties[i][k][STDEV];
                    } else {
                        ca_clssfcn->cov_matrix[i]
                                    [l+l_+(k__ * ca_clssfcn->cor_n_att)] =
                                      ca_clssfcn->properties[i]
                                                 [cur_ids[k_]][STDEV]
                                      * ca_clssfcn->properties[i][k][STDEV]
                                      * corval;
                    }
                    l++;
                    /*l_++;*/
                }
                if (l__ == 3) {
                    l_ += 4;
                }
                l = 0;
                l__++;
                k__++;
            }
        }
      }
    }
 escape:
    regfree (&get_att_num);
    regfree (&get_prop);
    regfree (&get_clss_wt);
    regfree (&get_corr_matrix);
    return ret;
}

/* ---------------- new_calpha_clssfcn-----------------------------
 * Allocate and set up a new Calpha Classification
 */
static struct calpha_clssfcn *
new_calpha_clssfcn ( const size_t n_class, const size_t n_att, float t_error,
                     float c_error)
{
    size_t i,j;
    struct calpha_clssfcn *ca_clssfcn =E_MALLOC(sizeof (struct calpha_clssfcn));
    ca_clssfcn->n_class = n_class;
    ca_clssfcn->n_att   = n_att;
    ca_clssfcn->cor_n_att = 0; /* will be filled later dependent on corvals */
    /* set the absolute error vals */
    ca_clssfcn->tau_error = t_error;
    ca_clssfcn->ca_dist_error = c_error;
    /* Alloc 3D-array - the one in matrix.c does strange things
       and initialise the correlation and defined values */
    /*(float ***) ca_clssfcn->properties;*/
    ca_clssfcn->properties = (float ***) malloc (sizeof (float **) * n_class);
    for (i = 0 ; i < n_class ; i++) {
        ca_clssfcn->properties[i] =(float **)malloc(sizeof (float *)*(n_att*2));
        for (j = 0 ; j < (n_att*2) ; j++) {
            ca_clssfcn->properties[i][j] = (float *) malloc
                                           (sizeof (float) * CPROPERTIES);
            ca_clssfcn->properties[i][j][DEFINED] = 0.0;
            ca_clssfcn->properties[i][j][CORR_PARTNER] = 0.0;
            ca_clssfcn->properties[i][j][CORR_PARTNER1] = 0.0;
            ca_clssfcn->properties[i][j][CORR_PARTNER2] = 0.0;
            ca_clssfcn->properties[i][j][CORR_PARTNER3] = 0.0;
            ca_clssfcn->properties[i][j][CORR_PARTNER4] = 0.0;
        }
    }
    ca_clssfcn->class_wt = E_MALLOC (n_class * sizeof (float));
    ca_clssfcn->rcs = E_MALLOC (n_class * sizeof (float));
    return ca_clssfcn;
}

/* ---------------- calpha_clssfcn_return --------------------------------
 * Returns the values of the attributes from the descriptors of a
 * specified class
 */
float
calpha_clssfcn_return( struct calpha_clssfcn* ca_clssfcn,
                   int i_class, int pos, int att_num)
{
    return ca_clssfcn->properties[i_class][pos][att_num];
}

/* ---------------- calpha_clssfcn_destroy ------------------------
 * This clean up function will be called by the perl interface.
 */
void
calpha_clssfcn_destroy(struct calpha_clssfcn* ca_clssfcn)
{
    size_t i,j;
    free (ca_clssfcn->class_wt);
    free (ca_clssfcn->rcs);
    for (i = 0 ; i < ca_clssfcn->n_class ; i++) {
        free(ca_clssfcn->cov_matrix[i]);
    }
    free(ca_clssfcn->cov_matrix);
    for (i = 0 ; i < ca_clssfcn->n_class ; i++) {
        for (j = 0 ; j < ca_clssfcn->n_att ; j++) {
            free(ca_clssfcn->properties[i][j]);
        }
        free(ca_clssfcn->properties[i]);
    }
    free(ca_clssfcn->properties);
    /*
    kill_3d_array ((void *)ca_clssfcn->properties);
    */
    free (ca_clssfcn);
}

/* ---------------- att_size_calpha -------------------------------
 * Return the number of attributes
 */
size_t
att_size_calpha (const struct calpha_clssfcn *ca_clssfcn)
{
    return ca_clssfcn->n_att;
}

/* ---------------- ac_size_calpha  ----------------------------------
 * Return the fragment length size associated with a
 * classification. This is *not* the number of classes. It will
 * typically be a number between 4 and 12.
 */
static size_t
ac_size_calpha (const struct calpha_clssfcn *ca_clssfcn)
{
    const char *this_sub = "ac_size";
    if (ca_clssfcn->n_att < 4) {
        err_printf(this_sub, "there must be an error, "
                             "fragment_length cant be smaller than 4");
        return 0;
    } else if (ca_clssfcn->n_att == 4) {
        /* special case */
        return 4;
    } else if (ca_clssfcn->n_att > 4) {
        /* special case */
        return ((ca_clssfcn->n_att+8)/3);
    } else {
        err_printf(this_sub, "there must be an error, "
                             "fragment_length cant be smaller than 4");
        return 0;
    }
}

/* ---------------- ac_nclass_calpha --------------------------
 * Return the number of classes in a classification.
 */
size_t
ac_nclass_calpha ( const struct calpha_clssfcn *ca_clssfcn)
{
    return ca_clssfcn->n_class;
}

#ifdef want_class_stats_for_interface
/* ---------------- ac_class_wt --------------------------------
 * Return the weight of a class in the classification.
 */
float
ac_class_wt_calpha ( const struct calpha_clssfcn *ca_clssfcn, int i_class)
{
    return ca_clssfcn->class_wt[i_class];
}

/* -------------- ac_class_influ -------------------------------
 * Return the influence value for an attribute in a given class
 */
float
ac_class_influ_calpha ( const struct calpha_clssfcn *ca_clssfcn, 
                        int attrib, int i_class) 
{
    return ca_clssfcn->properties[i_class][attrib][INFLUENCE];
}
#endif /* want_class_stats_for_interface */


/* ---------------- ac_read_calpha  ----------------------------
 * Given a pointer to a filename, read in the data.
 * This is an interface function, visible to the outside world.
 */
struct calpha_clssfcn *
ac_read_calpha (const char *fname, float t_error, float c_error, int corr_num)
{
    FILE *fp;
    size_t n_class, n_att;
    struct calpha_clssfcn *ca_clssfcn = NULL;
    const char *this_sub = "ac_read_calpha";
#   ifndef BUFSIZ
        enum {BUFSIZ = 512};
#   endif
    char buf [BUFSIZ];

    if ((fp = mfopen (fname, "r", this_sub)) == NULL) {
        return (NULL);
    }
    if ((n_class = get_n_class (fp, buf, BUFSIZ)) == 0) {
        goto end;
    }
    if ((n_att = calpha_get_n_att (fp, buf, BUFSIZ)) == 0) {
        goto end;
    }
    ca_clssfcn = new_calpha_clssfcn (n_class, n_att, t_error, c_error);
    if (read_classes_calpha (fp, buf, BUFSIZ, ca_clssfcn, corr_num)
        ==EXIT_FAILURE){
          goto end;
    }
 end:
    fclose (fp);
    return (ca_clssfcn);
}

/* --------------- calpha_strct_ComputeMembership ---------------
 * calpha_strct_ComputeMembership
 * Core of the original function computeMembershipStrct
 * mship should be initialized, normally to class_weight or 1
 * Mixed model
 */
static int
calpha_strct_ComputeMembership (float *mship, const float *test_vec,
                                const struct calpha_clssfcn *ca)
{
    const char this_sub[] = "computeMembershipStrct";
    size_t i, j , k, l;
    double *C_matrix;
    gsl_matrix_view C;
    gsl_permutation *p;
    gsl_matrix* C_inv;
    double C_det;
    double *y_array;
    gsl_vector_view y;
    gsl_vector* C_inv_y;
    double exponent = 1.0;
    size_t to_mall;
    int s;
    i = 0;
    if (mship == NULL){
        err_printf (this_sub, "mship is null\n");
        return EXIT_FAILURE;
    }
    if (test_vec == NULL){
       err_printf (this_sub, "test_vec is NULL\n");
       return EXIT_FAILURE;
    }
    to_mall = (ca->cor_n_att*ca->cor_n_att);
    C_matrix = E_CALLOC(to_mall, sizeof (C_matrix[0]));
    y_array = E_CALLOC(ca->cor_n_att, sizeof(double));
    p = gsl_permutation_alloc(ca->cor_n_att);
    C_inv = gsl_matrix_alloc(ca->cor_n_att, ca->cor_n_att);
    C_inv_y = gsl_vector_alloc(ca->cor_n_att);
    /* CALC mship for each class */
    for (i = 0; i < ca->n_class; i++) {
        for (k = 0; k < (ca->cor_n_att*ca->cor_n_att); k++) {
            C_matrix[k] = ca->cov_matrix[i][k];
        }
        C = gsl_matrix_view_array(C_matrix, ca->cor_n_att, ca->cor_n_att);
        gsl_linalg_LU_decomp(&C.matrix, p, &s);
        gsl_linalg_LU_invert(&C.matrix, p, C_inv);
        C_det = gsl_linalg_LU_det(&C.matrix, s);
        /* PRE MATRIX CALCULATIONS DONE */
        l = 0;
        for (j = 0 ; j < (ca->n_att + (ca->n_att - ca->cor_n_att)) ; j++) {
            if ((int)ca->properties[i][j][DEFINED]) {
                if ( ca->properties[i][j][DISTRIBUTION] == 0.0) {
                    /* SINGLE_NORMAL */
                    mship[i] *= (1/ca->properties[i][j][STDEV])*
                        (gsl_sf_erf_Q((test_vec[j] - fabs(ca->ca_dist_error)
                                  - ca->properties[i][j][MEAN])
                                  / ca->properties[i][j][STDEV]) -
                                    gsl_sf_erf_Q((test_vec[j]
                                         + fabs(ca->ca_dist_error)
                                         - ca->properties[i][j][MEAN])
                                         / ca->properties[i][j][STDEV]));
                } else if (ca->properties[i][j][DISTRIBUTION] == 1.0) {
                    /* MULTI_NORMAL */
                    if (test_vec[j] == BAD_ANGLE_FLOAT ) {
                        y_array[l] = 0.0;
                        err_printf(this_sub,"BAD_ANGLE detected");
                    } else {
                        y_array[l] = test_vec[j] - ca->properties[i][j][MEAN];
                    }
                    l++;
                } else {
                    err_printf(this_sub, "Unknown classmodel : %f\n",
                           ca->properties[i][j][DISTRIBUTION]);
                    return EXIT_FAILURE;
                }
            }
        }
        /*  CALC MATRIX  */
        y = gsl_vector_view_array(y_array,ca->cor_n_att);
        gsl_blas_dgemv(CblasNoTrans, 1.0, C_inv, &y.vector, 0.0, C_inv_y);
        gsl_blas_ddot(C_inv_y, &y.vector, &exponent);
        exponent *= -0.5;
        mship[i] *=  pow(2 * ca->tau_error, ca->cor_n_att) * exp(exponent) /
                                                             sqrt(fabs(C_det));
        /* Look for some errors */

        if (isnan(mship[i]) != 0) {
            err_printf(this_sub,"found nan in %lu\n",i);
            return EXIT_FAILURE;
        } else if ( exponent == 0) {
            err_printf(this_sub,"exponent 0 in mship[%lu]\n",i);
        }
    }
    /*  free MATRIX  */
    gsl_vector_free(C_inv_y);
    gsl_matrix_free(C_inv);
    free_if_not_null(y_array);
    free_if_not_null(C_matrix);
    gsl_permutation_free(p);
    return EXIT_SUCCESS;
}

/* ---------------- calpha_seq_strct_2_prob_vec -----------------
 * Calculates probabilities using combinations of sequence and
 * structure, sequence profile and structure, sequence only,
 * sequence profile only, and structure only.
 */

static struct prob_vec *
calpha_seq_strct_2_prob_vec (struct coord *structure,
                   const size_t size, const struct calpha_clssfcn *ca,
                   const enum yes_no norm)
{
    const char *this_sub = "calpha_seq_strct_2_prob_vec";
    size_t i, j, compnd_len;
    struct prob_vec *pvec;
    int frag_len = ac_size_calpha(ca);
    const size_t n_pvec = size - frag_len + 1;
    if ((pvec = new_pvec (frag_len, size, n_pvec,
                          ca->n_class)) != NULL) {
        /* was ca->n_att, size, n_pvec .. i put fraglen */
        /* i also had size - frag_len for n_pvec above */
        /* Initialize membership */
        for (i = 0 ; i < n_pvec ; i++) {
            for (j = 0 ; j < ca->n_class ; j++) {
                pvec->mship[i][j] = ca->class_wt[j];
            }
        }
        if (structure) { /* Calculate structure membership */
            float *fragment = NULL;
            coord_calc_tau(structure);
            fragment = E_CALLOC((ca->n_att + (ca->n_att - ca->cor_n_att) + 1),
                                sizeof(float));
            for (i = 0 ; i < n_pvec ; i++) { /* for every fragment   */
                calpha_getFragment(i, ca->n_att, ca->cor_n_att, structure,
                                   fragment);
                if (calpha_strct_ComputeMembership(pvec->mship[i], fragment,ca)
                                         == EXIT_FAILURE) {
                    prob_vec_destroy (pvec);
                    err_printf(this_sub,"Mship NULL!\n");
                    return NULL;
                }
            }
            free(fragment);
        }
        if (norm == YES) { /* Normalise the score */
            for (i = 0; i < n_pvec; i++){
                double sum = 0.0;
                for (j = 0; j < ca->n_class; j++)
                    sum += pvec->mship[i][j];
                for (j = 0; j < ca->n_class; j++)
                    pvec->mship[i][j] /= sum;
            }
            pvec->norm_type = PVEC_TRUE_PROB;
        }
    }
    if ( structure ) {
        compnd_len = structure->compnd_len;
    } else {
        compnd_len = 0;
    }
    if(compnd_len > 0) {           /* finally read compound */
        pvec->compnd = E_MALLOC(compnd_len*sizeof(char));
        memmove(pvec->compnd, structure->compnd, compnd_len);
    }
    pvec->compnd_len = compnd_len;
    return pvec;
}

/* ---------------- calpha_strct_2_prob_vec  --------------------
 * Interface function to calculate probabilities, only using the
 * Calpha structure terms
 */
struct prob_vec *
calpha_strct_2_prob_vec (struct coord *structure,
                  const struct calpha_clssfcn *ca, const int norm)
{
    const char *this_sub = "calpha_strct_2_prob_vec";
    err_printf(this_sub, "entering calpha_strct_2_prob_vec...\n");

    if (!structure){
        err_printf (this_sub, "No Structure Input!\n");
        return NULL;
    }
    if (!ca) {
        err_printf (this_sub, "No Classification Input!\n");
        return NULL;
    }
    return calpha_seq_strct_2_prob_vec (structure, structure->size,
                                        ca, norm);
}


/* ---------------- calpha_getFragment --------------------------
 * contructs a descriptor (tau angles) and Ca-distances of the
 * fragment of length and saves it in an array so that the indices
 * are similar to the indices used in the AutoClass output
 *
 * remember to free the fragment after usage!
 *
 * return value:
 * the descriptor - if all angles were defined
 * NULL           - if a bad angle was encountered
 */

void
calpha_getFragment(const size_t residue_num, const size_t n_att,
                   const size_t cor_n_att, struct coord *structure,
                   float *frag)
{
    /* Attributes are mapped to the index of the classification
       which makes it much easier */
    size_t i = 0, j = 0, t = 0;
    float angle;
    float distance;
    /* get all Tau angles */
    for (i = residue_num + 2, j = 3 ; i < ((cor_n_att/2) + residue_num + 2)
         ; i++, j+= 2) {
        angle = coord_tau(structure, i, 0); 
        if (angle < 0) {
            frag[j] = BAD_ANGLE_FLOAT;
            frag[j+1] = BAD_ANGLE_FLOAT;
        } else {
            if (angle > M_PI)    /* We have 0->2*pi. Make it -pi->pi*/
                angle -= (2*M_PI);
            frag[j] = sin(angle);
            frag[j+1] = cos(angle);
            j++;
        }
    }
    /* get all Ca-distances */
    t = j;
    for (i = residue_num , j = t - 1 ; j < (t+(n_att - cor_n_att - 1))
         ; i++, j++) {
        distance = coord_ca_ca_dist(structure,i,i+2,1);
        if (distance < 3.8 || distance > 7.7 ) {
            frag[j] = BAD_DISTANCE_FLOAT;
        } else {
            frag[j] = log(distance);
        }
    }
}
