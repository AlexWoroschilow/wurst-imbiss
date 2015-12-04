#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "fio.h"
#include "mprintf.h"
#include "matrix.h"
#include "coord.h"
#include "scor_set.h"
#include "scor_set_i.h"
#include "scoranlys_i.h"
#include "mcigraph_io_i.h"


/* --------- write_mci_mat -----------------------------------
 * takes a square matrix of floats and
 * write out a mci-style matrix file to the given filename 
 */

int write_mci_mat(const char *fname, size_t sz, float **cm) {
    FILE *tfile;
    const char *tfn = "write_mci_mat";
    size_t i,j;

    if ((tfile = mfopen(fname, "w", tfn))==NULL)
        return EXIT_FAILURE;

    mfprintf(tfile, "(mclheader\nmcltype matrix\ndimensions %ix%i\n)\n", 
             sz,sz);
    mfprintf(tfile, "(mclmatrix\nbegin\n");

    for (i=0; i<sz; i++) {
        mfprintf(tfile,"%i   ", i);
        for (j=0; j<sz;j++)
            mfprintf(tfile, "  %i:%f", j,cm[i][j]);
        mfprintf(tfile, " $\n");
    }

    mfprintf(tfile,")\n");
    if (fclose(tfile)) {
        mperror(tfn);
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}

