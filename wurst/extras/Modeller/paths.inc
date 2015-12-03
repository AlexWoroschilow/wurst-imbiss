# Here are some paths which can be used by the other scripts
# in this directory.

# There is a little trickiness. Unlike every other file, we do
# not 'use strict' here. We offer a pile of definitions. I am
# not sure if this constitutes name space pollution, when the
# callers are using 'use strict'.
# rcsid = $Id: paths.inc,v 1.1 2007/09/28 16:57:17 mmundry Exp $

use FindBin;
*MATRIX_DIR = \"$FindBin::Bin/../../matrix";
*MATRIX_FILE = \'blosum62.50';
*SEQ_DIR = \"$FindBin::Bin/../../seq";

*PHD_DIR = \"$ENV{HOME}/../torda/phd/results";


$COORD_DIR = "/bm/pdb_bin_pool";

#*COORD_DIR = \"$FindBin::Bin/../../struct";

*PARAM_DIR = \"$FindBin::Bin/../../params";
*FX9_PARAM_FILE = \'g9.param';
*RS_PARAM_FILE  = \'cc_allat_p891+0.param';
1;

