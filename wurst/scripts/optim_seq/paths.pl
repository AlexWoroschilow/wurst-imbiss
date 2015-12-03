# 6 March 2002
# $Id: paths.pl,v 1.1 2007/09/28 16:57:13 mmundry Exp $

*DFLT_BIN_DIR = ['/rsc/appenzeller/data1/pdblib6apr00',
                 '/rsc/appenzeller/data1/pdb_bin_pool' ];
*DFLT_PHD_DIR = ["$ENV{HOME}/phd/results"];

use FindBin;
my $where_I_am = $FindBin::Bin;
*MATRIX_DIR  = \"$where_I_am/../../matrix";
*MATRIX_FILE = \'blosum62.50';

*FX_PARAM_DIR  = \"$where_I_am/../../params";
*FX_PARAM_FILE = \'g7.param';
