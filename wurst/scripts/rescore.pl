# 9 April 2002
# Demonstrate the rescoring of a coordinate file.
# rcsid = $Id: rescore.pl,v 1.1 2007/09/28 16:57:02 mmundry Exp $
use strict;
use FindBin;
use lib "$FindBin::Bin../blib/arch";
use lib "$FindBin::Bin../blib/lib";

use vars qw ($rescore_paramfile $PARAM_DIR);
do "$FindBin::Bin/paths.inc";

my $scriptdir = $FindBin::Bin;

if ($@) {
    die "Error on include file for paths\n$@";
}

my $rescore_paramfile = "$PARAM_DIR/cc_allat_p891+0.param";
my $coord_file = "$scriptdir/../struct/1b4aA.bin";
my $pdb_file   = "$scriptdir/../struct/pdb1lea";

use POSIX qw(EXIT_SUCCESS EXIT_FAILURE);

use Wurst;
# ----------------------- mymain  -----------------------------------
sub mymain ()
{

    if (func_int() != 42) { die 'cannot find wurst function'; }
    my $params = param_rs_read ($rescore_paramfile) || die 'Boo';
    my $coord  = coord_read ($coord_file) || die;
    my $float  = score_rs( $coord, $params);
    print "For structure $coord_file, I scored $float\n";
    my $coord  = pdb_read ($pdb_file, '', '') || die;
    print "For structure $pdb_file, scored ", score_rs ($coord, $params), "\n";
    return EXIT_SUCCESS;
}
exit mymain();
