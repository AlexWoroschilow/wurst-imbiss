# 31 Oct 2001
# For playing with making models.
use FindBin;
use lib "$FindBin::Bin../blib/arch";
use lib "$FindBin::Bin../blib/lib";

use strict;
use POSIX qw(EXIT_SUCCESS EXIT_FAILURE);

use Wurst;
use vars qw ($MATRIX_DIR $MATRIX_FILE $SEQ_DIR);

do "$FindBin::Bin/paths.inc" || die "Fail on paths.inc: $!\n";
*SEQ_DIR = \'../seq';

my $gap_open  = 3;
my $gap_widen = 2;

use vars qw($COORD_DIR);

my $coord_name = '1aac_.bin';

# ----------------------- mymain       ------------------------------
sub mymain ()
{
    $SIG{TRAP} = \&catch_trap;

    my $tmp_seq = 'crap_seq';
    my $name = $COORD_DIR . '/' . $coord_name;

    my $coord = coord_read ( $name) || die "coord_read fail on $name";

    my $c_seq = coord_get_seq ($coord) || die "No sequence in $coord";

    my $seq_string =
          "aivvdiakm
           kyetpelhvkvgdtvtwinreamphnvh
           wawawawa
           fvawawgeaawawawawawawaw";

    $seq_string = "vvdiakmkyetpelhvkvgdtv";
    my $seq0 = seq_from_string ($seq_string);

    my $scr_mat = score_mat_new (seq_size ($seq0), seq_size($c_seq));

    my $matname = $MATRIX_DIR . '/' . $MATRIX_FILE;
    my $m = sub_mat_read ($matname) || die "Fail getting sub mat $matname";
    undef ($matname);
    score_smat ($scr_mat, $seq0, $c_seq, $m);
    
    my $pair_set =
      score_mat_sum_smpl (my $result_mat, $scr_mat,
                          $gap_open, $gap_widen,
                          $gap_open, $gap_widen,
                          $S_AND_W);
    undef ($m);

    print pair_set_string ($pair_set, $seq0, $c_seq);
    print pair_set_pretty_string ($pair_set, $seq0, $c_seq);

#   Second, re-alignment

    my $pair_set2 =
      score_mat_sum_smpl (my $result_mat2, $scr_mat,
                          $gap_open, $gap_widen,
                          $gap_open, $gap_widen,
                          $N_AND_W, $pair_set);
    print "Second alignment\n",
    $DB::single = 1;
    print pair_set_pretty_string ($pair_set2, $seq0, $c_seq);
    my $new_model = make_model ($pair_set, $seq0, $coord);
    coord_2_pdb ("tmp.$$", $new_model) || die 'coord_2_pdb fail ';
    print "Please remove tmp.$$ after you are happy\n";
    return (EXIT_SUCCESS);
}
exit mymain();
