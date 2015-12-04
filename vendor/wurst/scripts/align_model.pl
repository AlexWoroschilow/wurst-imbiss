# 21 dec 2001
# For playing with making models.
use strict;
eval 'use warnings';
if ($@) {
    print "Warnings not present\n";}

use FindBin;

my $paths_place = 'paths.inc';
do $paths_place || die "Failed to read $paths_place: $!";
if ($@) { die "Compilation fail on $paths_place:\n$@" };

use vars qw(@WURST_DIRS $N_AND_W);
my $tdir = "$FindBin::Bin/../src/Wurst";
*WURST_DIRS = ["$tdir/arch",
               "$tdir/lib" ];
# But after installation
# WURST_DIRS = [ "$ENV{HOME}/pl/lib" ];

use lib "$FindBin::Bin../blib/arch";
use lib "$FindBin::Bin../blib/lib";
#use lib "$ENV{HOME}";  # Where the wurst stuff lives after installation
use Wurst;

use POSIX qw(EXIT_SUCCESS EXIT_FAILURE);


use vars qw ($MATRIX_DIR $MATRIX_FILE $SEQ_DIR $PHD_DIR);

my $gap_open  = 5;
my $gap_widen = 2;

use vars qw($COORD_DIR);


my $coord_name = '1b4aA.bin';
my $seqcoord   = '1aoy_';

my $out_pdb = 'foo.pdb';

# ----------------------- mymain       ------------------------------
sub mymain ()
{
    my $name = $COORD_DIR . '/' . $coord_name;
    my $coord = coord_read ( $name) || die "coord_read fail on $name";
    my $c_seq = coord_get_seq ($coord) || die "No sequence in $coord";
    print "seq1 is \n", seq_print (coord_get_seq($coord)), "\n";
    my $x = $COORD_DIR . '/' . $seqcoord . '.bin';
    my $cseq  = coord_read ($x) || die "coord fail on $x ";
    my $seq0 = coord_get_seq ($cseq);
    print "seq0 is \n", seq_print ($seq0), "\n";

    my $phdfile = "$PHD_DIR/$seqcoord.phd";
    print "phdfile is $phdfile\n";
    my $sec_s_data = sec_s_data_read ($phdfile)
        || die "sec_s_data failed on $phdfile";
    print "Sec data 1\n", sec_s_data_string ($sec_s_data);

    my $scr_mat = score_mat_new (seq_size ($seq0), seq_size($c_seq));

    my $matname = $MATRIX_DIR . '/' . $MATRIX_FILE;
    my $m = sub_mat_read ($matname) || die "Fail getting sub mat $matname";

    score_smat ($scr_mat, $seq0, $c_seq, $m);
    
    my $pair_set =
      score_mat_sum_smpl (my $result_mat, $scr_mat,
                          $gap_open, $gap_widen,
                          $gap_open, $gap_widen,
                          $N_AND_W);
    undef ($matname);
    undef ($m);

#   print pair_set_pretty_string ($pair_set, $seq0, $c_seq, $sec_s_data,$coord);
    pair_set_extend ($pair_set, seq_size ($seq0), coord_size ($coord), $EXT_LONG);
    print "Now, extended a bit \n",

    pair_set_pretty_string  ($pair_set, $seq0, $c_seq);
#   $SIG{TRAP} = 'IGNORE'; kill 'TRAP', $$;
    my $new_model = make_model ($pair_set, $seq0, $coord);
    coord_2_pdb ("tmp.$$", $new_model) || die 'coord_2_pdb fail ';
    print "Please remove tmp.$$ after you are happy\n";
    return (EXIT_SUCCESS);
}
exit mymain();
