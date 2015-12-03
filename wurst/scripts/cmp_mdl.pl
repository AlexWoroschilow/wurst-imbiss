# 31 Oct 2001
# For playing with making models.
use FindBin;
use lib "$FindBin::Bin../src/Wurst/blib/arch";
use lib "$FindBin::Bin../src/Wurst/blib/lib";

use strict;
use POSIX qw(EXIT_SUCCESS EXIT_FAILURE);

use Wurst;
use vars qw ($MATRIX_DIR $MATRIX_FILE $SEQ_DIR);
*MATRIX_DIR = \'../matrix';
*MATRIX_FILE= \'blosum62.50';

*SEQ_DIR = \'../seq';

my $gap_open  = 3;
my $gap_widen = 2;

use vars qw($COORD_DIR);
*COORD_DIR = \"../bin_tmp";

my $coord_name = '1aa0_.bin';

my $out_pdb = 'foo.pdb';

# ----------------------- catch_trap   ------------------------------
sub catch_trap { return; }

# ----------------------- mymain       ------------------------------
sub mymain ()
{
    $SIG{TRAP} = \&catch_trap;

    my $tmp_seq = 'crap_seq';
    my $name = $COORD_DIR . '/' . $coord_name;

    my $coord = coord_read ( $name) || die "coord_read fail on $name";

    my $c_seq = coord_get_seq ($coord) || die "No sequence in $coord";

    my $seq_string =
          "dkatipsespfaaaevadgaivvdiakm
           kyetpelhvkvgdtvtwinreamphnvh
           fvagvlgeaalkgpmmkkeqaysltfteagtydyhctphpfmrgkvvve\n";
    $seq_string =
          "> moo cow
           ka cti p\n";

    my $seq0 = seq_from_string ($seq_string);


    my $scr_mat = score_mat_new (seq_size ($seq0), seq_size($c_seq));

    my $matname = $MATRIX_DIR . '/' . $MATRIX_FILE;
    my $m = sub_mat_read ($matname) || die "Fail getting sub mat $matname";

    score_smat ($scr_mat, $seq0, $c_seq, $m);


    my $result_mat;
    my $pair_set =
      score_mat_sum_smpl ($result_mat, $scr_mat,
                          $gap_open, $gap_widen,
                          $gap_open, $gap_widen,
                          $N_AND_W);
    undef ($matname);
    undef ($m);

    print pair_set_string ($pair_set, $seq0, $c_seq);
    kill 'TRAP', $$;
    my $new_model = make_model ($pair_set, $seq0, $coord);
#   coord_2_pdb ("tmp.$$", $new_model) || die 'coord_2_pdb fail ';
#   print "Please remove tmp.$$ after you are happy\n";
    my $frac;
    dme_thresh ($frac, $new_model, $coord, 0.4) || die 'dme_thresh error ';
    print "frac is $frac\n";
    return (EXIT_SUCCESS);
}
exit mymain();
