# 21 Dec 2001
# For trying out sequence to structure alignments.

use FindBin;
use lib "$FindBin::Bin../blib/arch";
use lib "$FindBin::Bin../blib/lib";

use strict;
use POSIX qw(EXIT_SUCCESS EXIT_FAILURE);

use Wurst;

# ----------------------- Constants    ------------------------------
use vars qw ($GAP_OPEN $GAP_WIDEN);
*GAP_OPEN  = \1.0;
*GAP_WIDEN = \0.5;

use vars qw ($STRUCT_DIR $MATRIX_DIR $MATRIX_FILE);
*STRUCT_DIR = \'../bin_tmp';
*MATRIX_DIR = \'../matrix';
*MATRIX_FILE= \'blosum62.50';



# Some files for playing with
my $pdb1 = '1b4aA.bin';        # get sequence from here
my $pdb2 = '1aa0_.bin';        # This is hand made for playing with

# ----------------------- mymain       ------------------------------
sub mymain ()
{
    $SIG {TRAP} = 'IGNORE';#    kill 'TRAP', $$;

    my $seq1 = seq_from_string ('lvvl');

    my $c1 = coord_read ("$STRUCT_DIR/$pdb2") || die "Open fail on $pdb2 ";
    my $seq2 = coord_get_seq ($c1);


#   Now get a substitution matrix
    my $matname = "$FindBin::Bin/$MATRIX_DIR/$MATRIX_FILE";
    my $subst_mat = sub_mat_read ($matname) || die "Fail reading $matname ";

    my $scr_mat = score_mat_new (seq_size ($seq1), coord_size ($c1)) || die ;

    $scr_mat = score_smat ($scr_mat, $seq1, $seq2, $subst_mat)
        || die "score error ";
#    score_cmpct ($scr_mat, $seq1, $c1, $param) || die "compact error ";
                              
    $DB::single = 1;

    print "seq size is ", seq_size ($seq1),
    " coord size is ", coord_size ($c1), "\n",
    "Sequence is\n", seq_print (coord_get_seq ($c1)), "\n\n";

    my $sec_pnlty = 0.1;
    my $pair_set = score_mat_sum_sec (
                                          my $result_mat, $scr_mat,
                                          $c1, $sec_pnlty,
                                          $GAP_OPEN, $GAP_WIDEN,
                                          $GAP_OPEN, $GAP_WIDEN,
                                          $S_AND_W);


    print "now, with the fancy sec struct penalties\n";
    pair_set_extend ($pair_set, seq_size ($seq1), seq_size ($seq2), $EXT_LONG);
    print pair_set_string ($pair_set, $seq1, coord_get_seq ($c1));
#   my $crap_coord = "tmp_$$.pdb";
#   coord_2_pdb ($crap_coord, make_model ($pair_set, $seq1, $c1));
#   print "Please rm $crap_coord after you are happy\n";
    return (EXIT_SUCCESS);

}

exit (mymain ());
