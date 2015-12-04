=pod

=head1 NAME

example_seq_align.pl

=head1 SYNOPSIS

example_seq_align.pl <aSequence> <aSequence> <path to matrixfile>

=head1 DESCRIPTION

This script demonstrates wurst's sequence alignment capabilities.

=head1 EXAMPLE

perl example_seq_align.pl acdefg gggaaa /wurstfolder/testdata/blosum62.mat

=cut

use lib "$FindBin::Bin../blib/arch";
use lib "$FindBin::Bin../blib/lib";

use strict;
use POSIX qw(EXIT_SUCCESS EXIT_FAILURE);
use Wurst;

my $gap_open1  = 40;
my $gap_widen1 = 11;
my $gap_open2  = $gap_open1 / 10;
my $gap_widen2 = $gap_widen1 / 10;


# ----------------------- string_to_seq  ----------------------------
# Write the string given as an argument to a temporary file,
# read it up using wurst's seq_read function and return the
# sequence.
sub string_to_seq ( $ )
{
    my $string = shift;
    my $tmpfile = POSIX::tmpnam();
    open (T, ">$tmpfile") || die "Open fail on $tmpfile";
    print T  $string;
    close (T) || die "close fail";
    my $t = seq_read ($tmpfile) || die "Fail on $tmpfile";
    unlink ($tmpfile) || die "Failed to unlink $tmpfile";
    return $t;
}


# -----------------------   do_align   ------------------------------
# Aligns two sequences
# Parameters:
# 1 sequence1
# 2 sequence2
# 3 alignment type (either $N_AND_W or $S_AND_W)
# 4 The Substitution Matrix
sub do_align ($ $ $ $)
{
    my ($s1, $s2, $al_type, $subst_mat) = @_;
    print "seq1 is: ", seq_print ($s1), 
          "gaps are: $gap_open1, $gap_widen1\n\n";
    print "seq2 is: ", seq_print ($s2),
          "gaps are: $gap_open2, $gap_widen2\n\n";
    
    my $scr_mat = score_mat_new (seq_size ($s1), seq_size ($s2)) || die
        "score_mat_new fail: $!";

    score_smat ($scr_mat, $s1, $s2, $subst_mat);
    $scr_mat = score_mat_shift ($scr_mat, 3);
    print score_mat_string ($scr_mat, $s1, $s2);

    my $result_mat;
    my $pair_set =
        score_mat_sum_smpl ($result_mat, $scr_mat,
                            $gap_open1, $gap_widen1,
                            $gap_open2, $gap_widen2,
                            $al_type);

    print pair_set_pretty_string ($pair_set, $s1, $s2);
    print "Short extension:\n";
    pair_set_extend ($pair_set, seq_size ($s1), seq_size ($s2), $EXT_SHORT);
    print pair_set_pretty_string ($pair_set, $s1, $s2);
    print "Long extension:\n";
    pair_set_extend ($pair_set, seq_size ($s1), seq_size ($s2), $EXT_LONG);
    print pair_set_pretty_string ($pair_set, $s1, $s2);
}


# ----------------------- mymain       ------------------------------
sub mymain ()
{
    if ($#ARGV !=2) {
     print "Usage:   perl $0 <aSequence> <aSequence> <path to matrixfile>\n";
     print "Example: perl $0 acdefg gggaaa /wurstfolder/matrix/blosum62.mat\n";
     return EXIT_FAILURE;
    }
    
    my $seq1=string_to_seq($ARGV[0]);
    my $seq2=string_to_seq($ARGV[1]);
    
    my $matname = $ARGV[2];
    my $subst_mat = sub_mat_read ($matname) || die
        "Fail getting substitution matrix: $matname";
    
    
    print ("Alignment Needleman and Wunsch:\n");
    do_align($seq1, $seq2, $N_AND_W, $subst_mat);
    print ("Alignment by Smith and Waterman:\n");
    do_align($seq1, $seq2, $S_AND_W, $subst_mat);
    
    return (EXIT_SUCCESS);
}
exit mymain();
