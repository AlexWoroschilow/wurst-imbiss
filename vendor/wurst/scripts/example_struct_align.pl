=pod

=head1 NAME

example_struct_align.pl

=head1 SYNOPSIS

example_struct_align.pl <class file> <pdb file> <pdb file>
                        <pdb file to create> <pdb file to create> <SW or NW>

=head1 DESCRIPTION

This script demonstrates wurst's structure alignment capabilities.

=head1 EXAMPLE

perl example_struct_align.pl wurstdir/extras/testfiles/classfile 1H1W.pdb
                             1PHK.pdb 1H1W_aligned.pdb 1PHK_aligned.pdb SW

=cut

use FindBin;
use lib "$FindBin::Bin/../blib/arch/";
use lib "$FindBin::Bin/../blib/lib/";

use Wurst;


my $m_shift=-0.1;
my $pgap_open=3.25;
my $pgap_widen=0.8942;
my $qgap_open=3.25;
my $qgap_widen=0.8942;


# ----------------------- mymain       ------------------------------
sub mymain ()
{
    if ($#ARGV !=5) {
     print "Usage:\n",
            "perl $0 <class file> <pdb file> <pdb file>",
                    "<pdb file to create> <pdb file to create> <SW or NW>\n";
     print "Example:\n",
           "perl $0 wurstdir/extras/testfiles/classfile 1H1W.pdb 1PHK.pdb",
                   "1H1W_aligned.pdb 1PHK_aligned.pdb SW\n";
     return EXIT_FAILURE;
    }

    my $classfcn = aa_strct_clssfcn_read ($ARGV[0], 0.4);
    my $coord1=pdb_read ($ARGV[1], '', '');
    my $coord2=pdb_read ($ARGV[2], '', '');
    my $alignmenttype=$S_AND_W;
    if ($ARGV[5] eq "NW")
    {
        $alignmenttype=$N_AND_W;
    }
    print "Calculating Probability vectors. This may take a while...\n";
    my $pvec1 = strct_2_prob_vec ($coord1, $classfcn);
    my $pvec2 = strct_2_prob_vec ($coord2, $classfcn);
    print "Finished\n";
    # create a score matrix
    my $matrix = score_mat_new (prob_vec_length($pvec1),
                                prob_vec_length($pvec2));
    # fill out score matrix
    score_pvec ($matrix, $pvec1, $pvec2);
    $matrix = score_mat_shift ($matrix, $m_shift);
    $sw_pair_set = score_mat_sum_smpl (
        my $crap_mat,   $matrix, $pgap_open, $pgap_widen,
        $qgap_open, $qgap_widen, $alignmenttype
    );
    # superimpose coordinates
    ($rmsd, $new_c1, $new_c2) =
        coord_rmsd ($sw_pair_set, $coord1, $coord2, 0);

    print "RMSD is: $rmsd\n";
    print "Writing aligned structures to $ARGV[3] and $ARGV[4]\n";
    coord_2_pdb ($ARGV[3], $new_c1);
    coord_2_pdb ($ARGV[4], $new_c2);
}
exit mymain();