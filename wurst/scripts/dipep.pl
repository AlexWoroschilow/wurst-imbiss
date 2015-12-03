# 8 Feb 2005
# First script to check the dipeptide alignment code.
# We use the rather old interface to temp file functions so the
# script will work with a very old perl (5.0).
# rcsid = $Id: dipep.pl,v 1.1 2007/09/28 16:57:01 mmundry Exp $

use FindBin;
use lib "$FindBin::Bin../blib/arch";
use lib "$FindBin::Bin../blib/lib";

use strict;
use POSIX qw(EXIT_SUCCESS EXIT_FAILURE);

use Wurst;

# ----------------------- Constants ---------------------------------
my $seq2 = seq_from_string ('rarara');
my $seq1 = seq_from_string ('grrawrqggc');

# ----------------------- dipeptide info ----------------------------
# This script is for testing and debugging, so let's bundle up
# the necessary components. We used to use a separate file, but
# lets put our test data in here.
sub
write_test_file ()
{
    use Fcntl;
    use POSIX qw (tmpnam);
    $DB::single = 1;
    my $text = "
# 10 feb 2005 for testing dipeptide information
aw aw 4
# this is just a comment
gc aa 2 1.1
aa gw 4 2.2
ac gg 1 3.3
ac de 3 4.4
ac df 3 5.5
ac de 1 6.6
";

    my $old_tmp = $ENV{TMPDIR};
    $ENV{TMPDIR} = '.';
    my $fname = tmpnam();
    $ENV{TMPDIR} = $old_tmp;
    if ( ! open (TMP, ">$fname")) {
        warn "Failed to make temp file $fname\n";
        return undef;
    }
    print TMP $text;
    close (TMP);
    return $fname;
}


# ----------------------- test_get_set ------------------------------
# Test the dpt_get_val/dpt_set_val functions
sub
test_get_set ($)
{
    my $dipep = shift;
    my $i;
    my $ret = 1;
    for ($i = 0; $i < dpt_get_n($dipep); $i++) {
        my $f = dpt_get_val ($dipep, $i);
        if (!$f) {
            warn "bug at i = $i. val should not be undefined\n"; $ret = undef}
        print "dpt list entry $i is $f\n";
    }
    my $f = dpt_get_val ($dipep, dpt_get_n($dipep));
    if ($f) {
        warn  "I got a value \'$f\' for a non-existing value.\n"; $ret = undef}
    else {
        print "Checked.. Getting beyond dipep data produces an error.\n"}

    for ($i = 0; $i < dpt_get_n($dipep); $i++) {
        if (!dpt_set_val ($dipep, $i, -$i)) {
            warn "bug at i = $i with dpt_set_val()\n"; $ret = undef} }

    print "Have set all the values in the dpt list.\n";
    if (!dpt_set_val ($dipep, dpt_get_n($dipep), 0.0)) {
        print "Checked.. setting beyond dipep data produces an error.\n"}
    else {
        warn "Bug: dpt_get_val did not flag invalid number\n"; $ret = undef}
    print "After setting all values, I now have \n", dpt_string ($dipep), "\n";
    return $ret;
}

# ----------------------- test_align --------------------------------
# To convince us that one can do alignments with the dipeptide
# code. To play here, try setting the gap open/widen values to
# big or small numbers.
sub
test_align ($ $ $)
{
    my ($dipep, $seq1, $seq2) = @_;
    my $scr_mat = score_mat_new (seq_size ($seq1), seq_size($seq2)) ||
        die "score_mat_new() fail: $!\n";
    score_dpt ($scr_mat, $seq1, $seq2, $dipep);
    print "After dipep score, matrix is\n",
           score_mat_string ($scr_mat, $seq1, $seq2), "\n";
    my $result_mat;
    my $al_type = 'S_AND_W';
    my ($gap_open1, $gap_widen1) = (1.0, 0.1);
    my ($gap_open2, $gap_widen2) = (1.0, 0.1);

    my $pair_set =
        score_mat_sum_smpl ($result_mat, $scr_mat,
                            $gap_open1, $gap_widen1,
                            $gap_open2, $gap_widen2,
                            $al_type);
    print "We have calculated the following alignment\n";
    print pair_set_pretty_string ($pair_set, $seq1, $seq2);
    pair_set_extend ($pair_set, seq_size ($seq1), seq_size ($seq2), $EXT_LONG);
    print "After sequence extension, we get\n";
    print pair_set_pretty_string ($pair_set, $seq1, $seq2);
    return 1;
}

# ----------------------- mymain  -----------------------------------
sub
mymain ()
{
    my $dipep;
    my $dipep_dat = write_test_file ();
    $SIG{TRAP} = 'IGNORE';
    if (! $dipep_dat) {
        return EXIT_FAILURE; }
    if ( ! ($dipep = dpt_read ($dipep_dat))) {
        warn "Fail reading $dipep_dat\n"; return EXIT_FAILURE; }
    unlink ($dipep_dat);
    print "Just read ", dpt_get_n ($dipep),
    " lines of information: \n", dpt_string ($dipep);
    test_align($dipep, $seq1, $seq2);
    if ( ! test_get_set($dipep)) {
        warn "Bugs in get/set routines\n"; return EXIT_FAILURE; }

    return EXIT_SUCCESS;
}

exit mymain();
