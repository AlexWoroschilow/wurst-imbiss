# 11 Oct 2001
# Exercise the code for reading sequences.
use FindBin;
use lib "$FindBin::Bin../src/Wurst/blib/arch";
use lib "$FindBin::Bin../src/Wurst/blib/lib";

use strict;
use POSIX qw(EXIT_SUCCESS EXIT_FAILURE);

use Wurst;


use vars qw( $SEQ_DIR);
*SEQ_DIR = \'../seq';
my $fname1 = $SEQ_DIR . '/' . 'oneseq.seq';
my $fname2 = $SEQ_DIR . '/' . 'twoseq.seq';
my $fname3 = $SEQ_DIR . '/' . 'pti.seq';

my $SEQ_DIR = '../seq';

sub catch_trap { return; }


# ----------------------- junk         ------------------------------
sub junk ()
{
    my $devnull = '/dev/null';
    open (CRAP, ">$devnull") || die "dev null failed $!";
    my $t = seq_read_many ($fname3) || die "Fail on fname3\n";
    my $t = seq_read_many ($fname1, $t) || die "Fail on fname1\n";
    print "t has ", seq_num ($t), " sequences\n";

    for (my $j = 0; $j < 3; $j++) {
        for (my $i = 0; $i < seq_num ($t); $i++) {
            my $x = seq_get_1 ($t, $i);
            print CRAP "hello seq $i has ", seq_size ($x), " members\n";
            my $copy = $x;
            print CRAP "copy is $copy, x is $x\n";
            undef $copy;
            print CRAP "size seq checks original seq, ", seq_size ($x);
        }
    }
    close (CRAP);
}

# ----------------------- mymain       ------------------------------
sub mymain ()
{
    $SIG{TRAP} = \&catch_trap;

    my $t = seq_read_many ($fname1) || die "Fail on $fname1";

    $t = seq_read_many ($fname2, $t) || die "broke reading $fname2";
    for ( my $i = 0; $i < 2; $i++) {
        junk();
    }
    my $seq0 = seq_get_1 ($t, 0);
    my $seq1 = seq_get_1 ($t, 1);

    return (EXIT_SUCCESS);
}
exit mymain();
