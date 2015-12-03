# This just tests reading of a sequence
# Pissy piss. This uses File::Temp which is not in older perl 5's.
use FindBin;
use lib "$FindBin::Bin../blib/arch";
use lib "$FindBin::Bin../blib/lib";
use POSIX qw(EXIT_SUCCESS EXIT_FAILURE);

use strict;

use File::Path;
use File::Temp qw/ tempfile tempdir /;
use Wurst;

my $seq1 = "
  > This is a nonsense sequence   
aaddaadd

";

my $seq2 = "
>more junk
aadgh
> another interesting sequence
xxxaa

";
my $seq3 = "> this is the third and most interesting sequence.
10 wwwww wwwww

";

# ----------------------- main --------------------------------------
sub main ()
{
    my $template = 'tmp_XXXXXXXX';

    my $dir = tempdir( CLEANUP => 0 );
    my ($fh, $filename) = tempfile( $template, SUFFIX => '.seq', DIR=> $dir);


    print "Writing temp sequence to $filename\n";
    print $fh $seq1;
    print "type something\n";
    my $crap;
    close ($fh);
    my $s = seq_read ($filename) || die "Failed to read seq1 from $filename";
    print "Sequence s is\n", seq_print ($s), "\n\n";
    my ($fh, $filename) = tempfile( $template, SUFFIX => '.seq', DIR=> $dir);
    print $fh $seq2; close ($fh);
    $s = seq_read_many ($filename) || die "Failed on seq2 from $filename";
    print "Second sequence \n", seq_print_many $s, "\n\n";
    rmtree ($dir, 0,1) || die "rmtree failed";
    $SIG{TRAP} = 'IGNORE'; kill 'TRAP', $$;
    my $s = seq_from_string ($seq3) || die "seq from string broken at";
    print "I have got a string from text\n", seq_print ($s), "\n\n";
    print "That is what it looks like, length ", seq_size ($s), "\n";
    return EXIT_SUCCESS;
}

exit main();
