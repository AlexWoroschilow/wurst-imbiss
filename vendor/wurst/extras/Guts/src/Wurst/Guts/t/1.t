# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl 1.t'

#########################

# change 'tests => 1' to 'tests => last_test_to_print';

use Test;
BEGIN { plan tests => 3 };
use Wurst::Guts;
use FindBin;
use lib "$FindBin::Bin/../../../../../../src/Wurst/blib/lib";
use lib "$FindBin::Bin/../../../../../../src/Wurst/blib/arch";
use Wurst;
use POSIX qw(EXIT_SUCCESS EXIT_FAILURE);
my $Usage = "BROKEN INSTALLATION!";
my $pstr = "t/testbinf.bin";
my $bfile = "t/testprof.chk";
my $c = coord_read( $pstr ) || die "Failed on $pstr\n$Usage";
my $sprof = blst_chk_read($bfile) || die "failed on $bfile\n$Usage";




ok(1); # If we made it this far, we're ok.

#########################

# Insert your test code below, the Test::More module is use()ed here so read
# its man page ( perldoc Test::More ) for help writing this test script.

