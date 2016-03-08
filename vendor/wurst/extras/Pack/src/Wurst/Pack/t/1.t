# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl 1.t'

#########################

# change 'tests => 1' to 'tests => last_test_to_print';

use Test;
BEGIN { plan tests => 5 };
use Wurst::Pack;
ok(1); # If we made it this far, we're ok.

#########################

# Insert your test code below, the Test::More module is use()ed here so read
# its man page ( perldoc Test::More ) for help writing this test script.
use FindBin;
use lib "$FindBin::Bin/../../../../../../src/Wurst/blib/lib";
use lib "$FindBin::Bin/../../../../../../src/Wurst/blib/arch";
use Wurst;
use Wurst::Pack qw(:all);
use POSIX qw(EXIT_SUCCESS EXIT_FAILURE);
use Storable qw (nfreeze thaw);
use Compress::Zlib;
my $Usage = "BROKEN INSTALLATION!";
my $seqnam = "t/test.seq";
my $sstr = "t/testsecs.phd";
my $pstr = "t/testbinf.bin";
my $bfile = "t/testprof.chk";
my $c = coord_read( $pstr ) || die "Failed on $pstr\n$Usage";
my $seq = seq_read($seqnam) || die "failed on $seqnam\n$Usage";
my $ssec = sec_s_data_read($sstr) || die "failed on $sstr\n$Usage";
my $sprof = blst_chk_read($bfile) || die "failed on $bfile\n$Usage";


my $m1 = ("Coord is called ".coord_name($c),
          "The sequence is ".(seq_print($seq)),
          "Its secondary structure :\n".(sec_s_data_string($ssec)),
          "And the profile is \n".seqprof_str($sprof)."\n");

my $a = compress ( nfreeze( [coord_pack($c), sec_s_pack($ssec), seq_print($seq), seqprof_pack($sprof)] ) );
my $b = $a;
print "length : ".(length $a)."\n";
$c = thaw(Compress::Zlib::uncompress($b));
my $c2 = coord_unpack($c->[0]);
my $ssec2 = sec_s_unpack($c->[1]);
my $seq2 = seq_from_string($c->[2]);
my $sprof2 = seqprof_unpack($c->[3]);
my $m2 = ("Coord is called ".coord_name($c2),
          "The sequence is ".(seq_print($seq2)),
          "Its secondary structure :\n".(sec_s_data_string($ssec2)),
          "And the profile is \n".seqprof_str($sprof2)."\n");
my @t = ("pack_coord", "seq_print","sec_s_pack","seqprof_pack");
while (scalar @t) {
    my $tt = shift @t;
    (print "Failed $tt\n") 
        unless (ok((shift @m2) eq (shift @m1)));
}

        
    
