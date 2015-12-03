#!/usr/bin/perl

use FindBin;
use lib "$FindBin::Bin/../blib/lib";
use lib "$FindBin::Bin/../blib/arch";
use lib "$FindBin::Bin/../../../../../../src/Wurst/blib/lib";
use lib "$FindBin::Bin/../../../../../../src/Wurst/blib/arch";

use Wurst;
use Wurst::Pack qw(:all);
use POSIX qw(EXIT_SUCCESS EXIT_FAILURE);
use Storable qw (nfreeze thaw);
use Compress::Zlib;

use sigtrap qw (BUS SEGV PIPE ABRT); # try and get a stacktrace out of us.

$SIG{TRAP} = 'IGNORE'; kill 'TRAP', $$;

my ($seqnam, $sstr, $pstr, $bfile) = @ARGV;
my $Usage = "Usage : <seq file> <dssp file> <.bin file> <seqprof>\n";
(($seqnam) and (-f $seqnam)) or die "no sequence file\n$Usage";
((-f $pstr) ? $pstr : ((-f $pstr.".bin") ? ($pstr = $pstr.".bin") : 0)) or ( die "no .bin file\n$Usage");
(-f $sstr) or die "no dssp file\n$Usage";
(-f $bfile) or die "no blast file\n$Usage";
my $c = coord_read( $pstr ) || die "Failed on $pstr\n$Usage";
my $seq = seq_read($seqnam) || die "failed on $seqnam\n$Usage";
my $ssec = sec_s_data_read($sstr) || die "failed on $sstr\n$Usage";
my $sprof = blst_chk_read($bfile) || die "failed on $bfile\n$Usage";

my $m1 = sprintf "%s", ("Coord is called ".coord_name($c)."\nThe sequence is ".(seq_print($seq))."\nIts secondary structure :\n".(sec_s_data_string($ssec))."\nAnd the profile is \n".seqprof_str($sprof)."\n");

my $a = compress ( nfreeze( [coord_pack($c), sec_s_pack($ssec), seq_print($seq), seqprof_pack($sprof)] ) );
my $b = $a;
print "length : ".(length $a)."\n";
$c = thaw(Compress::Zlib::uncompress($b));
my $c2 = coord_unpack($c->[0]);
my $ssec2 = sec_s_unpack($c->[1]);
my $seq2 = seq_from_string($c->[2]);
my $sprof2 = seqprof_unpack($c->[3]);

my $m2 = sprintf( "%s", "Coord is called ".coord_name($c2)."\nThe sequence is ".(seq_print($seq2))."\nIts secondary structure :\n".(sec_s_data_string($ssec2))."\nAnd the profile is \n".seqprof_str($sprof2)."\n");
if ($m1 eq $m2) {
    print "Pack test was a success.\n";
} else {
    print "Pack test failed. Different data\n$m1\n$m2\n";
}
$c=\0;
$sprof2 = \0;
$seq2=\0;
$ssec2 = \0;
$c2 = \0;
$seq=\0;
$sprof=\0;
exit(EXIT_SUCCESS);
