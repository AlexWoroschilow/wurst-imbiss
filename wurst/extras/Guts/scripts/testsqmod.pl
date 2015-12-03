#!/usr/bin/perl

use FindBin;
use lib "$FindBin::Bin/../src/Wurst/blib/arch";
use lib "$FindBin::Bin/../src/Wurst/blib/lib";

use strict;
use POSIX qw(EXIT_SUCCESS EXIT_FAILURE);

use Wurst;

use Seq_Str_Cmp;

$SIG{TRAP} = 'IGNORE'; kill 'TRAP', $$;
my ($fname1,$fname2); 

Set_Verbosity(100);

$fname1 = shift @ARGV;
$fname2 = shift @ARGV;
(-f $fname1 and -f $fname2) or die "Give me two binary filenames\n";
my $c1 = coord_read ($fname1) || die "Fail on $fname1";
my $c2 = coord_read ($fname2) || die "Fail on $fname2";

my @resset = x_sq_struct_align($c1,$c2);

# do stuff with results
my @legset = x_sq_struct_legend();
while (scalar @legset) {
  print "".(shift @legset)."\t".(shift @resset)."\n";
};

#Make models.
#Reconstruct the alignment from the coverage
#calculate %age id. and anything else.

exit ((scalar @resset > 4) ? EXIT_SUCCESS : EXIT_FAILURE);


