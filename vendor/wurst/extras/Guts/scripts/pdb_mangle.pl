#!/usr/bin/perl

#Simple PDB Mangler

my $Usage = "pdb_mangle.pl <pdbfile> <ferocity of mangling 5..1> [<alternative filename>]\n";

use FindBin;
use lib "$FindBin::Bin/../src/Wurst/blib/arch";
use lib "$FindBin::Bin/../src/Wurst/blib/lib";

use strict;
use POSIX qw(EXIT_SUCCESS EXIT_FAILURE);

use Wurst;

# ----------------------- constants ---------------------------------
use vars qw ($pdbdir @pdbfiles @chains);

my ($o_struct, $f_struct,$prop,$o_seq);
((scalar @ARGV >3) or (scalar @ARGV < 2)) and die $Usage;
(-f $ARGV[0]) or die "$ARGV[0] nonexistent.\n";
($o_struct = (pdb_read $ARGV[0], ' ', ' ')) or die "Failed to read $ARGV[0] as a pdb file.\n".$Usage;
$o_seq = coord_get_seq($o_struct);
my ($sz, $iter, $ex_len);
$prop = $ARGV[1];
$prop =~ s/[^1-5]//;
(($prop>0) and ($prop<6)) or die "Out of Range for Mangling.\n";
$sz = coord_size $o_struct;
$ex_len = int($sz/(12.0-$prop));
$iter = (6-$prop);
while (($iter-->0) and ($sz>$ex_len)) {
  $f_struct = coord_deletion $o_struct, int(($sz-2*$ex_len)*rand), int($ex_len*rand), int($ex_len*rand);
  $sz = coord_size $f_struct;
  $ex_len = int($sz/(12.0-$prop));
  $o_struct = $f_struct;
}
(scalar @ARGV>2) and coord_2_pdb $ARGV[2], $f_struct, $o_seq;
(scalar @ARGV==2) and coord_2_pdb "m.".$ARGV[0], $f_struct, $o_seq;
