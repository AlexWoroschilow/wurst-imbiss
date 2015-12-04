#!/usr/bin/perl

# PDB Cutter
#Simpler PDB Mangler

my $Usage = "pdb_cutter.pl <pdbfile> <i|d|a> <%age excision> <alternative filename>\n";

use FindBin;
use lib "$FindBin::Bin/../src/Wurst/blib/arch";
use lib "$FindBin::Bin/../src/Wurst/blib/lib";

use strict;
use POSIX qw(EXIT_SUCCESS EXIT_FAILURE);

use Wurst;

# ----------------------- constants ---------------------------------
use vars qw ($pdbdir @pdbfiles @chains);
my %Mode = ( I => 1, i => -1, A => 0, a => 0, D => 1, d => 1);
my ($o_struct, $f_struct,$prop,$mode);
((scalar @ARGV != 4)) and die $Usage;
(-f $ARGV[0]) or die "$ARGV[0] nonexistent.\n";
($o_struct = (pdb_read $ARGV[0], ' ', ' ')) or die "Failed to read $ARGV[0] as a pdb file.\n".$Usage;

my ($sz, $iter, $ex_len);
$mode = $ARGV[1];
$mode =~ s/[^iIdDaA]//;
($mode eq '') and die "Mode of cutting : I/i Insertion from sequence, \nI/i deletion from sequence, A/a affine gap (non-alignment excision).\n$Usage";
$mode = $Mode{$mode};
$prop = $ARGV[2];
$prop =~ s/[^0-9.]//;
$sz = coord_size $o_struct;
$ex_len = int($sz*($prop/100.0));
srand (-$prop);
$iter = int(($sz-2*$ex_len)*rand);

  ($mode == -1) and ($f_struct = coord_deletion $o_struct, $iter, 0, $ex_len);
  ($mode == 1) and ($f_struct = coord_deletion $o_struct, $iter, $ex_len, 0);
  ($mode == 0) and ($f_struct = coord_deletion $o_struct, $iter, $ex_len, $ex_len);

coord_2_pdb $ARGV[3], $f_struct, coord_get_seq ($f_struct);
print "[$ARGV[1], $iter, $ex_len]\n";

