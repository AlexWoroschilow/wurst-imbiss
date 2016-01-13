#!/usr/bin/perl
use strict;
use warnings;

my $file = $ARGV[0] or die "Need to get CSV file on the command line\n";

open( INPUT, "<$ARGV[0]" ) or die "Could not open '$file' $!\n";

my $first = 1;

while ( my $line = <INPUT> ) {
	chomp $line;
	if ( $first == 1 ) {
		$first = 0;
		next;
	}
	my @fields = split ";", $line;

	my $seq1     = $fields[2];
	my $seq2     = $fields[3];
	my $tm_score = $fields[4];
	print "$seq1\t$seq2\t$tm_score\n";
}
