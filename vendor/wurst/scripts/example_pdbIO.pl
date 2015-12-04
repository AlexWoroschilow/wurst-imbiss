=pod

=head1 NAME

example_pdbIO.pl

=head1 SYNOPSIS

example_pdbIO.pl <filename1> <filename1>

=head1 DESCRIPTION

This is a simple script demonstrating wurst's pdb reading
and writing capabilities.
It parses a pdb-file and compiles it into another file.

=head1 PARAMETERS

=over

=item filename1

A .pdb file that should be read

=item filename2

Name of a file to write in the script's output

=over

=head1 EXAMPLE


=item perl example_pdbIO.pl wurstfolder/testdata/1H1W.pdb /somewhere/out.pdb

=over

=cut


use FindBin;
use lib "$FindBin::Bin/../blib/arch/";
use lib "$FindBin::Bin/../blib/lib/";

use strict;
use POSIX qw(EXIT_SUCCESS EXIT_FAILURE);
use Wurst;

# ----------------------- mymain       ------------------------------
sub mymain ()
{ 
    if ($#ARGV != 1) {
        print "Usage: perl $0 </path/to/pdbfile.pdb> </path/to/pdbfile.pdb>\n";
        return EXIT_FAILURE;
    }
    # here a pdb file is parsed by wurst
    # and stored in the Coord data type
    my $Coord = pdb_read ($ARGV[0], '', '') || die "pdb_read fail: $!";
    # now the Coord data will be written back in a pdb file
    coord_2_pdb ($ARGV[1], $Coord) || die "Failed writing parsed PDB data";
    return EXIT_SUCCESS;
}
exit mymain();
