package Wurst::Pack;

# rcsid = "$Id: Pack.pm,v 1.1 2007/09/28 16:57:14 mmundry Exp $"

use 5.00004;
use strict;
use Carp;

require Exporter;
require DynaLoader;
use AutoLoader;
use vars qw($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS $AUTOLOAD);
@ISA = qw(Exporter
	DynaLoader);

%EXPORT_TAGS = ( 'all' => [ qw(
	coord_pack
	coord_unpack
	sec_s_pack
	sec_s_unpack
	seqprof_pack
	seqprof_unpack
) ] );

@EXPORT_OK = ( @{ $EXPORT_TAGS{'all'} } );

@EXPORT = qw(
	
);

$VERSION = '0.01';

sub AUTOLOAD {
    # This AUTOLOAD is used to 'autoload' constants from the constant()
    # XS function.

    my $constname;
    ($constname = $AUTOLOAD) =~ s/.*:://;
    croak "&Wurst::Pack::constant not defined" if $constname eq 'constant';
    my ($error, $val) = constant($constname);
    if ($error) { croak $error; }
    {
	no strict 'refs';
	# Fixed between 5.005_53 and 5.005_61
#XXX	if ($] >= 5.00561) {
#XXX	    *$AUTOLOAD = sub () { $val };
#XXX	}
#XXX	else {
	    *$AUTOLOAD = sub { $val };
#XXX	}
    }
    goto &$AUTOLOAD;
}

bootstrap Wurst::Pack $VERSION;

1;
__END__
=head1 NAME

Wurst::Pack -  Wurst objects across wires without nfs

=head1 SYNOPSIS

    use Wurst::Pack qw(:all);
    
=head1 ABSTRACT

    Functions of use when sending wurst through serialization 
    mediums like Storable (though you should compress the scalars
    they create).

=head1 DESCRIPTION
    
    These functions are slightly faster than writing to disk 
    and reading again.
    Generally, they are not safe across different native-binary 
    formats, and the format probably changed since you last
    checked out this module - BEWARE!

=head1 FUNCTIONS

=over

=item coord_pack COORD

Returns a scalar containing a binary form of the B<Coord> referred
to by COORD. This is convenient for storing a large number of
structures (such as a library) in a form that can be transmitted
to another process with a similar byte encoding form :


  my $packed_coord = coord_pack ($coord);

The structure can now be reliably stored for later recovery, via
a module such as Perl::Storable. See coord_unpack for decoding 
and limitations.

=item coord_unpack $packed_coord

$packed_coord is a scalar containing the data returned by a
prior call to B<coord_pack>. This function recovers a copy
of the original B<Coord> object that was packed.  There are 
certain limitations :


=over

- The B<Coord> is packed in the native binary format

- No serious checks are made to see if the scalar is garbage

- There is, however, a simple checksum.

=back

Typically, this is used to recover B<Coord> objects from
an array of scalars retrieved from a file, or sent via MPI.


=item sec_s_pack SEC_S

Returns a packed scalar, with a simple checksum, encapsulating
a Sec_s_data object.

=item sec_s_unpack Packed_SEC_S

Unpacks and returns a Sec_s_data object (as the Sec_s_dataPtr 
object).

  $x = sec_s_pack ( sec_s_data_read ("sec_data_file") );
  $y = sec_s_unpack($x);

=item seqprof_pack

Returns a packed scalar, with a very unreliable checksum,
encapsulating a sequence profile (B<Seqprof>)

=item seqprof_unpack

Recovers the B<Seqprof> object from a packed scalar.

    $x = seqprof_pack ( blst_chk_read ("blast_chktfile") );
    $y = seqprof_unpack($x);


=back

=head1 EXPORT

None by default. use qw (all:)

=head2 Exportable functions

=over
    
SCALAR coord_pack(B<Coord> coord)

B<Coord> coord_unpack(SCALAR pak_coord)

SCALAR sec_s_pack( B<Sec_s_data> msdat)

B<Sec_s_data> sec_s_unpack( SCALAR pb)

SCALAR seqprof_pack( B<Seqprof> prof)

B<Seqprof> seqprof_unpack( SCALAR pp)

=back

=head1 SEE ALSO

Wurst

=head1 AUTHOR

J.B. Procter

=head1 COPYleft AND LICENSE

2004 by James Procter

This library is part of Wurst.

=cut

