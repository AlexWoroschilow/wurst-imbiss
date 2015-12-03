package Wurst::Guts;

use 5.00004;
use strict;
use Carp;

require Exporter;
require DynaLoader;
use AutoLoader;
use vars qw($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS $AUTOLOAD);
@ISA = qw(Exporter
	DynaLoader);

# Items to export into callers namespace by default. Note: do not export
# names by default without a very good reason. Use EXPORT_OK instead.
# Do not simply export all your public functions/methods/constants.

# This allows declaration	use Wurst::Guts ':all';
# If you do not need this, moving things directly into @EXPORT or @EXPORT_OK
# will save memory.
%EXPORT_TAGS = ( 'all' => [ qw(
	coord_deletion
        coord_segment
        coord_merge
	mci_contact_rs
	mci_sc_n_rs
	pair_set_xchange
	scor_set_rs
	seq_deletion
        pair_set_trim
        pair_set_shift
        seqprof_merge
        seqprof_trim
        ps_cmp_and_model
) ] );

@EXPORT_OK = ( @{ $EXPORT_TAGS{'all'} } );

@EXPORT = qw(
             coord_segment
             seqprof_trim
             coord_merge
             seqprof_merge
);

$VERSION = '0.01';

sub AUTOLOAD {
    # This AUTOLOAD is used to 'autoload' constants from the constant()
    # XS function.

    my $constname;
    ($constname = $AUTOLOAD) =~ s/.*:://;
    croak "&Wurst::Guts::constant not defined" if $constname eq 'constant';
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

bootstrap Wurst::Guts $VERSION;

# Preloaded methods go here.

# Autoload methods go after =cut, and are processed by the autosplit program.

1;
__END__
# Below is stub documentation for your module. You'd better edit it!

=head1 NAME

Wurst::Guts - Perl extension for blah blah blah

=head1 SYNOPSIS

  use Wurst::Guts;
  blah blah blah

=head1 ABSTRACT

  This should be the abstract for Wurst::Guts.
  The abstract is used when making PPD (Perl Package Description) files.
  If you don't want an ABSTRACT you should also edit Makefile.PL to
  remove the ABSTRACT_FROM option.

=head1 DESCRIPTION

Stub documentation for Wurst::Guts, created by h2xs. It looks like the
author of the extension was negligent enough to leave the stub
unedited.

Blah blah blah.

=head2 EXPORT

None by default.

=head2 Exportable functions

  struct coord *coord_deletion(struct coord *, size_t, size_t,size_t)
  int mci_contact_rs( const char *fname, struct coord *coord, const float *P)
  int mci_sc_n_rs( const char *fname, struct scor_set *scrs, struct coord *coord, const float *P)
  struct pair_set *pair_set_xchange (struct pair_set *p_s)
  struct scor_set *scor_set_rs (struct coord *coord, const float *P)
  struct seq *seq_deletion(struct seq *, size_t, size_t,size_t)
 
=item coord_deletion COORD Start Coord_Length Seq_Length

Returns a pointer to a new B<Coord> (structure) where an
affine gap has been introduced, arising from the deletion
of some sequence or coordinate data from an existing structure.
The gap is specified by its B<Start> [0,size), and number
of residues excised from the structure (B<Coord_Length>)
and sequence (B<Seq_Length>).

Minimal checking ensures gap specifications are sensible.
This function is provided for the generation of synthetic
threading results and other forms of improper structure 
sets.

=item coord_merge COORD1 COORD2

returns a concatenation of COORD1 and COORD2.

=item seqprof_merge SEQPROF1 SEQPROF2

returns a concatenation of SEQPROF1 and SEQPROF2.


=item  pair_set_trim PAIR_SET Start 
 * extract a range of pairs aligning to a specified range of
 * sequence, and return them as a new pairset.
 * end=0 and we define the range from start->end of sequence
 * start = 0 and we define the range from start->end

  pair_set_shift

 * Adds a signed constant to the sequence entry of a pair_set 
 * object.


  seqprof_trim
 * extract a range of a blast sequence profile
 * as a new profile.


=item mci_contact_rs MCI_FILE COORD PARAMS

Writes the contact map as defined by the tanH score-function (as
used for rescoring a model) for the model in B<COORD>. Parameters 
are given by PARAMS, and the string in MCI_FILE is used as the 
filename.

=item mci_sc_n_rs MCI_FILE SCOR_SET COORD PARAMS

Writes a map of contact scores (from tanh PARAMS) modulated by the
local sequence structure fitness (SCOR_SET) for the alignment
that generated the model COORD. The calculation is rather ad-hoc:

 cm[i][j] = tanh_score(i,j)
            + Sum(i from 1 to N_res) (
               Sum (j from i+1 to N_res) ( 
                 local_fitness(i) + local_fitness(j)
              ) )
  


=item ps_cmp_and_model PAIRSET_A2B SCORE_MATRIX_A2B COORDA PAIR_SET_B2A SCORE_MATRIX_B2A COORDB

B<Returns> a list :

 CoordPtr(s): Cons(SeqAStrA), Cons(A2B), Cons(B2A), Cons(SeqBStrB),
 String:Fragments, Pair_setPtr:Conserved_Pairs.

Enables the comparison of structures via threading, as used in
the the B<Seq_Str_Cmp> module.

Given two structures (A and B), other functions can be used 
to generate a threading (A2B) of sequence A on structure B, 
and the converse (B2A).  This function can then be used to 
compare the two alignments, finding stretches of each 
structure which are consistently aligned by the threading 
method. 

The consistent structure is identified as a new pair set
(B<Conserved_Pairs>) where the scores indicate the 'simple' 
score matrix summation for the conserved stretches. Using 
B<pair_set_score> on this returns a list giving the score of 
sequence B on structure A, followed by the score of 
sequence A on structure B. This abuse of nomenclature allows 
standard rescoring on the conserved parts of each threading, 
via B<pair_set_xchange> (because the gapped score is used
in these functions).

In addition to the pair set, a string is returned, giving the
conserved pair alignments as a fragment list of the form :

    (startposA, startposB, length)   # length gt 1
    (startposA, startposB)           # length == 1

The excessive number of models that are returned allows the
geometric assessment of the conserved regions between the
structures.  The first two models are the conserved stretches
of sequence A, overlayed on their original coordinates, and
their co-locations in structure B. The second two are
respectively the co-location of the conserved region of
sequence B on structure A, and the original coordinates of
the conserved regions.

See the module B<Seq_Str_Cmp> in 'Wurst::scripts' for the 
one and only use of this function.


=item pair_set_xchange PAIR_SET

Returns a new pair set, where the entries of pairs have been
reversed. Wurst pair sets are assumed to indicate a sequence's
alignment to a structure. If one wishes to analyse the
structural properties of an alignment from the sequence's
perspective, then this function can be used to reverse the
original sense of a threading.

=item scor_set_rs MODEL PARAMS

Returns a B<Scor_set> object where the score at each site in
MODEL is the sum of tanh-forcefield scores involving that 
position. (see scoranlys.c for implementation)



=head1 SEE ALSO

Mention other useful documentation such as the documentation of
related modules or operating system documentation (such as man pages
in UNIX), or any relevant external documentation such as RFCs or
standards.

If you have a mailing list set up for your module, mention it here.

If you have a web site set up for your module, mention it here.

=head1 AUTHOR

James Procter, E<lt>procter@suse.deE<gt>

=head1 COPYRIGHT AND LICENSE

Copyright 2004 by James Procter

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself. 

=cut
