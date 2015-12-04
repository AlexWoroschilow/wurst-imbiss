#!/usr/bin/perl
# Eat a line from the MC output and print it as an Andrew format
# substitution matrix
# Since it is only one line, read from standard input ?

=pod

=head1 NAME

line2mat.pl [ I<file> ]

=head1 SYNOPSIS

B<line2mat.pl> S<[ B<-p> I<prec> ]> S<[ B<-w> I<width>]> S<[ B<-s> I<matshift>]> S<[ B<-t> I<matscale> ]>

=head1 DESCRIPTION

Take a line from a minimisation calculation like MC or simplex
and turn it into a substitution matrix. If a filename is given on
the command line, the first line will be read, otherwise the
standard output is read.

When written out, two manipulations are always done

=over

=item *

We reorder the amino acids so they come out in the same style as
the blosum matrices distributed with fasta and blast.

=item *

We invent the missing amino acids.

=over

=item B

The average of asp and asn types.

=item Z

The average of glu and gln types.

=item X

The average over all types

=item '*'

Translation stop. This is the only fake residue which is not
calculated in some way. It is simply set at -1.

=back

=back

=head1 OPTIONS

=over

=item B<-p> [I<num>]

The precision of of the matrix as written. This should be
S<C<-p 0>>
if the matrix will be used by blast or fasta.
This should be an integer.

=item B<-w> [I<num>]

Set the width of field in the output. By setting it wide enough,
the matrix can be written out in a form easy for humans to read.
This should be an integer.

=item B<-s> [I<num>]

Scale the matrix by I<num> by simple multiplication.

=item B<-t> [I<num>]

Add I<num> to each element of the matrix before scaling.

=back

=head1 EXAMPLE

 perl line2mat.pl -p 0 -w 3 -s 100 -t 0.6 < file_in > matrix.out

would add 0.6 to all the matrix elements and then multiply them
by 100. This would be useful for keeping some of our precision,
but in a form suitable for use in fasta.

=cut

use strict;
undef $@;
use warnings;
if (defined ($@)) { warn "Not using warning module $@\n" };
use warnings;
use POSIX qw(EXIT_SUCCESS EXIT_FAILURE);
# rcsid = $Id: line2mat.pl,v 1.1 2007/09/28 16:57:01 mmundry Exp $
# ----------------------- constants ---------------------------------
use vars qw ($N_LETTER);
*N_LETTER = \20;

# ----------------------- get_nums ----------------------------------
sub get_nums ($)
{
    my $infile = shift;
    my $line = <$infile>;
    if ( ! defined ($line) ) {
        warn "broke reading input line\n"; return (); }
    my @nums = split (' ', $line);
    my $diff = $#nums - 210 + 1;
    for ( my $i = 0; $i < $diff; $i++ ) {
        shift @nums ;}
    return @nums;
}

# ----------------------- num2mat  ----------------------------------
# Take the one dimensional array of numbers and spread it into a
# two dimensional array for printing out
sub num2mat (\@ \@)
{
    my ($nums, $mat) = @_;
    my $n = 0;
    for (my $i = 0; $i < $N_LETTER; $i++) {
        for (my $j = $i; $j < $N_LETTER; $j++) {
            $$mat[$i][$j] = $$nums [$n++]; } }
    for (my $i = 0; $i < $N_LETTER; $i++) {
        for (my $j = $i + 1; $j < $N_LETTER; $j++) {
            $$mat[$j][$i] = $$mat[$i][$j];} }
}

# ----------------------- mat2hash ----------------------------------
sub mat2hash (\% \@ \@)
{
    my ($mathash, $mat, $letters) = @_;
    for (my $i = 0; $i < 20; $i++) {
        for ( my $j = 0; $j < 20; $j++) {
            my $a = $$letters[$i];
            my $b = $$letters[$j];
            $$mathash {"$a$b"} = $$mat[$i][$j];
            $$mathash {"$b$a"} = $$mat[$j][$i];
        }
    }
}

# ----------------------- add_crap ----------------------------------
# Add information to the hash about the non-real residues
# There is the possibility of errors here, so return undef if we
# have problems.
# Remember B = asp/asn, Z = glu/gln, X = any, * = termination
sub add_crap (\% \@)
{
    my ($mathash, $letters) = @_;
    my @fakes = ('B', 'Z', 'X', '*');
    foreach my $i (@fakes) {
        if ( defined ( $$mathash {"$i$i"})) {
            warn "Elements with $i defined. Should not be\n"; return undef; }}
    my $ntot = 0;
    my $tot = 0.0;
    foreach my $a (@$letters) {
        my $X = 0.0;
        foreach my $b (@$letters) {
            $X += $$mathash {"$a$b"}; }
        $tot += $X; $ntot++;
        $$mathash{"${a}X"} = $$mathash{"X$a"} = $X / @$letters;
    }
    $$mathash {"XX"} = $tot / ($ntot * @$letters);
    foreach my $a (@$letters, 'X', 'B') {   # define the 'B' as av of B and N
        my $B =  ($$mathash {"${a}D"} + $$mathash {"${a}N"}) / 2.0;
        $$mathash {"B$a"} = $$mathash {"${a}B"} = $B;
    }
    foreach my $a (@$letters, 'X', 'B', 'Z') { # define 'Z' av of E and Q
        if (!defined ($$mathash {"${a}E"})) {
            $DB::single = 1; }
        if (!defined ($$mathash {"${a}Q"})) {
            $DB::single = 1; }
        my $Z =  ($$mathash {"${a}E"} + $$mathash {"${a}Q"}) / 2.0;
        $$mathash {"Z$a"} = $$mathash {"${a}Z"} = $Z;
    }
    foreach my $a (@$letters, 'X', 'B', 'Z') {
        $$mathash {"${a}*"} = $$mathash {"*${a}"} = -4; }
    $$mathash {'**'} = 1;
    return 1;
}

# ----------------------- write_out ---------------------------------
sub write_out ($ \% \@ $ $ $ $ )
{
    my ($outfile, $mathash, $outletters, $width, $prec, $matshift, $matscale)
        = @_;
    print $outfile '# Matrix written at ', scalar (localtime()), "\n";
    
    if ( defined ($matshift) ) {
        print $outfile "# Matrix shifted from original by $matshift\n"; }
    if ( defined ($matscale) ) {
        print $outfile "# Matrix scaled from original by $matscale\n"; }
    if ( !defined ($width)) {
        $width = 2; }
    if ( ! defined ($prec)) {      
        $prec  = 1; }

    print $outfile ' ';
    my $swidth;
    if ($width < ($prec)) {
        $swidth = $width + $prec; }
    else {
        $swidth = $width}
    for (my $i = 0; $i < @$outletters; $i++) {
        printf $outfile " %${swidth}s", $$outletters[$i]; }
    print $outfile "\n";
    foreach my $a (@$outletters) {
        print $outfile $a;
        foreach my $b (@$outletters) {
            if (!defined ($$mathash{"$a$b"})) {
                print STDERR "\n$a $b undef\n" };
            printf $outfile " %${width}.${prec}f", $$mathash{"$a$b"}; }
        print $outfile "\n";
    }
}

# ----------------------- scale_mat ---------------------------------
sub scale_mat ( \% \@ $ $)
{
    my ( $mathash, $outletters, $matshift, $matscale) = @_;
    if (defined ($matshift)) {
        foreach my $a (@$outletters) {
            foreach my $b (@$outletters) {
                $$mathash {"$a$b"} += $matshift; } } }

    if (defined ($matscale)) {
        foreach my $a (@$outletters) {
            foreach my $b (@$outletters) {
                $$mathash {"$a$b"} *= $matscale; } } }
}

# ----------------------- usage    ----------------------------------
sub usage ()
{
    print STDERR "\
$0 [-p precision] [-s matrix_scale] [-t matrix_shift ] [-w width] [infile]
width and precision refer to printf statement in matrix printing\n";
}
# ----------------------- mymain   ----------------------------------
sub mymain ()
{
    use Getopt::Std;

    my ($precision, $matscale, $matshift, $width);
    my %opts;
    if (! getopts ('p:s:t:w:', \%opts)) {
        usage(); return EXIT_FAILURE; }

    if (defined ($opts {p})) { $precision = $opts {p};}    
    if (defined ($opts {s})) { $matscale  = $opts {s};}
    if (defined ($opts {t})) { $matshift  = $opts {t};}
    if (defined ($opts {w})) { $width     = $opts {w};}
    
    my $infile;

    if ( @ARGV == 0 ) {
        $infile = \*STDIN; }
    else {
        if ( !open ( INFILE , "<$ARGV[0] ")) {
            warn "open fail on $ARGV[0]\n"; return EXIT_FAILURE;}
        $infile = \*INFILE;
    }

    my (@letters, %letters, @outletters);
    @letters = ( 'G', 'A', 'V', 'L', 'I', 'F', 'P', 'S', 'T', 'C',
                 'M', 'W', 'Y', 'N', 'Q', 'D', 'E', 'K', 'R', 'H');
    @outletters = (
        'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I',
        'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V',
        'B', 'Z', 'X', '*'
    );

    %letters = (
        G => 0,        A => 1,        V => 2,        L => 3,
        I => 4,        F => 5,        P => 6,        S => 7,
        T => 8,        C => 9,        M => 10,       W => 11,
        Y => 12,       N => 13,       Q => 14,       D => 15,
        E => 16,       K => 17,       R => 18,       H => 19
    );

    my @nums = get_nums ($infile);
    if ( ! @nums) {
        return EXIT_FAILURE;}
    close ($infile);
    my @mat;
    num2mat (@nums, @mat);
    if (@nums < 209) {
        warn "Too few numbers found\n";
        return EXIT_FAILURE;
    }
    my $outfile = \*STDOUT;
    my %mathash;
    mat2hash (%mathash, @mat, @letters);
    if ( ! add_crap (%mathash, @letters )) {
        warn "Problem in add_crap\n"; return EXIT_FAILURE; }

#    write_out ($outfile, @nums, @mat, @letters);
    scale_mat ( %mathash, @outletters, $matshift, $matscale);
    write_out ($outfile, %mathash, @outletters, $width, $precision, $matshift, $matscale);
    $DB::single = 1;
    return EXIT_SUCCESS;
}

# ----------------------- main      ---------------------------------
exit (mymain());
