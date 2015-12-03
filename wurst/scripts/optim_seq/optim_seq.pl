#!/usr/bin/perl
# 15 Dec 2001
# For optimising alignment parameters.
# rcsid = "$Id: optim_seq.pl,v 1.1 2007/09/28 16:57:14 mmundry Exp $"
=pod

=head1 NAME

optim_seq.pl - Optimise parameters for alignment

=head1 SYNOPSIS

perl S<[ B<-v> ] optim_seq.pl pairfile param_file align_func_file

=head1 DESCRIPTION

The idea is to optimise parameters for sequence to structure alignment.
The basic idea is that we have a list of pairs of proteins.
Within each pair, the proteins have some amount of structural
similarity.

Align the sequence of the first on to the structure of the
second and vice versa. From each alignment, construct a model.
Compare this model against the correct answer for the sequence.
Adjust the parameters so as to get the best possible alignment.

There are two ways to approach this

=over

=item *

We could run over the whole set of proteins, read up everyone's
coordinates once. At each step of the optimisation, use the stored
coordinates and sequences.  This needs enough memory to store lots of
structures and their sequences.

=item *

Alternatively, we could be very economical with memory. Read a
pair of proteins, do the alignment, calculate the cost
function from this pair.  Throw away the coordinates and move
on to the next pair.  This means we do not have to store a lot
of coordinates.

=back

=head2 OPTIONS

This script does not encourage options since other big scripts
have degenerated into a mess of cryptic one-letter codes
necessary to get things to run.  The list as it stands is

=over

=item B<-v>

Be more verbose.  You might see the following..

=over 8

=item *

After reading the list of proteins in the pairfile, print them
out. This is useful to find the line where an error occurred.

=back

=back

=head2 INPUT

There are some input files

=over

=item pair file

This argument is compulsory.
The format looks like

 # A comment here.
 1abc_ 2xyxA
 3qrs_ 4abc_

Where the first four letters are the PDB code and the fifth
character is a chain name.

=over 8

=item *

The protein names may be separated by any amount of white
space.

=item *

Anything after a hash C<#> is treated as a comment and
ignored.

=item *

Leading and trailing white space are stripped and thrown away.

=back

=item parameter file

The rules are different. Each line begins by a keyword,
followed by a comma separated list of values. There is no
comma after the keyword name.

The comma is important. It allows us to have white space embedded
in a value. In this example,

  # comments are allowed and ignored
  names gap_open, gap_widen, other param
  lower 0.1, 0.001, -3
  upper 10, 100, 3
  start 0.5, 0.8, 1

the last name will be "S<other param>".
So, the keywords are

=over 8

=item names

These are names for the parameters. The code does not
understand them, but will print them out at the end of a
calculation.

=item start

Starting values or initial guess. This is compulsory.

=item lower

Lower bounds. This is optional. If used, the  number of values
must be the same as for the starting array.

=item upper

Upper bounds. This is the same as for lower bounds.

=back

=item align_func_file

This script calls functions like B<init_align_cost> from
somewhere else, depending on what is being optimised. These
functions are provided by a file which does the work. This file
(without the .pm extension) must be provided on the command line.

=back

=cut

use strict;

eval 'use warnings';
if ($@) {
    print "Perl warnings not in use\n"; }

use POSIX qw(EXIT_SUCCESS EXIT_FAILURE);

# ----------------------- Defaults and constants --------------------
use vars qw ($DFLT_MAX_ITER $DFLT_MAX_RESTART $DFLT_F_TOL );

*DFLT_MAX_ITER    = \500;
*DFLT_MAX_RESTART = \5;
*DFLT_F_TOL       = \10e-5;

# ----------------------- Prototypes    -----------------------------
sub init_align_cost (\@ \% $);

# ----------------------- get_line     ------------------------------
# Given a reference to a file handle, read from it and return a line.
# Might as well allocate the line here.
# Hop over blank lines.
# Delete anything after a hash (#) comment marker.
# Delete leading and trailing white space.
# Return undef() when we hit the end.
sub get_line ( $ )
{
    my $f = shift;
    while ( defined (my $t = <$f>)) {
        $t =~ s/#.*//;     # Remove comments
        $t =~ s/^\s*//;    # and leading blanks
        $t =~ s/\s+$//;    # and trailing blanks
        if (length ($t) == 0) {              # Blank ? keep going
            next; }
        return $t;
    }
    return undef();
}

# ----------------------- get_pair_list -----------------------------
sub get_pair_list ($)
{
    if (! open (PAIRS, "<$_[0]") ) {
        print STDERR "Open fail on $_[0]: $!\n"; return 0 }


    my @pairlist;
    while ( defined (my $line = get_line (\*PAIRS))) {
        my @t = split ' ', $line;
        if ($#t < 1) {
            print STDERR "Broken line in pairlist \"$line\"\n";
            return 0;
        }
        if ($#t >1) {    # Trim off any other info that could be in file
            $#t = 1; }
        push (@pairlist, @t);
    }

    close (PAIRS);
    return @pairlist;
}

# ----------------------- parameters_sane ---------------------------
# Check if we have essential parameters and anything else
# before starting a run.
sub parameters_sane ($)
{
    my $s_args = shift;
    my $result = EXIT_SUCCESS;
    my $n_value = $#{$$s_args{ini_pt}};
    if ( defined $$s_args{names}) {
        if ((my $t = $#{$$s_args{names}}) != $n_value) {
            print STDERR "Wrong number of names ", $t + 1, "\n";
            $result = EXIT_FAILURE;
        }
    }

    if ( defined $$s_args{lower}) {
        if ((my $t = $#{$$s_args{lower}}) != $n_value) {
            print STDERR "Wrong number of lower bounds ", $t + 1, "\n";
            $result = EXIT_FAILURE;
        }
    }

    if ( defined $$s_args{upper}) {
        if ((my $t = $#{$$s_args{upper}}) != $n_value) {
            print STDERR "Wrong number of upper bounds ", $t + 1, "\n";
            $result = EXIT_FAILURE;
        }
    }

    if ( $result == EXIT_FAILURE ) {
        print STDERR "Reading parameters, from the initial point array,\n",
        "I would like ", $n_value + 1, " values\n";
    }
    return $result;
}

# ----------------------- get_param_list ----------------------------
sub get_param_list ($ $ $)
{
    my ($filename, $s_args, $fixed) = @_;

    my $result = EXIT_SUCCESS;
    if (! open (PARAMS, "<$filename" )) {
        print STDERR "Open fail on $filename: $!\n"; return EXIT_FAILURE }

    while (defined (my $line = get_line (\*PARAMS))) {
        my ($name, $val) = split '[\s,]+', $line, 2;
#       my @val = split ' *, *', $val; # old, space tolerant version
        my @val = split '[\s,]+', $val;
        $name =~ tr/A-Z/a-z/;
      SWITCH: {
          if ($name eq 'names') { $$s_args{names}  = \@val; last SWITCH; }
          if ($name eq 'start') { $$s_args{ini_pt} = \@val; last SWITCH; }
          if ($name eq 'upper') { $$s_args{upper}  = \@val; last SWITCH; }
          if ($name eq 'lower') { $$s_args{lower}  = \@val; last SWITCH; }
          if ($name eq 'o_file')  { $$s_args{$name} = $val; last SWITCH;}
          if ($name eq 'max_iter'){ $$s_args{$name} = $val; last SWITCH;}
          if ($name eq 'max_restart') {
              $$s_args{$name} = $val; last SWITCH;}
          if ($name eq 'f_tol'){ $$s_args{$name} = $val; last SWITCH;}
          if ($name eq 'scatter'){ $$s_args{$name} = $val; last SWITCH;}
          if ($name eq 'fixed')   { $$fixed{$val[0]} = $val[1]; last SWITCH;}
          {print STDERR "Param file unknown line \"$line\"\n";}
      }
    }

    $result = parameters_sane ($s_args);  # Check if input params are sensible
    return $result;
}

# ----------------------- usage        ------------------------------
sub usage {
    print STDERR "usage: $0 [ -v ] pair_list_file param_file align_fnc_file\n";
    exit (EXIT_FAILURE);
}


# ----------------------- check_and_default -------------------------
# We will repeat this operation a few times, so it gets its
# own function.
# We have our argument hash, a possible parameter and its
# default.
# If the parameter has not been set, set it to its default.
sub check_and_default ($ $ $)
{
    my ($s_args, $param, $dflt) = @_;
    if ( ! defined ($$s_args {$param})) {
        $$s_args {$param} = $dflt; }
}

# ----------------------- crap_func    ------------------------------
sub crap_func ($)
{
    my ($a, $b) = @_;
    return ( ($a -5) * ($a - 5) + ($b + 3) * ($b + 3));
}


# ----------------------- optimise     ------------------------------
# Return EXIT_FAILURE or EXIT_SUCCESS after we try to optimise the
# function.

sub optimise ($ $) {
    use FindBin;
    use lib "$FindBin::Bin/../../src/Simplex";
    use Simplex;

    my $s_args = shift;
    my $verbosity = shift;

    check_and_default ( $s_args, 'max_iter', $DFLT_MAX_ITER);
    check_and_default ( $s_args, 'max_restart', $DFLT_MAX_RESTART);
    check_and_default ( $s_args, 'f_tol', $DFLT_F_TOL);

    $$s_args{func} = \&align_cost;

    if ($verbosity >= 1) {
        my $fmt = "%15s %15s \n";
        print "Run starting at ", scalar(localtime()), "  With:\n";
        printf ($fmt, 'Max iterations', $$s_args{max_iter});
        printf ($fmt, 'Max restarts',   $$s_args{max_restart});
        printf ($fmt, 'Output file', $$s_args{o_file});
        printf ($fmt, 'tolerance', $$s_args{f_tol});

        print "Parameters to minimise: \n";
        my $n = $#{$$s_args{names}};
        my ($low, $up) = undef;
        if (defined $$s_args{lower}) {
            $low = 1;}
        if (defined $$s_args{upper}) {
            $up = 1;}
        $fmt = "%15s";
        my $gfmt = "%15g";
        printf ($fmt, 'parameter');
        printf ($fmt, 'initial');
        if ($low) {
            printf ($fmt, 'lower'); }
        if ($up) {
            printf ($fmt, 'upper'); }
        print "\n";
        for (my $i = 0; $i <= $n; $i++) {
            printf ($fmt,  $ {$$s_args{names}}[$i]);
            printf ($gfmt, $ {$$s_args{ini_pt}}[$i]);
            if ($low) {
                printf ($gfmt, $ {$$s_args{lower}}[$i]); }
            if ($up) {
                printf ($gfmt, $ {$$s_args{upper}}[$i]); }
            print "\n";
        }
    }

    return (simplex ($s_args));
}

# ----------------------- catch_int    ------------------------------
# The main thing is, if we get a KILL or TERM, to call exit and get
# out of here. This means there is a better chance of closing files
# wherever we were up to.
sub catch_kill
{
    print STDERR "Exiting on signal\n";
    exit EXIT_FAILURE;
}

# ----------------------- mymain       ------------------------------
sub mymain
{
    use Getopt::Std;
    use FindBin;
    use lib "$FindBin::Bin";

    my $verbosity = 1;
    my %opts;
    if ( ! getopts ('v:', \%opts)) {
        usage(); }
    if (defined ($opts {v})) {
        $verbosity = $opts {v}; }
    undef %opts;


    if ( $#ARGV < 2) {
        print STDERR "Not enough arguments\n";
        usage();
    }

    my $pair_file = $ARGV[0];
    my $param_file = $ARGV[1];
    my $align_func_file = $ARGV[2];
    eval "use $align_func_file";
    if ($@) {
        die "Failed on eval of $align_func_file\n$@";
    }


    my %s_args;               # This is where simplex args will go
    my %fixed;                # Where we put fixed paramters for the alignments


    $s_args{verbosity} = $verbosity;
    print "Pairs coming from $pair_file\n";
    my @pairlist = get_pair_list($pair_file);
    if ( ! @pairlist) {
        print STDERR "Fail to get pairs\n";
        return EXIT_FAILURE;
    }

    if ($verbosity > 4) {
        for (my $i = 0; $i <= $#pairlist; $i++) {
            print "$pairlist[$i]";
            if ($i % 2 != 0) {
                print "\n";
            } else {
                print ' ';
            }
        }
    }

    print "Initial parameters from $param_file\n";
    if (  get_param_list ($param_file, \%s_args, \%fixed) == EXIT_FAILURE ) {
        print STDERR "Broke reading parameters from $param_file\n";
        return EXIT_FAILURE;
    }
    if (init_align_cost (@pairlist, %fixed, $s_args{names}) == EXIT_FAILURE){
        print STDERR "Error setting up for cost function\n";
        return EXIT_FAILURE;
    }
    my $result;
    $SIG{INT } = \&catch_kill;
    $SIG{QUIT} = \&catch_kill;
    $SIG{TERM} = \&catch_kill;
    if (($result = optimise (\%s_args, $verbosity)) == EXIT_FAILURE) {
        print STDERR "Optimisation FAILED !!!!\n";
        print STDERR "Dodgy values...";
    }

    print "Final results\n",
          "-------------\n";
    print "Num cycles $$result{ncycle}\n",
          "Final func value $$result{value}\n",
          "Best parameters were \n";

    my $best = $$result{best};
    my $nn;
    if (defined ($nn = $s_args{names})) {
        for (my $i = 0; $i < $#{$$result{best}} + 1; $i++) {
            printf "%20s %10g\n", $$nn[$i], $$best[$i]; } }
    else {
        for (my $i = 0; $i < $#{$$result{best}} + 1; $i++) {
            printf "%10g\n", $$best[$i];}}

    return $$result{success};
}

exit (mymain);
