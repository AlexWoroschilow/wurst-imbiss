# 13 Dec 2001
# rcsid = $Id: simplex.pl,v 1.1 2007/09/28 16:57:02 mmundry Exp $

# This is for testing the simplex code in Simple.pm.
# If you use the thing definition of toyfunc() below, then you
# might also like to take the following lines and feed them to
# gnuplot.
=pod
f(x) = (x + 8) * (x+ 8) + (x-40) * (x-40) + 30 * x * sin(x)

set xran[-20:20]

set arrow from -4.27916,-500 to -4.27916,5000
set arrow from 5.18,-500 to 5.18,5000
set arrow from 17.32,-500 to 17.32,5000
plot f(x)
pause 15

=cut
# They will draw a nice picture of the function below, at least
# with respect to the $a variable.

use lib "$FindBin::Bin../blib/arch";
use lib "$FindBin::Bin../blib/lib";
use Simplex;
use strict;
use POSIX qw(EXIT_SUCCESS EXIT_FAILURE);
# ----------------------- toyfunc -----------------------------------
sub toyfunc ( @)
{
    my ($a, $b, $c) = @_;
    my $r = ($a + 8) * ($a + 8) + ($a - 40) * ($a - 40) + 30 * $a * sin($a)
        + $b * $b * $b * $b
        + ($c - 5.5) * ($c - 5.5);
    return $r;
}

# ----------------------- mymain       ------------------------------
sub mymain
{
    my $fref = \&toyfunc;
    my @guess = (2, 14, -2);
    my @lower = (-10, -3000, -5);
    my @upper = (20, 230, 20);
    my $max_iter = 500;
    my $max_restart = 5;
    my $f_tol = 10e-5;

    my %result;
    my %s_arg = (
        func        => \&toyfunc,
        ini_pt      => \@guess,
        lower       => \@lower,
        upper       => \@upper,
        max_iter    => $max_iter,
        max_restart => $max_restart,
        o_file      => 'splx_out',
        scatter     => 0.20,
        f_tol       => $f_tol,
        result      => \%result
                 
    );
    
    my $result = simplex (\%s_arg);
    
    if ($$result {success} == $SPLX_SUCCESS) {
        print "Simplex happy \n";
    } elsif ($$result {success} == $SPLX_TOO_MANY) {
        print "Simplex too many\n";
    } elsif ($$result {success} == $SPLX_BROKEN) {
        print STDERR "Simplex broken";
    } else {
        die "Programming bug ";
    }
    for (my $i = 0; $i < $s_arg{n_dim}; $i++) {
        printf '%4g ', "${${$s_arg{result}}{best}}[$i]"; }
    print "\n";

#   Alternative, readable form
    my $success = $$result{success};
    my @best = @{${$s_arg{result}}{best}};
    my $best_value = $$result {value};
    print "best value of $best_value at \n@best\n";
    return (EXIT_SUCCESS);
}
exit (mymain);

