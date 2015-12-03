# 4 March 2002
# This provides the cost function for the alignment optimisation.
# There are some conventions.
# We provide a cost function, so lower is better.
# We calculate the amount of distance matrix preserved after
# applying a threshold. This is then used as input to an
# activation function (logistic). The sign is changed so it
# runs from -1..0.  -1 means very good.
# If we reject a model, we should return zero, which means bad.
#
# This will use Thomas's fx9 set, but with a similarity matrix
# and shifting of both fx9 and similarity matrices. If
# sec_s_weight is non-zero, it will use that too.
# $Id: Align_fx9_sim.pm,v 1.1 2007/09/28 16:57:14 mmundry Exp $
use strict;
package Align_fx9_sim;
BEGIN {
    use vars qw(@ISA @EXPORT);
    use Exporter ();
    @ISA = qw(Exporter);
    @EXPORT = qw(&align_cost &init_align_cost);
}

use lib "$ENV{HOME}/pl/lib";
use Wurst;
use POSIX qw(EXIT_SUCCESS EXIT_FAILURE);
use vars qw ($first);
$first = 1;
# ----------------------- Defaults and Constants --------------------
use vars qw(@DFLT_BIN_DIR @DFLT_PHD_DIR);
use vars qw ($MATRIX_DIR $MATRIX_FILE );
use vars qw ($FX_PARAM_DIR  $FX9_PARAM_FILE );
my $path_file = 'paths.pl';
if ( ! do $path_file ) {
    print STDERR "Error on $path_file: $!\n";
    exit (EXIT_FAILURE);
}

use vars qw ($MIN_MODEL $SWITCH_POINT $X_SCALE);
*MIN_MODEL = \0.5;     # Model must account for this fraction of residues
*SWITCH_POINT = \0.8;  # Activation function switches at this point
*X_SCALE      = \20.0; # Scaling for fraction for activation function

# ----------------------- Global variables --------------------------

use vars qw(@coord $sub_mat $fx_params @sec_info);
use vars qw(@sec_scr_mats @fx_scr_mats);
use vars qw($count);
$count = 0;

# ----------------------- Global parameters -------------------------
# These are separated from the variables above, since they
# will come in from the outside world.  They will be
# parameters which may or may not be optimised.
use vars qw( $pgap_open $pgap_widen $qgap_open $qgap_widen
             $sec_pnlty $sec_s_scale $sim_scale
             $sim_mat_bottom $fx_mat_shift);
$pgap_open      = 10;
$pgap_widen     = 1;
$qgap_open      = 10;
$qgap_widen     = 1;
$sec_pnlty      = 5.0;
$sec_s_scale    = 1.0;
$sim_scale      = 1.0;
$sim_mat_bottom = -3.7;
$fx_mat_shift   = 1;

use vars qw( $align_type);
$align_type = $S_AND_W;
use vars qw( @param_name);


# ----------------------- init_align_cost ---------------------------
# There are some things to do before calculating alignment costs.
# This routine should be called once, before align_cost().
# Things we need
#  * protein list
#  * substitution matrix
#  * set up fixed parameters - do this via soft references
# Philosophy:
# We do not bail out when we first encounter an error. It is
# more useful to check as many files as possible and print out
# as many mistakes as we can find.  We do, however, return an
# error so the caller knows to give up.
sub init_align_cost (\@ \% $)
{
    my ($pairlist, $fixed, $param_name) = @_;
    my $err = 0;
    my @bin_dirs = @DFLT_BIN_DIR;
    my @phd_dirs = @DFLT_PHD_DIR;
    my $tmp_p = "$FX_PARAM_DIR/$FX9_PARAM_FILE";

    PROT: foreach my $p (@$pairlist) {
        my $found = 0;
        foreach my $d (@phd_dirs) {
            my $path = "$d/$p.phd";
            if ( -e $path) {
                $found = 1;
                if (! (my $s = sec_s_data_read ($path))) {
                    print STDERR "sec_s_data error on $path\n";
                    $err = 1;
                } else {
                    push (@sec_info, ($s));
                    next PROT;
                }
            }
        }
        if ( ! $found) {
            print STDERR "Could not find $p.phd\n";
            $err = 1;
        }
    }
    PROTBIN: foreach my $p (@$pairlist) {
        my $found = 0;
        foreach my $d (@bin_dirs) {
            my $path = "$d/$p.bin";
            if ( -e $path) {
                $found = 1;
                if (! (my $c = coord_read ($path))) {
                    print STDERR "coord_read error on $path\n";
                    $err = 1;
                } else {
                    push (@coord, ($c));
                    next PROTBIN;
                }
            }
        }
        if ( ! $found) {
            print STDERR "Could not find $p\n";
            $err = 1;
        }
    }
    if ( !($fx_params = param_fx_read ($tmp_p))) {
        print STDERR "Fail reading parameters from $tmp_p\n";
        $err++;
    }

    @param_name = @$param_name;
    {
        no strict 'refs';
        foreach my $key (keys (%$fixed)) {
            foreach my $p (@param_name) {
                if ($key eq $p) {
                    print STDERR "$key is both fixed and variable. Fix\n";
                    $err++;
                }
            }
            $$key = $$fixed{$key}; }
        for (my $i = 0; $i <= $#param_name; $i++) {
            if (! defined ($ {$param_name[$i]})) {
                print STDERR "Tried to set variable $param_name[$i]\n",
                "This is unknown. Check list of parameters at start of file\n",
                "Make sure each is declared and initially defined\n";
                $err++;
            }
        }
        my $fmt  = "%15s %15s\n";
        my $gfmt = "%15s %15g\n";
        print "Parameters fixed:\n";
        printf ($fmt, 'parameter', 'value');
        no strict 'refs';
        foreach my $key (keys (%$fixed)) {
            printf ($gfmt, $key, $$key);
        }

        use strict;
    }

#   Special handling of this one with a symbolic name
    if ($align_type eq "n_and_w") {
        $align_type = $N_AND_W; }
    elsif ($align_type eq "s_and_w") {
        $align_type = $S_AND_W; }

    my $matname = $MATRIX_DIR . '/' . $MATRIX_FILE;
    if ( !($sub_mat = sub_mat_read ($matname))) {
        print STDERR "Fail reading substitution matrix from $matname\n";
        $err++;
    }

    if ($err) {
        undef @coord;
        return (EXIT_FAILURE);
    } else {
        return EXIT_SUCCESS;
    }
}

# ----------------------- mat_drivvel  ------------------------------
sub mat_drivvel ($)
{
    my $mat = shift;
    my ($min, $max, $av, $std_dev);
    score_mat_info ($mat, $min, $max, $av, $std_dev);
    print "Min: $min Max: $max Av: $av Std_dev $std_dev\n";
}

# ----------------------- store_mat    ------------------------------
sub store_mat ($ $ $ $)
{
    my ($n, $sec_info, $c1, $c2) = @_;
    my $seq1 = coord_get_seq ($c1);
    my $seq2 = coord_get_seq ($c2);

    my $fx_scr_mat  = score_mat_new (seq_size ($seq1), seq_size($seq2));
    score_fx   ($fx_scr_mat, $seq1, $c2, $fx_params);
    $fx_scr_mats[$n]  = $fx_scr_mat;

    if ($sec_s_scale != 0.0) {
        my $sec_scr_mat = score_mat_new (seq_size ($seq1), seq_size($seq2));
        score_sec  ($sec_scr_mat, $sec_info, $c2); 
        $sec_scr_mats[$n] = $sec_scr_mat;
    }
}

# ----------------------- get_model    ------------------------------
sub get_model ($ $ $ $)
{
    my ($n, $sec_info, $c1, $c2) = @_;
    if (0 == 1) {
        print "get_model working on ", coord_name ($c1),
        ", size ", coord_size ($c1); print " and ", coord_name ($c2),
        ", size ", coord_size ($c2); print "\n";
    }
    my ($seq1, $seq2);

    $seq1 = coord_get_seq ($c1);
    $seq2 = coord_get_seq ($c2);
    my $sim_scr_mat = score_mat_new (seq_size ($seq1), seq_size($seq2));
    my $total_mat   = score_mat_new (seq_size ($seq1), seq_size($seq2));
    my $tmp_sub_mat = sub_mat_shift ($sub_mat, $sim_mat_bottom);
    score_smat ($sim_scr_mat, $seq1, $seq2, $tmp_sub_mat);
    my $fx_scr_mat  = $fx_scr_mats[$n];
    my $fx_tmp_mat  = score_mat_shift ($fx_scr_mat, $fx_mat_shift);
    $total_mat    = score_mat_add ($fx_tmp_mat, $sim_scr_mat, $sim_scale);
    if ($sec_s_scale != 0.0) {
        my $sec_scr_mat = $sec_scr_mats[$n];
        $total_mat  = score_mat_add ($total_mat, $sec_scr_mat, $sec_s_scale);
    }

    my $pair_set =
        score_mat_sum_sec (my $result_mat, $total_mat,
                           $c2, $sec_pnlty,
                           $pgap_open, $pgap_widen,
                           $qgap_open, $qgap_widen,
                           $align_type);
    my $new_model = make_model ($pair_set, $seq1, $c2);
    return $new_model;
}

# ----------------------- activate     ------------------------------
# This is modelled after the activation (logistic) function
# used in neural nets.
sub activate ($)
{
    my $x = shift;
    my $offset = $SWITCH_POINT;
    my $xscl = $X_SCALE;
    my $a = $offset - $x;
    $a *= $xscl;
    $a = exp ($a);
    return ( 1.0 / ( 1.0 + $a));
}


# ----------------------- pair_err     ------------------------------
# For one pair of structures, get the error.
sub pair_err ( $ $ $ $ )
{
    my ($n, $sec_info, $c1, $c2) = @_;
    my $tmp_c = get_model ( $n, $sec_info, $c1, $c2);
    if ( !$tmp_c) {
        return 0; }
    my $smallsize = coord_size ($c1);
    {
        my $size2 = coord_size ($c2);
        if ($size2 < $smallsize) {
            $smallsize = $size2; }
    }
    my $model_size = coord_size ($tmp_c);
    if (( $model_size / $smallsize ) < $MIN_MODEL) {
        return 0;
    }


    my $fraction = 0.0;
    if ( ! dme_thresh ( $fraction, $tmp_c, $c1, 4.0)) {
        print STDERR "dme error comparing coordinates from..\nModel ",
        coord_name ($tmp_c), "\n";
        print " of size ", coord_size ($tmp_c), " and \n";
        print coord_name ($c1), " of size ", coord_size ($c1), "\n";
        return 1;
    }

    my $x = - activate ($fraction);
    return ($x);
}

# ----------------------- print_two    ------------------------------
# Print names and stuff about a pair of proteins.
# This is very verbose, but fun to have a look at.
sub print_two ($ $ $)
{
    my ($c1, $c2, $err) = @_;
    print coord_name ($c1), ', ', coord_name ($c2),
    ",  error $err, sizes @{[coord_size($c1)]} and @{[coord_size($c2)]}\n";
}

# ----------------------- align_cost   ------------------------------
sub align_cost (@)
{
    my @values = @_;
    my $num_pair = ($#coord + 1) / 2;

    {
        no strict 'refs';
        use strict 'vars';
        for (my $i = 0; $i <= $#param_name; $i++) {
            $ {$param_name[$i] } = $values[$i];
        }
    }

    use strict 'refs';
    my $err = 0;
    if ($pgap_widen >= $pgap_open) {
      print STDERR
      "Step warn because pgap_widen >= pgap_open,($qgap_widen,$qgap_open)\n";
    }
    if ($qgap_widen >= $qgap_open) {
      print STDERR
      "Step warn because qgap_widen >= qgap_open,($qgap_widen,$qgap_open)\n";
#      return $err;
    }



    my @err_values;
    $err_values[$num_pair - 1] = 0;
    if ($first) {
        $first = 0;
        for (my $i = 0; $i < $num_pair; $i++) {
            my $a = $i * 2;
            my $b = $a + 1;
            store_mat ($a, $sec_info[$a], $coord[$a], $coord[$b]);
            store_mat ($b, $sec_info[$b], $coord[$b], $coord[$a]);
        }
#       Now, we can free some memory only used in initial scoring       
        if ( defined ($fx_params)) {
            print "Getting rid of fx_params\n";
            undef $fx_params;
        }
        undef (@sec_info);
    }
    for (my $i = 0; $i < $num_pair; $i++) {
        my $a = $i * 2;
        my $b = $a + 1;
        my $x = pair_err ($a, $sec_info[$a], $coord[$a], $coord[$b]);
        $err += $x;
        $err_values [$a] = $x;
        $x    = pair_err ($b, $sec_info[$b], $coord[$b], $coord[$a]);
        $err += $x;
        $err_values [$b] = $x;
    }

    if ($count++ % 50 == 0) {
        my ($best, $worst) = 0;
        for ( my $i = 1; $i < $num_pair; $i++) {
            if ($err_values[$i] < $err_values[$best]) {
                $best = $i; }
            if ($err_values[$i] > $err_values[$worst]) {
                $worst = $i; } }
        my $a = $best * 2;
        my $b = $a + 1;
        print "For fun at f eval $count\n";
        print "Best   proteins ";
        print_two ($coord[$a], $coord[$b], $err_values[$best]);
        $a = $worst * 2;
        $b = $a + 1;
        print "Worst  proteins ";
        print_two ($coord[$a], $coord[$b], $err_values[$worst]);
    }
    return ( $err);
}

return 1;
