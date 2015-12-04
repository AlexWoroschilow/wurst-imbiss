# 3 Nov 2005
# A new script for optimising penalties associated with the 
# pure sequence based fragments.
# 19 Dec 2003

# Most of the action happens in Wutil.pm
# Plan..
# At the start, read up every protein and store the
# sequence-fragment profile.
# rcsid = $Id: seqfrag_sfrag.pl,v 1.1 2007/09/28 16:57:00 mmundry Exp $

=pod

=head1 NAME

something.pl

=head1 SYNOPSIS

 something.pl -S<[ B<-m> I<matrix>]

=options

=over

=item B<-v >I<verbosity>

Set verbosity level. This is not used much at the moment, but is
there for future hooks.

=back


=cut

use strict;
undef $@;
use warnings;
use FindBin;
use lib "$ENV{HOME}/../torda/pl/lib";  # Where wurst lives after installation
#use lib "$FindBin::Bin/../src/Wurst/blib/arch";
#use lib "$FindBin::Bin/../src/Wurst/blib/lib";

use Wurst;
use lib "$FindBin::Bin../blib/arch";
use lib "$FindBin::Bin../blib/lib";
#use forks;
#use forks::shared;

use POSIX qw(EXIT_SUCCESS EXIT_FAILURE);

# ----------------------- get_pair_set_mat_precalc ------------------
# This is for alignments based on pure sequence fragments, but the
# framework should function for any time where we are working with
# pre-calculated score matrices.
# 
sub get_pair_set_mat_precalc ( $ $ \%)
{
    my ($c1, $c2, $aln_param) = @_;

    my $up   = 1.2;   # Errors / range to be applied to gap opening
    my $down = 0.8;

    my $index = coord_name($c1) . coord_name ($c2);

    my $fx_mat = $$aln_param{f_smat}{$index};

    my $sim_scr_mat = score_mat_new (coord_size ($c1), coord_size($c2));

    my $total_mat = score_mat_shift ($fx_mat, $$aln_param{m_shift});

    my ($po, $pw, $qo, $qw);
    $po = $$aln_param {pgap_open};
    $pw = $$aln_param {pgap_widen};
    $qo = $po;
    $qw = $pw;

    my $pair_set1 =
        score_mat_sum_smpl (my $result_mat, $total_mat,
                            $po, $pw,
                            $qo, $qw,
                            $$aln_param{align_type});
    my $pair_set2 =
        score_mat_sum_smpl (my $crap_mat, $total_mat,
                            $po * $up, $pw,
                            $qo * $up, $qw,
                            $$aln_param{align_type});
    my $pair_set3 =
        score_mat_sum_smpl (my $crap2_mat, $total_mat,
                            $po * $down, $pw,
                            $qo * $down, $pw,
                            $$aln_param{align_type});

    return ($pair_set1, $pair_set2, $pair_set3);
}


# ----------------------- inner_cost --------------------------------
# This is the inner loop, turned into a function for threading
sub inner_cost ( $ $ $ $ $)
{
    my ($coords, $sub_mat, $aln_param, $first, $last) = @_;
    my $err = 0.0;
    my $profs = $$aln_param{profs};
    my @p_s;
    for (my $i = $first; $i < $last; $i++) {
        my $c = $$coords[$i];
        my $tmp_err = 0;
        my $c1 = $ {$c}[0];
        my $c2 = $ {$c}[1];
        @p_s = get_pair_set_mat_precalc ($c1, $c2, %$aln_param );
        if (! @p_s) { die "broke 1\n"};
        if ( @p_s == 0) { die "broke 2\n"};
        foreach my $p_s (@p_s) {
            $err += pair_err ($p_s, $c1, $c2); }

        @p_s = get_pair_set_mat_precalc ( $c2, $c1, %$aln_param );
        foreach my $p_s (@p_s) {
            $err += pair_err ($p_s, $c2, $c1); }
    }

    $err /= (@p_s * 2);       # normalise for number of alignments

    return ($err);
}

# ----------------------- align_cost --------------------------------
sub align_cost (\@ \% $)
{
    use FindBin;
    use lib "$FindBin::Bin";
    use Pair_err;
    use Wutil qw ( get_pair_set_str_seq_prof_prof array2mat zero_shift_mat);
    my $var_param = shift;
    my $fix_param = shift;
    my $names     = shift;  # ref to array of parameter names
    use vars qw($zero_shift);
    *zero_shift = \undef();
#   These are fixed parameters
    my ($coords, $sub_mat, $profs, $f_smat, $sub_mat_vary, $class_file);
    my ($align_type);
#   These guys can vary
    my ($pgap_open, $pgap_widen,
        $qgap_open, $qgap_widen, $m_shift);
    undef $@;

#   We pick variables out of the hash. First, we have to make sure the
#   variable exists. We do a harmless assignment of "1" to the
#   variable. If this does not generate an error, we then assign the real
#   value. This convolution is necessary to let us assign "undef" to
#   variables.
    for (my $i = 0; $i < @$names; $i++) {
        my $var = $$names[$i];
        my $code = "\$$var = \$\$var_param[$i]";
        my $r = eval $code;
        if (( ! defined ($r)) || ($@)) {
            die "Cost function failed on variable parameter $var\n$@\n";}
    }


    for my $var (keys %$fix_param) {
        my $code = "\$$var = 1";
        my $r = eval $code;
        if (( ! defined ($r)) || ($@)) {
            die "Cost function failed on fixed parameter $var\n$@\n"; }
        $code = "\$$var = \$\$fix_param { $var }";
        $r = eval $code;
    }

#   Now, we get ready for an alignment.
#   Let's build a hash, %aln_param which will contain everything we
#   need for an alignment. It's members, however, may come from
#   %fix_param or %var_param.

    my $aln_param;
    $aln_param = {
        'm_shift'    =>  $m_shift,
        'pgap_open'  =>  $pgap_open,
        'pgap_widen' =>  $pgap_widen,
        'f_smat'     =>  $f_smat
    };

#   Take our array, and put it into the substitution matrix, if we are
#   varying the substitution matrix.
    if ( ! exists ($$fix_param {sub_mat_vary})) {
        die "Program bug. In this version, sub_mat_vary vary must be defined."}
    if ( $$fix_param {sub_mat_vary} ) {
        my $start = $#$names + 1;
        my $last =  $#$var_param;
        my @t = @$var_param[$start .. $last];
        array2mat ( $sub_mat, \@t);
    }

#   Now, we may have to do a zero shift of the matrix, and
#   don't forget to put it back before returning

    if (defined ($zero_shift)) {
        zero_shift_mat ($sub_mat, $zero_shift);}

    my $err = 0.0;
    my @p_s;
    my $pchunk = int (@$coords / 2) + 1;

    my $np = 2; # num threads
    my $first = 0;

    while ($first < @$coords) {
        for (my $i=0; $i < $np && $first < @$coords; $i++, $first += $pchunk) {
            my $last = $first + $pchunk;
            if ($last > @$coords) {$last = @$coords; }
            $err +=inner_cost ($coords, $sub_mat, $aln_param, $first, $last); 
        }
    }
#     while ($first < @$coords) {
#         my @thr;
#         for (my $i=0; $i < $np && $first < @$coords; $i++, $first += $pchunk) {
#             my $last = $first + $pchunk;
#             if ($last > @$coords) {$last = @$coords; }
#             $thr [$i] = threads->new
#                 (\&inner_cost, $coords, $sub_mat, $aln_param, $first, $last);
#             if ( ! $thr[$i]) {
#                 die "Create thread failure on $i, $!\n";}
#         }
#         for ( my $i = 0 ; $i < @thr; $i++) {
#             $err += $thr[$i]->join; }
#     }

    if (defined ($zero_shift)) {
        zero_shift_mat ($sub_mat, -$zero_shift);}

    $err /= @$coords;   # and the number of proteins in set
    return $err;
}


# ----------------------- init_smat_sfrag_sfrag ---------------------

sub init_smat_sfrag_sfrag ( \% \@ \%)
{
    my ($f_smat, $coords, $fixed) = @_;
    my $class_file = $$fixed{class_file};
    my $aa_clssfcn;
    if ( ! (-f $class_file)) {
        print STDERR "\"$class_file\" is not a normal file. ";
        print STDERR "Disaster approaching\n";
    }
    if ( ! ($aa_clssfcn = ac_read ( $class_file ))) {
        print STDERR "Failed reading classification\n";
        return undef;
    }
    for my $c (@$coords) {
        my ($pvec1, $pvec2);
        my $c1 = $ {$c}[0];
        my $c2 = $ {$c}[1];
        my $name1 = coord_name ($c1);
        my $name2 = coord_name ($c2);
        if (! ($pvec1 = seq_2_prob_vec(coord_get_seq($c1), $aa_clssfcn))) {
            warn "Fail prob vec on $name1\n", return undef;}
        if (! ($pvec2 = seq_2_prob_vec(coord_get_seq($c2), $aa_clssfcn))) {
            warn "Fail prob vec on $name2\n", return undef;}
        my $size1 = coord_size ($c1);
        my $size2 = coord_size ($c2);
        my $mat1 = score_mat_new ($size1, $size2);
        my $mat2 = score_mat_new ($size2, $size1);
        score_pvec ($mat1, $pvec1, $pvec2);
        score_pvec ($mat2, $pvec2, $pvec1);
        $$f_smat {"$name1$name2"} = $mat1;
        $$f_smat {"$name2$name1"} = $mat2;
    }
    return 1;
}

# ----------------------- mymain    ---------------------------------
sub mymain ()
{
    use Getopt::Std;
    use Wutil qw ( array2mat
                  init_main
                  kill_handlers opt_simplex split_coord);

    wurstli_hello();

    my $verbosity = 1;
    my %opts;

    my $seed = undef;

    if ( ! getopts ('r:v:', \%opts)) {
        usage(); }
    if (defined ($opts {v})) {
        $verbosity = $opts {v}; }
    if (defined ($opts {r})) {
        $seed = $opts {r}; }

    undef %opts;
    if (defined ($seed)) {
        srand ($seed) }
    else {
        srand (2593); }

#   Anything which is going to persist gets declared up here
    my (%s_args, %fixed, @coords, @test_coords, %test_fix_param);
    my %clssfcns;                   # Where we will read up profiles
    my %f_smat;

    $fixed {f_smat}           = \%f_smat;
    $test_fix_param{f_smat}   = \%f_smat;


    $fixed {coords} = \@coords;
    $test_fix_param {coords} = \@test_coords;
    if ( ! init_main( %s_args, %fixed, %test_fix_param)) {
        warn "Failed to initialise\n"; return EXIT_FAILURE;}


    if ( ! (init_smat_sfrag_sfrag (%f_smat, @coords, %fixed))) {
        warn "Failed setting up score mats\n"; return EXIT_FAILURE; }
    if ( ! (init_smat_sfrag_sfrag (%f_smat, @test_coords, %fixed))) {
        warn "Failed setting up score mats\n"; return EXIT_FAILURE; }

    keys(%f_smat) = (@coords + @test_coords) * 2;

    if ($fixed {sub_mat_vary}) {        # Copy the sub matrix to simplex params
        mat2array ($s_args {ini_pt}, $fixed {sub_mat}); }
    kill_handlers();
    {
        my $c = $#{$fixed{coords}} + 1;
        my $t = $#{$test_fix_param {coords}} + 1;
                   
        print "My coordinate list has ", $c + $t, " pairs\n",
        "with $c for parameterising and $t reserved for testing\n";
    }
    my $r = opt_simplex(%s_args, %fixed);
    if ($fixed {sub_mat_vary}) {        # Copy optimised values back to matrix
        array2mat ($fixed {sub_mat}, $$r{best});
        print "Final matrix...\n", sub_mat_string ($fixed {sub_mat});
    }
    use Simplex2;
    if (! defined $$r{success}) {
        print STDERR "Seriously broken in simplex (result undef)\n";
        return EXIT_FAILURE;
    }
    if ($$r{success} == $SPLX_BROKEN) {
        return EXIT_FAILURE; }
    if ($$r{success} == $SPLX_TOO_MANY) {
        print "Run did not converge\n";
    } elsif ($$r{success} == $SPLX_SUCCESS) {
        print "Run claims to have converged\n";
    } else {
        print "Unknown return code from simplex\n"; }
    use Sys::Hostname;
    my $host = hostname();
    print "Run on $host\n";

    return EXIT_SUCCESS;
}

# ----------------------- main      ---------------------------------
exit (mymain());
