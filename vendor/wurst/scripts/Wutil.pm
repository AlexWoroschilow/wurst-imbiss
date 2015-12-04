# 21 July 2003
# Certain functions keep cropping up in the various scripts.
# Here, we gather some of the common things so as to avoid some
# of the code duplication.

package Wutil;
#rcsid = "$Id:"
=pod

=head1 NAME

Wutil.pm

=head1 SYNOPSIS

 use wutil;

=head1 DESCRIPTION

Various common functions for our wurst routines.  Generally, we
operate by passing a reference to arrays or whatever.

For function arguments, follow the convention of

  func (destination, source)

and, where sensible, return a reference to the destination. This
allows the caller to write slightly more concise code.

Also, we ask the caller to declare arrays, rather than rely on
perl growing and giving birth to them. It may not be important,
but the philosophy is that usually the caller allocates space. In
perl, this may not mean too much.

=head2 array2mat

  $sub_mat = array2mat ($sub_mat, $array_ref).

Given an array of 210 elements, copy it into a substitution
matrix. Both array and substitution matrix are provided.

=head2 catch_kill

If we are sent a signal, call exit and die.

=head2 check_pdb_files

 if (check_pdb_files(@pairlist, @bin_dir_list, @coords) == EXIT_FAILURE) {
     die;}

Take a two-dimensional array of protein files in I<@pairlist>, a
set of directories in I<bin_dir_list> and an (empty) array in
I<@coords>. Fill the array with coordinate structures.

=head2 get_pair_set

For use inside the score function.
Given a set of parameters and things to align, return a pair set.


=head2 get_param_list

 if (  get_param_list ($param_file, \%s_args, \%fixed) == EXIT_FAILURE ) {
     die;}

Read parameters for a minimisation run.

=head2 mat2array

 my @array;
 $array_ref = mat2array (\@array, $sub_mat);

Take a substitution matrix and append the elements to the array.
B<Note>: We append the elements to the S<array !> If the array has
two elements to start, it will be returned with S<2 + 210>
elements. This lets us put some special numbers at the start of
an array (like named parameters), then append the bulk to it.

If you use this, remember to pop magic elements from the array
before turning it back into a matrix.

=head2 split_coord

    split_coord (@coord_array, @test_array);

Take our array of coord structures and move every N'th member to
the I<test_array>.

=head2 get_pair_list

 my @pairlist = get_pair_list ($pair_file);

Get a list of protein pairs from I<$pair_file>. Return a
two-dimensional array.

=head2 init_main

Even the initialisation seems to be common between the various
optimisation runs. Now, it is all going into this function.

=head2 zero_shift_mat

 zero_shift_mat ($sub_mat, $zero_shift);

Apply the value in $zero_shift to the substitution matrix given by $sub_mat.

=cut

use strict;
use POSIX qw (EXIT_SUCCESS EXIT_FAILURE);
use vars qw (@EXPORT_OK @ISA);
@ISA = qw (Exporter);
@EXPORT_OK = qw (add_diag array2mat bounded_mover check_pdb_files get_pair_list
                 get_pair_list
                 get_pair_set
                 get_pair_set_prof_seq
                 get_pair_set_struct
                 get_pair_set_str_seq
                 get_pair_set_str_seq_prof
                 get_pair_set_prof_seq_prof_seq
                 get_pair_set_str_seq_prof_prof
                 get_param_list init_main init_profiles
                 init_smat_fx_prof
                 kill_handlers
                 mat2array opt_simplex rand_int rand_int_non_zero
                 split_coord
                 use_diag zero_shift_mat );


use Wurst;
{
# ----------------------- almost global -----------------------------
# This array must be visible to two functions
    my @letters = ( 'g', 'a', 'v', 'l', 'i', 'f', 'p', 's', 't', 'c',
                    'm', 'w', 'y', 'n', 'q', 'd', 'e', 'k', 'r', 'h');

# ----------------------- get_blast_prof ----------------------------
# Little helper function to return the blast profile.
sub  get_blast_prof ($)
{
    my $cname = shift;
    use FindBin;
    use lib "$FindBin::Bin";
    use Paths qw ($PROFILE_DIR);
    my $n1 = $cname;
    substr ($n1, 0, 4) = lc (substr ($n1, 0, 4));
    $n1 = "$PROFILE_DIR/${n1}.prof";
    my $prof = blst_chk_read ($n1);
    if ( ! $prof) {
        warn "Fail reading blast profile for $n1\n"; return undef; }
    return $prof;
}

# ----------------------- check_prof_and_seq ------------------------
# Check that a sequence matches a profile. We could check residue by
# residue. Let us at least check the size. This gets its own function
# since it should be rather verbose.
sub check_prof_and_seq ( $ $ $)
{
    my ($seq, $prof, $name) = @_;
    my $size1 = seq_size ($seq);
    my $size2 = seq_size ( seqprof_get_seq ($prof));
    if ($size1 != $size2) {
        print STDERR "Working on $name. Size from profile $size2. ",
                     "size from sequence $size1\n";
        print STDERR "From sequence...\n", seq_print($seq), "\n";
        print STDERR "From profile...\n",
                      seq_print(seqprof_get_seq ($prof)), "\n";
        return undef;
    }
    return 1;
}


# ----------------------- init_smat_fx_prof -------------------------
# Given a list of coordinate pairs, calculate the score matrix.
# Put each into hash, indexed by the protein names.
# WARNING: This version uses a sequence profile ! Not a sequence.
# It is best to use a hash to store results. A simple array was used
# in the past, but the hash lets us not worry about multiple sets
# of coordinates or preserving the order. This function can also be
# called more than once, simply adding to the hash elements.
sub init_smat_fx_prof ( \% \@)
{
    my ($f_smat, $coords) = @_;
    use FindBin;
    use lib "$FindBin::Bin";
    use Paths qw ($FX_PARAM_DIR $FX9_PARAM_FILE);
    my $fx_params;
    if ( !defined $FX_PARAM_DIR) {
        warn "FX_PARAM_DIR not defined\n"; return undef; }
    if ( ! defined ($FX9_PARAM_FILE)) {
        warn "FX9_PARAM_FILE not defined\n"; return undef; }
    my $tmp_p = "$FX_PARAM_DIR/$FX9_PARAM_FILE";
    if ( ! (-f $tmp_p)) {
        warn "Cannot find fx param file $tmp_p\n"; return undef; }
    if ( !($fx_params = param_fx_read ($tmp_p))) {
        warn "Fail reading parameters from $tmp_p\n"; return undef; }

    foreach my $c (@$coords) {
        my ($prof1, $prof2);
        my $c1 = $ {$c}[0];
        my $c2 = $ {$c}[1];
        my $seq1 = coord_get_seq ($c1);
        my $seq2 = coord_get_seq ($c2);
        my $name1 = coord_name ($c1);
        my $name2 = coord_name ($c2);
        my $fx_mat1 = score_mat_new (seq_size($seq1), seq_size ($seq2));
        my $fx_mat2 = score_mat_new (seq_size($seq2), seq_size ($seq1));
        if (! ( $prof1 = get_blast_prof ($name1))) {return undef;}
        if (! ( $prof2 = get_blast_prof ($name2))) {return undef;}
        if (! ( check_prof_and_seq ($seq1, $prof1, $name1))) {
            return undef; }
        if (! ( check_prof_and_seq ($seq2, $prof2, $name2))) {
            return undef; }
        undef ($seq2); undef ($seq1);
        score_fx_prof ($fx_mat1, $prof1, $c2, $fx_params);
        score_fx_prof ($fx_mat2, $prof2, $c1, $fx_params);
        $$f_smat {"$name1$name2"} = $fx_mat1;
        $$f_smat {"$name2$name1"} = $fx_mat2;
    }
    undef $fx_params;                        # Not necessary, but let's force it
    return 1;
}

# ----------------------- array2mat     -----------------------------
# pack a long array of parameters into a substitution matrix.

    sub array2mat ($ $)
    {
        my ($sub_mat, $ary) = @_;
        my $name = 'array2mat';
        if ((my $t = $#$ary) != 209) {
            warn "$name: Expected array of 210 elements, got $t\n";
            return undef;
        }
        my $k = 0;
        for (my $i = 0; $i <= $#letters; $i++) {
            for (my $j = $i; $j <= $#letters; $j++) {
                my $a = $letters[$i];
                my $b = $letters[$j];
                sub_mat_set_by_c ($sub_mat, $a, $b, $$ary[$k]);
                sub_mat_set_by_c ($sub_mat, $b, $a, $$ary[$k++]);
            }
        }
        return $sub_mat;
    }

# ----------------------- mat2array     -----------------------------
# We are given a reference to an array. Fill it with the elements
# from the substitution matrix.
    sub mat2array ($ $)
    {
        my ($ary, $sub_mat) = @_;
        my $myname = 'mat2array';
        if (! $ary) {
            warn "$myname: Please give me a reference to an  existing array\n";
            return undef;
        }
        if (! $sub_mat) {
            warn "$myname: Please give me a valid substitution matrix\n";
            return undef;
        }
        for (my $i = 0; $i <= $#letters; $i++) {
            for (my $j = $i; $j <= $#letters; $j++) {
                my $a = $letters[$i];
                my $b = $letters[$j];
                my $ini = sub_mat_get_by_c ($sub_mat, $a, $b);
                push (@$ary, $ini);
            }
        }

        return $ary;
    }
}

# ----------------------- check_pdb_files ---------------------------
# Have a look for binary protein files.
# For this version, I am going to use a two dimensional array
# since it makes other parts of the code prettier.
# First arg is a set of protein names.
# Second is the list of directories.
# Return undef on failure
sub check_pdb_files (\@ \@ \@)
{
    my $pairlist  = shift;
    my $bin_dirs  = shift;
    my $coords    = shift;
    my $ret = 1;

    my $i = 0;
    PROTBIN: foreach my $p (@$pairlist) {
        my $found = 0;
        foreach my $d (@$bin_dirs) {
            my $path = "$d/$p.bin";
            if ( -e $path) {
                $found = 1;
                if (! (my $c = coord_read ($path))) {
                    print STDERR "coord_read error on $path\n";
                    $ret = undef;
                } else {
                    $$coords [ int ($i / 2) ] [ $i % 2] = $c;
                    $i++;
                    next PROTBIN;
                }
            }
        }
        if ( ! $found) {
            print STDERR "Could not find $p\n";
            $ret = undef;
        }
    }
    if (!$ret) {
        print STDERR "Directory list for files was @$bin_dirs\n";}
    return ($ret);
}

# ----------------------- split_coord -------------------------------
# Given an array ($coord), remove every n'th member and put it in
# the second array.  We do not care about rounding errors. We just
# want to remove some selection of element.
sub split_coord (\@ \@)
{
    my $coords = shift;
    my $test_coords = shift;
    use Al_defaults qw($TEST_FRAC);
    if ($TEST_FRAC < 1e-6) {
        $test_coords = undef; return ; }
    my $num_save = int (($#$coords + 1) * $TEST_FRAC);
    my $freq = int (1.0 / $TEST_FRAC) - 1;
    for (my $i = $freq; $i < $#$coords; $i += $freq) {
        my $t = splice (@$coords, $i, 1);
        push (@$test_coords, $t);
    }
}

# ----------------------- get_pair_list -----------------------------
sub get_pair_list ($)
{
    use Get_line;
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
    return EXIT_SUCCESS;
}

# ----------------------- get_param_list ----------------------------
sub get_param_list ($ $ $)
{
    my ($filename, $s_args, $fixed) = @_;

    my $result = EXIT_SUCCESS;
    if (! open (PARAMS, "<$filename" )) {
        print STDERR "Open fail on $filename: $!\n";
        print STDERR "Working dir is ", `pwd`;
        return EXIT_FAILURE }

    while (defined (my $line = get_line (\*PARAMS))) {
        my ($name, $val) = split '[\s,]+', $line, 2;
#       my @val = split ' *, *', $val; # old, space tolerant version
        my @val = split '[\s,]+', $val;
        $name =~ tr/A-Z/a-z/;
      SWITCH: {
          if ($name eq 'dipep_file')  {$$s_args{$name} = $val; last SWITCH; }
          if ($name eq 'f_tol')       {$$s_args{$name} = $val; last SWITCH; }
          if ($name eq 'fixed')     {$$fixed{$val[0]} = $val[1]; last SWITCH; }
          if ($name eq 'ini_range') {$$s_args{ini_range} = \@val; last SWITCH;}
          if ($name eq 'lower')       {$$s_args{lower}  = \@val; last SWITCH; }
          if ($name eq 'max_iter')    {$$s_args{$name} = $val; last SWITCH; }
          if ($name eq 'max_restart') {$$s_args{$name} = $val; last SWITCH; }
          if ($name eq 'names')       {$$s_args{names}  = \@val; last SWITCH; }
          if ($name eq 'num_step')    {$$s_args{$name}  = $val; last SWITCH; }
          if ($name eq 'o_file')      {$$s_args{$name} = $val; last SWITCH; }
          if ($name eq 'start')       {$$s_args{ini_pt} = \@val; last SWITCH; }
          if ($name eq 'sub_mat')     {$$s_args{$name} = $val; last SWITCH; }
          if ($name eq 't_ini')       {$$s_args{$name} = $val; last SWITCH; }
          if ($name eq 't_final')     {$$s_args{$name} = $val; last SWITCH; }
          if ($name eq 'upper')       {$$s_args{upper}  = \@val; last SWITCH; }
          if ($name eq 'scatter')     {$$s_args{$name} = $val; last SWITCH; }
          if ($name eq 'write_all')   {$$s_args{$name} = $val; last SWITCH; }
          {print STDERR "Param file unknown line \"$line\"\n"; }
      }
    }

    $result = parameters_sane ($s_args);  # Check if input params are sensible
    return $result;
}

# ----------------------- get_pair_set  -----------------------------
# We are given two structures, a substitution matrix and
# (later), even more stuff like secondary structure or
# information for sequence to structure alignments. Here, we
# all the stuff to get an alignment.  Maybe, we will have to
# retrieve a pre-calculated score matrix (in the case of
# expensive score functions.
# Now, the tricky, ad hoc bit.
# Lets calculate three pair sets, based on increasing and
# decreasing the gap opening penalty.
sub get_pair_set ( $ $ $ \%)
{
    my ($c1, $c2, $sub_mat, $aln_param) = @_;
    my ($rmat);
    my $up   = 1.2;
    my $down = 0.8;
    my $seq1 = coord_get_seq ($c1);
    my $seq2 = coord_get_seq ($c2);
    my $sim_scr_mat = score_mat_new (seq_size ($seq1), seq_size($seq2));

    score_smat ($sim_scr_mat, $seq1, $seq2, $sub_mat);
    $rmat = undef;
    my ($pair_set1, $pair_set2, $pair_set3);
    my ($po, $pw, $qo, $qw);
    $po = $$aln_param {'pgap_open'};
    $pw = $$aln_param {'pgap_widen'};
    $qo = $$aln_param {'qgap_open'};
    $qw = $$aln_param {'qgap_widen'};

    $pair_set1 = score_mat_sum_smpl ($rmat, $sim_scr_mat,
                                     $po, $pw,
                                     $qo, $qw,
                                     $$aln_param {al_type});

    $pair_set2 = score_mat_sum_smpl ($rmat, $sim_scr_mat,
                                     $po * $up, $pw,
                                     $qo * $up, $qw,
                                     $$aln_param {al_type});

    $pair_set3 = score_mat_sum_smpl ($rmat, $sim_scr_mat,
                                     $po * $down, $pw,
                                     $qo * $down, $qw,
                                     $$aln_param {al_type});


    return ($pair_set1, $pair_set2, $pair_set3);
}

# ----------------------- get_pair_prof_seq -------------------------
# This is like get_pair_set, but for profiles against sequences.
# For the first runs, we will just work with one value for the gap
# penalties. When we are a bit closer to sensible values, go
# back to the version with up and down variation.
sub get_pair_set_prof_seq ( $ $ $ \%)
{
    my ($prof, $seq, $sub_mat, $aln_param) = @_;
    my ($rmat);
    my $up   = 1.2;
    my $down = 0.8;

    my $sim_scr_mat =
        score_mat_new (seq_size (seqprof_get_seq ($prof)), seq_size($seq));

    if ( !score_sprof ($sim_scr_mat, $prof, $seq, $sub_mat)) {
        return undef; }
    $rmat = undef;
    my ($pair_set1, $pair_set2, $pair_set3);
    my ($po, $pw, $qo, $qw);
    $po = $$aln_param {'pgap_open'};
    $pw = $$aln_param {'pgap_widen'};
    $qo = $$aln_param {'qgap_open'};
    $qw = $$aln_param {'qgap_widen'};

    $pair_set1 = score_mat_sum_smpl ($rmat, $sim_scr_mat,
                                     $po, $pw,
                                     $qo, $qw,
                                     $$aln_param {al_type});

    $pair_set2 = score_mat_sum_smpl ($rmat, $sim_scr_mat,
                                     $po * $up, $pw,
                                     $qo * $up, $qw,
                                     $$aln_param {al_type});

    $pair_set3 = score_mat_sum_smpl ($rmat, $sim_scr_mat,
                                     $po * $down, $pw,
                                     $qo * $down, $qw,
                                     $$aln_param {al_type});


    return ($pair_set1, $pair_set2, $pair_set3);
}

# ----------------------- get_pair_set_prof_seq_prof_seq ------------
# This does profile - profile alignments without any structural
# information
sub  get_pair_set_prof_seq_prof_seq ( $ $ $ \%)
{
    my ($prof1, $prof2, $sub_mat, $aln_param) = @_;
    my ($rmat);
    my $up   = 1.2;
    my $down = 0.8;

    my $sim_scr_mat =
        score_mat_new (seq_size (seqprof_get_seq ($prof1)),
                       seq_size (seqprof_get_seq ($prof2)));

    if ( !score_prof_prof ($sim_scr_mat, $prof1, $prof2, $sub_mat)) {
        return undef; }
    $rmat = undef;
    my ($pair_set1, $pair_set2, $pair_set3);
    my ($po, $pw, $qo, $qw);
    $po = $$aln_param {'pgap_open'};
    $pw = $$aln_param {'pgap_widen'};
    $qo = $$aln_param {'qgap_open'};
    $qw = $$aln_param {'qgap_widen'};

    $pair_set1 = score_mat_sum_smpl ($rmat, $sim_scr_mat,
                                     $po, $pw,
                                     $qo, $qw,
                                     $$aln_param {al_type});

    $pair_set2 = score_mat_sum_smpl ($rmat, $sim_scr_mat,
                                     $po * $up, $pw,
                                     $qo * $up, $qw,
                                     $$aln_param {al_type});

    $pair_set3 = score_mat_sum_smpl ($rmat, $sim_scr_mat,
                                     $po * $down, $pw,
                                     $qo * $down, $qw,
                                     $$aln_param {al_type});

    return ($pair_set1, $pair_set2, $pair_set3);
}




# ----------------------- get_pair_set_str_seq-----------------------
# This does an alignment based on a the fx9 function and a substitution
# matrix. The fx9 function matrix is precalculated and retrieved from
# the f_smat hash.
sub get_pair_set_str_seq ( $ $ \%)
{
    my ($c1, $c2, $aln_param) = @_;
    my $seq1 = coord_get_seq ($c1);
    my $seq2 = coord_get_seq ($c2);

    my $up   = 1.2;
    my $down = 0.8;

    my $index = coord_name($c1) . coord_name ($c2);

    my $fx_mat = $$aln_param{f_smat}{$index};
    my $sim_scr_mat = score_mat_new (seq_size ($seq1), seq_size($seq2));
    my $total_mat   = score_mat_new (seq_size ($seq1), seq_size($seq2));
    score_smat ($sim_scr_mat, $seq1, $seq2, $$aln_param{sub_mat});

#   Now, I should apply a shift to the sim_scr_mat (next version)

    my $fx_wt = 1.0 - $$aln_param{seq_wt};
    my $tmp_fx_mat = score_mat_shift ($fx_mat, $$aln_param{fx_shift});
    $total_mat = score_mat_add ($total_mat, $tmp_fx_mat, $fx_wt);
    undef ($tmp_fx_mat);
    $total_mat = score_mat_add ($total_mat, $sim_scr_mat, $$aln_param{seq_wt});


    my ($po, $pw, $qo, $qw);
    $po = $$aln_param {pgap_open};
    $pw = $$aln_param {pgap_widen};
    $qo = $$aln_param {qgap_open};
    $qw = $$aln_param {qgap_widen};

    my $pair_set1 =
        score_mat_sum_smpl (my $result_mat, $total_mat,
                            $po, $pw,
                            $qo, $qw,
                            $$aln_param{align_type});
    my $pair_set2 =
        score_mat_sum_smpl (my $result_mat, $total_mat,
                            $po * $up, $pw,
                            $qo * $up, $qw,
                            $$aln_param{align_type});
    my $pair_set3 =
        score_mat_sum_smpl (my $result_mat, $total_mat,
                            $po * $down, $pw,
                            $qo * $down, $pw,
                            $$aln_param{align_type});

    return ($pair_set1, $pair_set2, $pair_set3);
}

# ----------------------- get_pair_set_str_seq_prof -----------------
# Do an alignment based on both sequence and structure, but using
# profiles in both cases.
# For the sequence to structure, there is a precalculated score
# matrix, so it is not seen here.
# This function began its life on 19 Dec.
sub get_pair_set_str_seq_prof ( $ $ \%)
{
    my ($c1, $c2, $aln_param) = @_;
    my $sub_mat = $$aln_param{sub_mat};
    my $profs = $$aln_param{profs};
    my $prof1 = $$profs{ coord_name ($c1)};
    my $seq2 = coord_get_seq ($c2);

    my $up   = 1.2;
    my $down = 0.8;

    my $index = coord_name($c1) . coord_name ($c2);

    my $fx_mat = $$aln_param{f_smat}{$index};

    if ( ! (defined ($prof1))) {
        my $names = 'profile not defined on ' . coord_name ($c1) .
            ' structure was ' . coord_name ($c2) . "\n";
        die $names;
    }
    if ( ! (defined ($fx_mat))) {
        my $names = coord_name ($c1) . ' ' . coord_name ($c2);
        die "fx_mat not defined on $names\n";
    }

    my $sim_scr_mat = score_mat_new (coord_size ($c1), coord_size($c2));
    my $total_mat   = score_mat_new (coord_size ($c1), coord_size($c2));
    score_sprof ($sim_scr_mat, $prof1, $seq2, $$aln_param{sub_mat});

#   Now, I should apply a shift to the sim_scr_mat (next version)

    my $fx_wt = 1.0 - $$aln_param{seq_wt};
    my $tmp_fx_mat = score_mat_shift ($fx_mat, $$aln_param{fx_shift});
    $total_mat = score_mat_add ($total_mat, $tmp_fx_mat, $fx_wt);
    undef ($tmp_fx_mat);
    $total_mat = score_mat_add ($total_mat, $sim_scr_mat, $$aln_param{seq_wt});


    my ($po, $pw, $qo, $qw);
    $po = $$aln_param {pgap_open};
    $pw = $$aln_param {pgap_widen};
    $qo = $$aln_param {qgap_open};
    $qw = $$aln_param {qgap_widen};

    my $pair_set1 =
        score_mat_sum_smpl (my $result_mat, $total_mat,
                            $po, $pw,
                            $qo, $qw,
                            $$aln_param{align_type});
    my $pair_set2 =
        score_mat_sum_smpl (my $result_mat, $total_mat,
                            $po * $up, $pw,
                            $qo * $up, $qw,
                            $$aln_param{align_type});
    my $pair_set3 =
        score_mat_sum_smpl (my $result_mat, $total_mat,
                            $po * $down, $pw,
                            $qo * $down, $pw,
                            $$aln_param{align_type});

    return ($pair_set1, $pair_set2, $pair_set3);
}

# ----------------------- get_pair_set_str_seq_prof_prof ------------
# Do an alignment based on both sequence and structure, but using
# profiles in both cases.
# For the sequence term, use profile-profile alignments
# For the sequence to structure, there is a precalculated score
# matrix, so it is not seen here.
# This function began its life on 19 Dec 2003.
sub get_pair_set_str_seq_prof_prof ( $ $ \% \@)
{
    my ($c1, $c2, $aln_param, $p_s) = @_;
    my $sub_mat = $$aln_param{sub_mat};
    my $profs = $$aln_param{profs};
    my $prof1 = $$profs{ coord_name ($c1)};
    my $prof2 = $$profs{ coord_name ($c2)};

    my @p_s;
    if (! (defined ($$p_s[0]))) {
        my $nothing =  undef();
        @p_s = ($nothing, $nothing, $nothing);
    } else {
        @p_s = @$p_s;
    }


    my $up   = 1.2;
    my $down = 0.8;

    my $index = coord_name($c1) . coord_name ($c2);

    my $fx_mat = $$aln_param{f_smat}{$index};

    if ( ! (defined ($prof1))) {
        my $names = 'profile not defined on ' . coord_name ($c1) .
            ' structure was ' . coord_name ($c2) . "\n";
        die $names;
    }
    if ( ! (defined ($fx_mat))) {
        my $names = coord_name ($c1) . ' ' . coord_name ($c2);
        die "fx_mat not defined on $names\n";
    }

    my $sim_scr_mat = score_mat_new (coord_size ($c1), coord_size($c2));
    my $total_mat   = score_mat_new (coord_size ($c1), coord_size($c2));
    score_prof_prof ($sim_scr_mat, $prof1, $prof2, $$aln_param{sub_mat});

#   Now, I should apply a shift to the sim_scr_mat (next version)

    my $fx_wt = 1.0 - $$aln_param{seq_wt};
    my $tmp_fx_mat = score_mat_shift ($fx_mat, $$aln_param{fx_shift});
    $total_mat = score_mat_add ($total_mat, $tmp_fx_mat, $fx_wt);
    undef ($tmp_fx_mat);
    $total_mat = score_mat_add ($total_mat, $sim_scr_mat, $$aln_param{seq_wt});


    my ($po, $pw, $qo, $qw);
    $po = $$aln_param {pgap_open};
    $pw = $$aln_param {pgap_widen};
    $qo = $$aln_param {qgap_open};
    $qw = $$aln_param {qgap_widen};

    my $pair_set1 =
        score_mat_sum_smpl (my $result_mat, $total_mat,
                            $po, $pw,
                            $qo, $qw,
                            $$aln_param{align_type}, $p_s[0]);
    my $pair_set2 =
        score_mat_sum_smpl (my $result_mat, $total_mat,
                            $po * $up, $pw,
                            $qo * $up, $qw,
                            $$aln_param{align_type}, $p_s[1]);
    my $pair_set3 =
        score_mat_sum_smpl (my $result_mat, $total_mat,
                            $po * $down, $pw,
                            $qo * $down, $pw,
                            $$aln_param{align_type}, $p_s[2]);
    return ($pair_set1, $pair_set2, $pair_set3);
}


# ----------------------- get_pair_set_struct -----------------------
# This does an alignment and returns an array of pairsets.
sub get_pair_set_struct ( $ $ $ \%)
{
    my ($c1, $c2, $f_smat, $aln_param) = @_;
    my $seq1 = coord_get_seq ($c1);
    my $seq2 = coord_get_seq ($c2);

    my $up   = 1.2;
    my $down = 0.8;

    my $index = coord_name($c1) . coord_name ($c2);


    my $fx_mat = $$f_smat{$index};

    my ($po, $pw, $qo, $qw);
    $po = $$aln_param {pgap_open};
    $pw = $$aln_param {pgap_widen};
    $qo = $$aln_param {qgap_open};
    $qw = $$aln_param {qgap_widen};


    my $total_mat = score_mat_shift ($fx_mat, $$aln_param{fx_shift});

    my $pair_set1 =
        score_mat_sum_smpl (my $result_mat, $total_mat,
                            $po, $pw,
                            $qo, $qw,
                            $$aln_param{align_type});
    my $pair_set2 =
        score_mat_sum_smpl (my $result_mat, $total_mat,
                            $po * $up, $pw,
                            $qo * $up, $qw,
                            $$aln_param{align_type});
    my $pair_set3 =
        score_mat_sum_smpl (my $result_mat, $total_mat,
                            $po * $down, $pw,
                            $qo * $down, $pw,
                            $$aln_param{align_type});

    return ($pair_set1, $pair_set2, $pair_set3);
}

# ----------------------- zero_shift_mat ----------------------------
# This shifts a substitution matrix, not an alignment score matrix.
sub zero_shift_mat ($ $)
{
    my ($sub_mat, $shift) = @_;
    for (my $i = 0; $i < 20; $i++) {
        for (my $j = $i; $j < 20 ; $j++) {
            my $t = sub_mat_get_by_i ($sub_mat, $i, $j) + $shift;
            sub_mat_set_by_i ($sub_mat, $i, $j, $t);
        }
    }
}

# ----------------------- catch_kill     ----------------------------
# The main thing is, if we get a KILL or TERM, to call exit so our
# buffers may get flushed and we can at least see output from slow,
# long running simplex things.
sub catch_kill
{
    print STDERR "Signal received\n";
    exit EXIT_FAILURE;
}

# ----------------------- kill_handlers  ----------------------------
# set up signal catchers so we can call exit() and die gracefully.
sub kill_handlers ()
{
    $SIG{INT } = \&catch_kill;
    $SIG{QUIT} = \&catch_kill;
    $SIG{TERM} = \&catch_kill;
}

# ----------------------- add_diag       ----------------------------
# Add the diagonal elements to the variable list
sub add_diag (\% \%)
{
    my ($s_args, $fixed) = @_;
    use Al_defaults qw (@letters);

    my $ini_pt = $$s_args{ini_pt};
    my $lower  = $$s_args{lower};
    my $upper  = $$s_args{upper};


    for (my $i = 0; $i <= $#letters; $i++) {
        my $a = $letters[$i];
        my $ini = sub_mat_get_by_c ($$fixed {sub_mat}, $a, $a);
        push (@{$$s_args{ini_pt}}, $ini);
        push (@{$$s_args{lower}},  -7);
        push (@{$$s_args{upper}}, 50);
    }
}

# ----------------------- use_diag       ----------------------------
# This is pretty tricky. We have the diagonal members of the
# substitution matrix sitting in a parameter
# array. Unfortunately, there are some other parameters up
# front. We want to hop over them and just put the rest into a
# substitution matrix.  If we agree to the rule that the diagonal
# parameters are at the tail of the matrix, then we can start
# there and walk backwards.

sub use_diag (\@ $)
{
    my ($var_param, $sub_mat) = @_;
    use Al_defaults qw (@letters);

    my $i = $#$var_param;

    for (my $j = $#letters; $j >= 0; $j--) {
        my $a = $letters[$j];
        if (! sub_mat_set_by_c ($sub_mat, $a, $a, $$var_param[$i--])) {
            die "use_diag broken. Prog bug. die evil trash\n"; }
    }
}


# ----------------------- opt_simplex -------------------------------
# Do the optimisation calculation.
sub opt_simplex ( \% \%)
{
    my $s_args      = shift;
    my $fix_param   = shift;
    use Paths qw ($MATRIX_DIR $MATRIX_FILE);
    use Al_defaults qw (@letters);


    my $coords      = $$fix_param {coords};
    my $test_coords = $$fix_param {test_coords};

    my %result;


    $$s_args {func}      = \&main::align_cost;
    $$s_args {result}    = \%result;
    $$s_args {fix_param} = $fix_param;

    if ( ! exists ($$fix_param{sub_mat_vary})) {
        warn "Go to the caller of opt_simplex() and define sub_mat_vary"; }

    if ($$fix_param{sub_mat_vary}) {
        print "Initial matrix\n", sub_mat_string ($$fix_param {sub_mat}) }
    use lib '../src/Simplex';
    use Simplex2;

    my $ret = simplex2 ($s_args);

    if ($$ret{success} == $SPLX_BROKEN) {
        warn "Simplex broke\n";
    } else {
        if ($$ret{success} == $SPLX_TOO_MANY) {
            warn "Simplex did not converge\n"}
        print "ncycle = $$ret{ncycle} and $$ret{restart} restarts\n",
              "best value $$ret{value} with..\n";

        while (@{$$s_args{names}}) {
            my $name = shift (@{$$s_args{names}});
            my $val  = shift (@{$$ret {best}});
            printf "%10s %10.4g\n", $name, $val;
        }
#       This is drastic, but safe. We have emptied the "names array".
#       To prevent future code trying to use it, let us wipe it out.
        undef @{$$s_args{names}};
        delete $$s_args{names};
    }

    my ($user, $system, $crap, $crap2) = times();
    printf "I took %d:%d min user and %.0f:%.0f min sys time\n",
    $user / 60, $user % 60, $system / 60, $system % 60;

    return $ret;
}


# ----------------------- rand_int_non_zero -------------------------
# Give me a random, non-zero integer between (and including) $a
# and $b
sub rand_int_non_zero ($ $)
{
    use Wutil qw (rand_int);
    my ($a, $b) = @_;
    my $ret = 0;
    while ($ret == 0) {
        $ret = rand_int ($a, $b); }
    return $ret;
}

# ----------------------- rand_int  ---------------------------------
# Give me a random integer between (and including)
# two args. Given (-1, 2), it will return -1, 0, 1, or 2.
sub rand_int ($  $)
{
    my ($min, $max) = @_;
    if ($max < $min) {
        my $t = $min;
        $min = $max;
        $max = $t;
    }
    if ($max == $min) {
        warn "silly min/max $min, $max given to rand_int\n"; return undef;}
    my $diff = $max - $min + 1;
    my $t = rand $diff;
    $t = int ($t);
    return ($t + $min);
}

# ----------------------- bounded_mover -----------------------------
# Move function for Monte Carlo
# We move 'to_move' randomly chosen parameters by up to
# 'max_step'. Reject moves outside specified upper or lower
# bounds by asking for a new random number. If there are no
# specified bounds, use defaults.
sub bounded_mover (\@ $)
{
    my ($params, $movehash) = @_;
    use Wutil qw (rand_int rand_int_non_zero);
    use Al_defaults qw ($MAX_PARAM $MIN_PARAM);
    my $mpar = $$movehash {to_move};
    my $max_step = $$movehash {max_step};
    my $lower = $$movehash{lower};
    my $upper = $$movehash{upper};

#   If we are allowed to move up to $mpar parameters, change to move
#   between 1 and $mpar parameters
    $mpar = rand_int (1, $mpar);

    for (my $n = 0; $n < $mpar; $n++) {
        my $ndx  = rand_int (0, $#$params);  # Pick param randomly
        my $ok = undef;
        my $t;
        while (! $ok) {
            my $step = rand_int_non_zero(-$max_step, $max_step);
            $t = $$params[$ndx] + $step; # possible move
            my $bound;
            if (! defined ($bound = $$lower[$ndx])) {
                $bound = $MIN_PARAM; }
            if ($t <= $bound) {
                next ; }
            undef $bound;
            if (! defined ($bound = $$upper[$ndx])) {
                $bound = $MAX_PARAM; }
            if ($t >= $bound) {
                next ; }
            $ok = 1;
        }
        $$params[$ndx] = $t;
    }
}

# ----------------------- usage     ---------------------------------
sub
usage ()
{
    print STDERR
        "$0: earlier versions had default arguments.\n",
        "This one does not.\n",
        "usage: $0 [options] pair_file parameter_file\n";
    return (undef);
}

# ----------------------- init_main      ----------------------------
# Initialisation steps that are common to the various scripts.
# Return undef if we break.
sub init_main( \% \% \%)
{
    my ($s_args, $fixed, $test_fix_param) = @_;

    use Paths qw ($MATRIX_DIR $MATRIX_FILE);
    use Paths qw ($DFLT_PAIR_FILE $DFLT_PARAM_FILE);
    use Paths qw (@DFLT_BIN_DIR);

    my $param_file;
    my $pair_file;

    if ($#ARGV != 1) {
        return (usage()); }

    if ( ! -f ($pair_file = $ARGV[0])) {
        warn "Pair_file \"$pair_file\" not found\n"; return undef }
    if ( ! -f ($param_file = $ARGV[1])) {
        warn "Parameter file \"$param_file\" not found\n"; return undef}

    if ( get_param_list ($param_file, $s_args, $fixed) == EXIT_FAILURE ) {
        warn "Broke reading parameters from $param_file\n"; return undef; }

    my $matname;
    if (! defined ($matname = $$s_args { sub_mat })) {
        if ( ! -d $MATRIX_DIR) {
            warn "\"$MATRIX_DIR\" matrix dir not found.\n"; return undef; } }

#   The substitution matrix might be local or maybe we have to find it.
    my $sub_mat;
    if ($matname) {
        if (! -f $matname) {
            my $t = "$MATRIX_DIR/$matname";
            if ( ! -f $t ) {
                warn "Cannot find sub matrix $matname\n"; return undef;}
            $matname = $t;
        }
        if ( !($sub_mat = sub_mat_read ($matname))) {
            warn "Fail reading sub matrix from $matname\n"; return undef;}
        $$fixed {sub_mat} = $sub_mat;
    } else {
        warn "No matrix name defined. Proceeding carefully.\n";
        delete $$fixed {sub_mat};
    }


#   There is a bit of a trick here.
#   In the parameter file, "sub_mat_vary" must be yes/no. Here, we check the
#   value and, if we are happy, we set it to 1/undef and overwrite the
#   previous value.

    if (exists ($$fixed{sub_mat_vary})) {
        print "Varying of submat has been explicitly set ";
        if (lc ($$fixed{sub_mat_vary}) eq 'no') {
            print "to \'no\'\n";
            $$fixed{sub_mat_vary} = undef;
        }
        elsif (lc ($$fixed{sub_mat_vary}) eq 'yes') {
            print "to \'yes\'\n";
            $$fixed {sub_mat_vary} = 1;
        }
        else {
            print "but to \'$$fixed{sub_mat_vary}\'\n";    return undef; }
    } else {
        $$fixed {sub_mat_vary} = undef;
    }

#   From here on, $dipep being defined indicates that we are going to use
#   dipeptide information. After we read the dipeptide information into the
#   dipep structure, put it into the $$fixed hash.

    my $dipep;
    if (defined ((my $dipep_file = $$s_args {dipep_file}))) {
        if ( ! ($dipep = dpt_read ($dipep_file))) {
            warn "Failed reading dipep info from $dipep_file: $!\n";
            return undef;
        } else {
            print "Read ", dpt_get_n ($dipep), " lines from $dipep_file\n";
            $$fixed {dipep} = $dipep;
        }
    } else {
        print "Not using dipeptide information.\n";
    }

#   The next section is for reading the coordinates (coord) data and calling
#   split_coord() to put some of them in the array for testing.
    {
        my @pairlist;
        if ((@pairlist = get_pair_list ($pair_file)) < 2) {
            warn "Fail to get pairs from $pair_file\n"; return undef; }

        my $coords = $$fixed {coords};
        if ( !check_pdb_files(@pairlist, @DFLT_BIN_DIR, @$coords)) {
            return undef; }


        my $test_coords = $$test_fix_param {coords};
        split_coord (@$coords, @$test_coords);

        if ( defined (@$test_coords)) {
            %$test_fix_param = %$fixed;
            $$test_fix_param {coords} = $test_coords;
            $$s_args {test_fix_param} = $test_fix_param;
        } else {
            print "no test array being used\n";
            undef $test_coords;
            undef $$test_fix_param {coords};
            undef %$test_fix_param;
        }
    }

    foreach my $i (@{$$s_args{names}} ) {
        if ( defined ($$fixed {$i})) {
            warn "$i is given as both fixed and variable parameter. Fix\n";
            return undef;
        }
    }
    return 1;
}

# ----------------------- init_profiles -----------------------------
# Read up and store all of the sequence profiles.
sub init_profiles ( \% \@)
{
    my ($profs, $coords) = @_;
    use FindBin;
    use lib "$FindBin::Bin";
    use Wutil qw(get_pair_list);
    use Paths qw ($PROFILE_DIR);

    my $error = undef;
    my $last = 'none yet';
    for my $x (@$coords) {
        for my $c ($ {$x} [0], $ {$x} [1]) {
            my $orig = coord_name ($c);
            my $name = $orig;
            substr ($name, 0, 4) = lc (substr($orig, 0, 4));
            $name = "$PROFILE_DIR/$name.prof";
            my $t;
            if (! ($t = blst_chk_read ("$name"))) {
                warn "profile reading failed on $name last ok was $last\n";
                $error++;
            }
            $last = $name;
            $$profs{$orig} = $t;
        }
    }
    if ($error) {
        return undef;
    } else {
        return 1; }
}

1;
