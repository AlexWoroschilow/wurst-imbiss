#!/usr/bin/perl
#$ -S /usr/bin/perl
##$ -clear
#$ -w e
##$ -l arch=glinux -l short=0
#$ -p -50
#$ -cwd
#$ -j y
##$ -m e -M torda@zbh.uni-hamburg.de
#$ -q hpc.q
# 18 March 2002
# Most of the other scripts here are for testing with
# deliberately crazy values.
# This one uses the installed, working binary and serious values
# for all parameters.
# This one goes back to the parameters which do *not* favour
# extended alignments.
# rcsid = $Id: libsrch.pl,v 1.1 2007/09/28 16:57:01 mmundry Exp $

=pod

=head1 NAME

libsrch.pl - Given a sequence, align it to a library of templates

=head1 SYNOPSIS

libsrch.pl [options] I<seq_file> I<struct_lib_list> I< S<[ phd_file ]> >


=head1 DESCRIPTION

Given a sequence, align it to every member of a library of
templates.  The sequence is given by I<seq_file>.  The library is
a list of protein names listed in <struct_lib_list>. The last
argument is optional and contains the name of a file with
secondary structure predictions. If it is not present, the script
will look in the directory containing I<seq_file> and strip off
anything that looks like a file extension. Then, it will append
B<.phd> and try to open that (so /boo/bar/1abc.seq gives
/boo/bar/1abc.phd).

This is the may 2006 version which uses

=over

=item * Thomas Huber's fx9 score function with profiles

=item * profile to profile alignments

=item * the fragment to fragment score function

=back

=head2 FILE FORMAT

The list of files which make up the template library is in a
simple format. The script will try to read anything that looks

like a four-letter protein name + chain id from the first
column. Leading white space is ignored. A valid form would look
like

   1abc_
   2qrsB
   1xyz  This text after first column is ignored

=head2 Changing library and templates.

Typically, a first run will be made with whatever library we are
using. However, one will often want to add extra .bin files
for a particular sequence. To do that,

=over

=item *

Add the new file names to the list of proteins and give it a name
like F<mylist>.

=item *

Make a directory with a name like I<templates> and put the extra
F<.bin> files in there.

=item *

Run the script with the B<-t> option like:

  perl libsrch.pl -t templates blahblah.seq mylist

=back

=head2 OPTIONS

=over

=item B<-a> I<N>

Print out details of the best I<N> alignments.

=item B<-d> I<modeldir>

Save final models in I<modeldir> instead of the default
directory, B<modeldir>.

=item B<-h> I<N>

After doing alignments, a default number (maybe 50) will be
printed out. Alternatively, you can ask for the I<N> best scoring
alignments to be printed out.

=item B<-i> I<N>

For the best I<N> alignments, write each, in modeller/PIR format,
to a file in F<pir_dir>. Currently, the name F<pir_dir> is
fixed. If this option is not set, the directory F<pir_dir> will
be created and the default number of alignments will be the same
as the number of models constructed. To turn off this output, set
I<N> to zero.

=item B<-m> I<N>

After alignments, I<N> models will be built, otherwise, a small
default number will appear. Set I<N> to zero if you do not want
any models.

=item B<-p>

Keep the temporary blast profile and put it in a name like F<seq.prof>.

=item B<-t> I<dir1[,dir2,...]>

Add I<dir1> to the list of places to look for template
files. This is a comma separated list, so you can add more
directories.

=back

=head1 OUTPUT

In all output

=over

=item *

B<SW> refers to the result of the second Smith and Waterman.

=item *

B<NW> refers to the result of the Needleman and Wunsch.

=item *

B<cvr> or B<cover> refers to "coverage".  This is the fraction of
sequence residues which are aligned to part of a
structure. Continuing, B<S<sw cvr>> refers to coverage from the
Smith and Waterman calculation.

=item *

The script prints out the coverage in a somewhat pictorial form
which might look like

   ----XXXXX---XXX

where the X's mean a residue was aligned.

=back

=head1 MODELS

Models will be written in PDB format for the best few
sequences. They will get written to a directory called
F<modeldir>. Maybe this should be made an option.

=head1 OPERATION

Currently the script does

=over

=item Smith and Waterman step 1

This is a locally optimal alignment.

=item Smith and Waterman step 2

This is another locally optimal alignment, but forced to pass
through the same path as the first one. It provides a small
extension to the alignment.

=item Needleman and Wunsch

This is a globally optimal alignment, but forced to pass through
the preceding Smith and Waterman.

=back

=head1 QUESTIONS and CHANGES

=item *

The selection of which scores to print is ver arbitrary.

=item *

the coverage picture is very ugly. It could be beautified.

=item *

The coverage picture corresponds to the Smith and
Waterman. Perhaps it should be the Needleman and
Wunsch. Obviously, both are possible, but just a bit ugly.

=cut

use warnings;
use lib "$FindBin::Bin../blib/arch";
use lib "$FindBin::Bin../blib/lib";

use FindBin;
#use lib "$FindBin::Bin/../src/Wurst/blib/arch";
#use lib "$FindBin::Bin/../src/Wurst/blib/lib";

use Wurst;

use vars qw ($MATRIX_DIR $MATRIX_FILE $PARAM_DIR
             $RS_PARAM_FILE $FX9_PARAM_FILE );
use vars qw ($TOPWURST $TOPSTRUCT);


# ----------------------- Paths   -----------------------------------
# These used to be elegantly kept in a separate file. Unfortunately,
# it is rather inconvenient to have the information somewhere else
# and especially so when we are running in a batch system which
# moves scripts around.

*TOPSTRUCT = \'/projects/bm';
*TOPWURST  = \"$ENV{HOME}/c/wurst/default/wurst";
# $TOPWURST="$FindBin::Bin/..";         # Best, but not under batch sys
$MATRIX_DIR = "$TOPWURST/matrix";
$PARAM_DIR  = "$TOPWURST/params";

use strict;

use POSIX qw(EXIT_SUCCESS EXIT_FAILURE);


# ----------------------- Defaults  ---------------------------------
# These are numbers you might reasonably want to change.
# They should (will) be changeable by options.
use vars qw ($N_BRIEF_LIST $N_MODELS $N_ALIGNMENTS $N_PIR);
$N_BRIEF_LIST   = 100;
$N_MODELS       =  25;
$N_ALIGNMENTS   =  50;
$N_PIR          =  50;
use vars qw ($modeldir $DFLT_MODELDIR $pir_dir $DFLT_PIRDIR);
*DFLT_MODELDIR = \ 'modeldir';
$modeldir = $DFLT_MODELDIR;
*DFLT_PIRDIR   = \ 'pir_dir';
$pir_dir  = $DFLT_PIRDIR;


# ----------------------- Alignment Constants -----------------------
# These are declared globally, and set by set_params().
use vars qw (
             $align_type
             $sw1_pgap_open
             $sw1_qgap_open
             $wt1 $wt2 $wt3
             $sw1_pgap_widen
             $sw1_qgap_widen
             $sw1_sec_pnlty
             $off1
             );

use vars qw ( $MATRIX_FILE );

# These parameters will be used for extending alignments via a
# Needleman and Wunsch

use vars qw (
             $nw_pgap_open
             $nw_qgap_open
             $nw_pgap_widen
             $nw_qgap_widen
             $nw_sec_pnlty
);
use vars qw ( $SW_RS_GEOG $SW_RS_NSG $SW_RS_OPEN $SW_RS_WIDEN $SW_SCR2);
( $SW_RS_GEOG, $SW_RS_NSG, $SW_RS_OPEN, $SW_RS_WIDEN, $SW_SCR2) =
( -0.678,        -0.021,    -2.04,       -0.188, 0.0654);
use vars qw( @DFLT_STRUCT_DIRS  @PROFILE_DIRS @DFLT_VEC_DIRS
             $CLASS_DIR $CLASS_FILE $ABS_ERROR
             $phd_suffix $bin_suffix $prof_suffix $vec_suffix);
*DFLT_STRUCT_DIRS = ["$TOPSTRUCT/pdb90_bin"];
*PROFILE_DIRS     = ["$TOPSTRUCT/pdb90_prof"];
*DFLT_VEC_DIRS    = ["$TOPSTRUCT/pdb90_vec_5mer8"];
*bin_suffix       = \'.bin';
*prof_suffix      = \'.prof';
*phd_suffix       = \'.phd';
*vec_suffix       = \'.vec';
*MATRIX_FILE      = \'741.out';
*CLASS_DIR        = *PARAM_DIR;
*CLASS_FILE       = \'5mer_seq_angle4/r8.influ-o-data-1';
*FX9_PARAM_FILE = \'g9.param';
*RS_PARAM_FILE    = \'cc_allat_p891+0.param';
*MATRIX_FILE      = \'blosum62.50';
*ABS_ERROR        = \ 0.4; 

# ----------------------- set_params    -----------------------------
# This gets its own function because it can be more complicated
# if, in the future, we have a version depending on various
# options like whether or not we have secondary structure
# information.
# pgap controls penalties in the sequence, qgap in the structure.

sub set_params ()
{
    *sw1_pgap_open  =  \  1.382;
    *sw1_pgap_widen =  \  0.3089;
    *sw1_qgap_open  =  \ $sw1_pgap_open;
    *sw1_qgap_widen =  \ $sw1_qgap_widen;
    *wt1            =  \  0.2399;
    *wt2            =  \  0.6372;
    *wt3            =  \ (1 - ($wt1 + $wt2));
    *off1           =  \  0.4949;
    *nw_pgap_open  =   \  $sw1_pgap_open;
    *nw_qgap_open  =   \  $sw1_qgap_open;
    *nw_pgap_widen =   \  $sw1_pgap_widen;
    *nw_qgap_widen =   \  $sw1_qgap_widen;
    *nw_sec_pnlty  =   \  0.0;
}



# ----------------------- get_prot_list -----------------------------
# Go to the given filename and get a list of proteins from it.
sub get_prot_list ($)
{
    my $f = shift;
    my @a;
    if ( ! open (F, "<$f")) {
        print STDERR "Open fail on $f: $!\n";
        return undef;
    }

    while (my $line = <F>) {
        chomp ($line);
        my @words = split (' ', $line);
        if (! defined $words[0]) { next;}
        $line = $words[0];
        $line =~ s/#.*//;            # Toss comments away
        $line =~ s/\..*//;           # Toss filetypes away
        $line =~ s/^ +//;            # Leading and
        $line =~ s/ +$//;            # trailing spaces.
        if ($line eq '') {
            next; }
        substr ($line, 0, 4) = lc (substr ($line, 0, 4)); # 1AGC2 to 1agc2
        if (length ($line) == 4) {  # Convert 1abc to 1abc_
            $line .= '_'; }
        push (@a, $line);
    }
    close (F);
    return (@a);
}

# ----------------------- get_path  ---------------------------------
# We have a filename and a list of directories where it could
# be. Return the path if we can find it, otherwise return undef.
sub get_path (\@ $)
{
    my ($dirs, $fname) = @_;
    foreach my $d (@$dirs) {
        my $p = "$d/$fname";
        if ( -f $p) {
            return $p; }
    }
    return undef;
}

# ----------------------- check_dirs --------------------------------
# Given an array of directory names, check if each one
# exists. Print something if it is missing, but do not give
# up. It could be that there is some crap in the command line
# args, but all the important directories are really there.
# This function is potentially destructive !
# If a directory does not seem to exist, we actually remove it
# from the array we were passed.  This saves some futile lookups
# later on.
sub check_dirs (\@)
{
    my $a = shift;
    my $last = @$a;
    for (my $i = 0; $i < $last; $i++) {
        if ( ! -d $$a[$i]) {
            print STDERR "$$a[$i] is not a valid directory. Removing\n";
            splice @$a, $i, 1;
            $last--;
            $i--;
        }
    }
}

# ----------------------- check_files -------------------------------
# We are given an array of directories and and array of protein
# names and an extension.
# Check if all the files seem to be there.
sub check_files (\@ \@ $)
{
    my ($dirs, $fnames, $ext) = @_;
    my $errors = 0;
    foreach my $f (@$fnames) {
        my $name = "$f$ext";
        if (! get_path (@$dirs, $name)) {
            $errors++;
            print STDERR "Cannot find $name\n";
        }
    }
    return $errors;
}

# ----------------------- usage   -----------------------------------
sub usage ()
{
    print STDERR "Usage: $0 seq_file struct_file \n";
    exit (EXIT_FAILURE);
}

# ----------------------- get_scores --------------------------------
sub get_scores ($ $ $ $ $)
{
    my ($pair_set, $coord, $seq, $rescore_params, $to_use) = @_;

    my ($scr_tot, $coverage, $score1, $score2, $geo_gap, $score1_gap);
    my ($str1, $crap);
    my ($open_cost, $widen_cost, $nseq_gap);
    ($score1_gap, $score1) = pair_set_score( $pair_set );
    ($str1, $crap) =
        pair_set_coverage ($pair_set, seq_size ($seq), coord_size ($coord));
    $coverage = ($str1 =~ tr/1//);       # This is coverage as an integer
    $coverage = $coverage / seq_size ($seq); #and as fraction of sequence

    my ($k_gap_geo, $k_seq_gap, $k_str_gap, $k_str_wdn, $k_scr2);
    if ($to_use eq 's_and_w') {
        ( $k_gap_geo,  $k_seq_gap,  $k_str_gap,  $k_str_wdn, $k_scr2) =
        ( $SW_RS_GEOG, $SW_RS_NSG,  $SW_RS_OPEN, $SW_RS_WIDEN, $SW_SCR2);
    } else {
        ( $k_gap_geo,  $k_seq_gap,  $k_str_gap,  $k_str_wdn, $k_scr2) =
        ( $SW_RS_GEOG, $SW_RS_NSG,  $SW_RS_OPEN, $SW_RS_WIDEN, $SW_SCR2);
    }

    if ($coverage  < .05 ) {
        $score2 = 0;
        $geo_gap = 0;
        $nseq_gap = 0;
        $open_cost = 0;
        $widen_cost = 0;
    } else {
        my $model = make_model ($pair_set, $seq, $coord); # Build model
        $score2 = score_rs ($model, $rescore_params);
        ($crap, $geo_gap, $crap, $nseq_gap) = coord_geo_gap ($model, 1, 10);
        ($open_cost, $widen_cost) = pair_set_gap($pair_set, 1, 1);
    }

    $scr_tot = $score1 +
               $k_gap_geo * $geo_gap +
               $k_seq_gap * $nseq_gap +
               $k_str_gap * $open_cost +
               $k_str_wdn * $widen_cost +
               $k_scr2    * $score2;
    return ( $scr_tot, $coverage, $score1, $score1_gap,
             $score2, $nseq_gap, $open_cost);
}


# ----------------------- bad_exit ----------------------------------
# This will run in a server, so if something goes wrong, we
# should at least mail back an indication.  The single parameter
# should be the error message returned by the function which was
# unhappy.
# Should we print to stderr or stdout ?
# This should not matter since we have grabbed both file handles.
sub bad_exit ( $ )
{
    my $msg = shift;
    print STDERR "Error: \"$msg\"\n";
    exit (EXIT_FAILURE);
}

# ----------------------- get_alt_scores ---------------------------------
# calculates scores on random paths through the scoring matrix
# parameters: number_of_paths/scores, scoring_matrix, pair_set_of_optimal_path
# return: the scores
sub get_alt_scores($ $ $)
{
   my ($num_scrs, $scr_mat, $pair_set) = @_;
   my @scr_fin;

   for (my $i = 0; $i < $num_scrs; $i++) {
       $scr_fin[$i] = find_alt_path_score_simple ($scr_mat, $pair_set); }

   return \@scr_fin;
}

# ----------------------- normalize_alt_scores ------------------------------
#
sub normalize_alt_scores($)
{
   my ($scrs_ref) = @_;
   my $mean = 0.0;

   foreach my $scr (@{$scrs_ref}) {
       $mean += $scr }
   $mean /= @$scrs_ref;

   my $deviation = 0.0;

   foreach my $scr (@$scrs_ref) {
       my $tmp = $scr - $mean;
       $deviation += ($tmp * $tmp);
   }
   $deviation /= (@$scrs_ref - 1);
   $deviation = sqrt($deviation);

   return ($mean, $deviation);
}

# ----------------------- line_splitter  ------------------------------------
# Split a long line (first argument) into newline separated pieces, each of
# length "len" (second argument). Return the new string.
# Removes trailing white space.
sub line_splitter ($ $)
{
    my ($in, $len) = @_;
    my $out;
    my $s;
    do {
        $s = substr ($in, 0, $len);
        $out .= "$s\n";
        substr ($in, 0, $len) = '';
    } while ($s);
    $out =~ s/\s+$//;   # Trailing white space.
    return $out;
}


# ----------------------- pir_from_raw_model --------------------------------
# This is taken from Jim's scripts. It returns a piece of text which is the
# wurst alignment in pir form.
sub pir_from_raw_model ( $ $ $ $ $ $ ) {
    my ($snam, $seq, $tnam, $tstrn, $coord, $pairset ) = @_;
    my @p;
    my $t_model_seq;
    my $coord_size = coord_size ($coord);
    my ($s_aligned, $t_aligned) = pair_set_coverage($pairset, seq_size($seq), $coord_size);

    my @txt_algnment = split "\n",(pair_set_string($pairset, $seq, coord_get_seq($coord)));
    my ($al_seq, $al_tem, $sseq, @sq);
    $sseq = seq_print($seq);
    $sseq =~ s/>.+\s+//;
    $sseq =~ s/\s//g;
    chomp $sseq;
    @sq = split '',$sseq;
    my $tseq = seq_print(coord_get_seq($coord));
    $tseq =~ s/>.+\s+//;
    $tseq =~ s/\s//g;
    chomp $tseq;
    my @tsq = split '',$tseq;
    # reconstruct pairs
    my (@cov_s,@cov_q);
    @cov_s = split '', $s_aligned;
    @cov_q = split '', $t_aligned;
    my $s = shift @cov_s;
    my $q = shift @cov_q;
    while (scalar @sq or scalar @tsq) {
        if ($s) {
            if ($q) {
                $al_seq.=shift @sq;
                $al_tem.=shift @tsq;
                $s = shift @cov_s;
                $q = shift @cov_q;
            } else {
                $al_seq.="-";
                $al_tem.=shift @tsq;
                $q=shift @cov_q;
            }
        } elsif ($q) {
            $al_tem.="-";
            $al_seq.=shift @sq;
            $s=shift @cov_s;
        } else {
            my ($tins,$sins);
            while ((scalar @sq) and (not $s)) {
                $sins.=(shift @sq);
                $s = shift @cov_s;
            };
            while ((scalar @tsq) and not $q) {
                $tins.=(shift @tsq);
                $q = shift @cov_q;
            }
            if (defined($tins)) {
                $al_seq.=("-"x(length $tins));
                $al_tem.=$tins;
            }
            if (defined($sins)) {
                $al_seq.=$sins;
                $al_tem.=("-"x(length $sins));
            }
        }
    }
    $al_seq = line_splitter ($al_seq, 75);
    $al_tem = line_splitter ($al_tem, 75); 
    # this case is really easy but we do have to use the original template model
    my $m_start = model_pdb_num ($coord, 0);
    my $m_end   = model_pdb_num ($coord, $coord_size - 1);
    my ($a_start, $a_end) = (1, seq_size ($seq));
    # Sequence comes first
    push @p, ">P1; $snam";
    push @p, "sequence:$a_start:1: :$a_end:\@:.:.:.:.";
    push @p, "".(uc $al_seq)."*","";
    my $chain = substr ($tnam, 4, 1);
    if ($chain eq '_') {
        $chain = '@'; }
    # template sequence last so we can use ALIGN_WHAT=last in top script
    # to expand alignment to whole PDB file
    push @p, ">P1; $tnam";
    push @p, "structure:$tstrn:$m_start:$chain:$m_end:$chain:.:.:.:.";
    push @p, "".(uc $al_tem)."*","";
    return (join "\n", @p,"");
}


# ----------------------- do_align ----------------------------------
# This does the alignment. Although it takes secondary structure
# data as an argument, we are not yet using this in our server
# script.
sub do_align ($ $ $ $ $ $ $ $)
{
    my ($profile, $template, $tmplt_prf_pth, $tmplt_vec_pth,
        $fx_params, $submat,
        $seqvec, $rescore_params)
        = @_;

    my $coord  = coord_read ($template) || bad_exit ("Fail on $template");
    my $tmplt_prof = blst_chk_read ($tmplt_prf_pth) || bad_exit("Fail on profile");
    my $tmplt_vec  = prob_vec_read ($tmplt_vec_pth) || bad_exit("Fail on vector");
#   Now we have all files we need. Start scoring.
#   Begin by giving us empty score matrices.

    my $seq        = seqprof_get_seq ($profile);
    my $seq_size   = seq_size ($seq);
    my $coord_size = coord_size ($coord);

    my $fx_mat        = score_mat_new ($seq_size, $coord_size);
    my $sim_scr_mat   = score_mat_new ($seq_size, $coord_size);
    my $seq_strct_mat = score_mat_new ($seq_size, $coord_size);

    if ( ! score_fx_prof ($fx_mat, $profile, $coord, $fx_params)) {
        print STDERR "fx_prof error on ", coord_name ($coord), "\n"; return undef}
    $fx_mat = score_mat_scale ($fx_mat, $wt1);
    $fx_mat = score_mat_shift ($fx_mat, $off1);


    if ( ! score_prof_prof ($sim_scr_mat, $profile, $tmplt_prof, $submat)) {
        print STDERR "prof_prof error on ", coord_name ($coord), "\n"; return undef ;}

    $sim_scr_mat = score_mat_scale ($sim_scr_mat, $wt2);

    if ( ! score_pvec ($seq_strct_mat, $seqvec, $tmplt_vec) ) {
        print STDERR "seq_strct error on ", coord_name ($coord), "\n"; return undef ;}
    $seq_strct_mat = score_mat_scale ($seq_strct_mat, $wt3);
    my $total_mat = score_mat_add ($fx_mat, $sim_scr_mat, 1.0);
    $total_mat = score_mat_add ($total_mat, $seq_strct_mat, 1.0);


#   This actually does the alignment.
    my ($sw_scr_tot, $sw_coverage, $sw_score1, $sw_score1_gap, $sw_score2,
        $sw_seq_gap, $sw_strct_gap);
    my ($nw_scr_tot, $nw_coverage, $nw_score1, $nw_score1_gap, $nw_score2,
        $nw_seq_gap, $nw_strct_gap);
    my $sw_pair_set =
        score_mat_sum_smpl (my $result_mat, $total_mat,
                            $sw1_pgap_open, $sw1_pgap_widen,
                            $sw1_qgap_open, $sw1_qgap_widen,
                            $S_AND_W);

    ($sw_scr_tot, $sw_coverage, $sw_score1, $sw_score1_gap, $sw_score2,
        $sw_seq_gap, $sw_strct_gap) =
    get_scores ($sw_pair_set, $coord, $seq, $rescore_params, 's_and_w');

    my $num_scrs = 3000;
    my $alt_scrs_ref = get_alt_scores($num_scrs, $total_mat, $sw_pair_set);

    my ($mean, $deviation) = normalize_alt_scores($alt_scrs_ref);
    undef ($alt_scrs_ref);
    my $z_scr;
    if ($deviation != 0) {
        $z_scr = ($sw_scr_tot - $mean) / $deviation; }
    else {
        $z_scr = 0.0; }                          # Should not really happen
#   If the alignment is tiny, one can get a ridiculous z-score
    if ($sw_coverage < 0.03) {                   # silently wipe these guys out
        $z_scr = 0.0; }
    if ($sw_coverage * seq_size ($seq) < 10) { # and wipe out anthing
        $z_scr = 0.0; }                             # less than 10 residues


    my $nw_pair_set =
        score_mat_sum_sec (my $result_mat2, $total_mat,
                           $coord, $nw_sec_pnlty,
                           $nw_pgap_open, $nw_pgap_widen,
                           $nw_qgap_open, $nw_qgap_widen,
                           $N_AND_W, $sw_pair_set);
    ($nw_scr_tot, $nw_coverage, $nw_score1, $nw_score1_gap, $nw_score2,
     $nw_seq_gap, $nw_strct_gap) =
    get_scores ($nw_pair_set, $coord, $seq, $rescore_params, 'n_and_w');
    my @r =
        ($sw_scr_tot, $nw_pair_set, $nw_scr_tot,
         $sw_coverage, $nw_coverage, $sw_score1, $sw_score2,
         $nw_score1, $nw_score2, $sw_pair_set,
         $z_scr);
    return (\@r);
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

# ----------------------- do_lib  -----------------------------------
# Walk over a library, doing alignments and saving interesting
# scores. The definition of interesting is a bit arbitrary.
# There is one very non-obvious coding trick.  We need to be able
# to pass the score information into the sorting functions. We
# could put everything into a big, two-dimensional array, but we
# can avoid copying data. Instead, we invent a package and put
# results into @r::r. The downside is that we have to manually
# free it up at the end by calling undef().
# The $formatflag argument is used for cases like the livebench
# server for whom we have to do a bit of file re-writing.
sub do_lib (\@ \@ $ $ $ $)
{
    my ($structlist, $struct_dirs, $profile_file, $title, $seqvec,
        $formatflag) = @_;

    my $pfile = "$PARAM_DIR/$RS_PARAM_FILE";
    my $rescore_params = param_rs_read ($pfile) || bad_exit( "Rescore params");

    my $matname = "$MATRIX_DIR" . '/' . "$MATRIX_FILE";
    my $submat = sub_mat_read ($matname) || bad_exit( "Fail on $matname");
#   zero_shifting of the substitution matrix was done here.
#   Maybe should be tried again later

    my ($tmp_p, $fx_params);
    $tmp_p = "$PARAM_DIR/$FX9_PARAM_FILE";
    $fx_params = param_fx_read ($tmp_p) || bad_exit ( "Fail on $tmp_p");

    my $profile = blst_chk_read ($profile_file);
    if ( ! $profile) {
        warn "Failed reading profile $profile_file\n"; return EXIT_FAILURE; }
    my $seq     = seqprof_get_seq ($profile);
    for (my $i = 0; $i < @$structlist ; $i++) {
        my ($template, $tmplt_prf_pth, $tmplt_vec_pth);
        $template   = get_path (@$struct_dirs, $$structlist[$i] . $bin_suffix);
        $tmplt_prf_pth = get_path (@PROFILE_DIRS, $$structlist[$i] . $prof_suffix);
        if ( ! $tmplt_prf_pth) {
            warn "NO profile for $$structlist[$i]\n"; return EXIT_FAILURE; }
        $tmplt_vec_pth = get_path (@DFLT_VEC_DIRS, "$$structlist[$i]$vec_suffix");
        if ( ! $tmplt_vec_pth) {
            warn "NO vect file for $$structlist[$i]\n"; return EXIT_FAILURE; }
        $r::r[$i] = do_align ( $profile, $template, $tmplt_prf_pth,
                               $tmplt_vec_pth,
                               $fx_params, $submat, $seqvec, $rescore_params);
    }

    my @indices;
    for ( my $i = 0; $i < @$structlist; $i++) {
        $indices[$i] = $i; }


    @indices = sort {
        $r::r[$b][10] <=> $r::r[$a][10];
    } @indices;

    print "\n",
"sw refers to Smith and Waterman alignment, nw refers to Needleman & Wunsch\n",
"z scr : the z-score of the alignment with 1000 alternative alignments\n",
"sw scr: the combined results of score function + gaps for sw alignment\n",
"sw cvr: coverage (fraction of sequence) accounted for by sw alignment\n",
"nw cvr: coverage                        accounted for by nw alignment\n",
"sw1   : score of sw alignment in first score function\n",
"sw2   : score of sw alignment in second score function\n",
"____________ Summary of best templates   _________________________________\n";

    my $todo = (@$structlist > $N_BRIEF_LIST ? $N_BRIEF_LIST : @$structlist);
    printf "%8s %8s %8s %6s %6s %8s %8s\n",
    'struct',  'z scr', 'sw scr', 'sw cvr', 'nw cvr', 'sw1', 'sw2';
    for (my $i = 0; $i < $todo; $i++) {
        my $idx = $indices [$i];
        my $a = $r::r[$idx];
        printf "%8s %8.2g %8.1f %6.2f %6.2g %8.2g %8.2g\n",
               $$structlist[$idx],
               $$a[10], $$a[0], $$a[3], $$a[4],
               $$a[5], $$a[6];
    }
    print "\n";


#   If the sequence is not too big, this output is informative.
#   If the output is too big, the long lines cause errors !
    if ( seq_size ($seq) > 550) {
        print "Sequence length ", seq_size ($seq), "too big.\n",
        "Not going to print out \"Summary of coverage of query sequence\"\n";
    } else {
        print
"____________ Summary of coverage of query sequence _______________________\n";

        for (my $i = 0; $i < $todo; $i++) {
            my $idx = $indices [$i];
            my $a = $r::r[$idx];
            my $csize =
                coord_size ( coord_read (
                  get_path (@$struct_dirs, $$structlist[$idx] . $bin_suffix)));
            my ($str, $crap) =
                pair_set_coverage ($$a[9], seq_size ($seq), $csize);
            print "S & W coverage with $$structlist[$idx]\n";
            $str =~ s/1/X/g;
            $str =~ s/0/\-/g;
            print $str, "\n";
        }
    }


    $todo = (@$structlist > $N_ALIGNMENTS ? $N_ALIGNMENTS : @$structlist);
    print "\n",
"____________ Best detailed alignments       ______________________________\n";
    for (my $i = 0; $i < $todo; $i++) {
        my $idx = $indices [$i];
        my $f = $$structlist[$idx];
        my $structfile = get_path (@$struct_dirs, $f . $bin_suffix);
        my $pair_set = $r::r[$idx][9];
        my $coord = coord_read ($structfile) || bad_exit ("from coord_read");
        my ($score, $sw_cvr, $nw_cvr);
        my $a = $r::r[$idx];

        print
"__________________________________________________________________________\n",
              "Alignment to $$structlist[$idx]\n";
        ($score, $sw_cvr, $nw_cvr) = ($$a[10], $$a[3], $$a[4]);
        printf "z-score is %.4g sw cover: %.2f nw cover %.2f\n",
                $score, $sw_cvr, $nw_cvr;

        print pair_set_pretty_string ($pair_set, $seq, coord_get_seq ($coord))
        }

    if ($N_PIR > 0) {
        if ( -d $pir_dir) {
            print "Writing pir files to existing dir: $pir_dir\n"; }
        else {
            if ( ! mkdir ($pir_dir, 0777)) {
                bad_exit ( "Fail creat pir_dir ($pir_dir): $!"); } } }

#   -------------- write pir format files
    $todo = (@$structlist > $N_PIR ? $N_PIR : @$structlist);
    for (my $i = 0; $i < $todo; $i++) {
        my $idx = $indices [$i];
        my $f = $$structlist[$idx];
        my $structfile = get_path (@$struct_dirs, $f . $bin_suffix);
        my $pair_set = $r::r[$idx][1];   # N & W pair set
        my $coord = coord_read ($structfile) || bad_exit("coord_read");
        my $coord_size = coord_size ($coord);
        my $seq_size   = seq_size ($seq);
        my $blah = pir_from_raw_model ($title, $seq, $f, substr ($f, 0, 4),
                                       $coord, $pair_set);
        if ( ! open (PIR, ">$pir_dir/$f.pir")) {
            print STDERR "Open fail on $pir_dir/$f.pir: $!\n"; bad_exit('fatal') ; }
        print PIR $blah;
        close (PIR);
    }
#   -------------- build models

    $todo = (@$structlist > $N_MODELS ? $N_MODELS : @$structlist);
    if ( -d $modeldir) {
        print "Directory $modeldir exists. Adding new models.\n"; }
    else {
        if ( ! mkdir("$modeldir",0777)) {
            bad_exit( "Fail create modeldir ($modeldir): $!"); } }

    print "Writing models for $todo structures\n";
    for (my $i = 0; $i < $todo; $i++) {
        my $idx = $indices [$i];
        my $f = $$structlist[$idx];
        my $structfile = get_path (@$struct_dirs, $f . $bin_suffix);
        my $pair_set = $r::r[$idx][1];
        my $coord = coord_read ($structfile) || bad_exit("coord_read");
        my $new_model = make_model ($pair_set, $seq, $coord); # Build model
        my $pdbfile = "${f}";
        $pdbfile =~ s/\.bin//;
        $pdbfile = "${modeldir}/${pdbfile}.pdb";
        coord_2_pdb ($pdbfile, $new_model, $seq);       # and write coordinates
        my $a = $r::r[$idx];
        my $score = $$a[10];
        if (defined ($formatflag)) {
            if ($formatflag eq 'meta') {
                polish_fmt ($pdbfile, $title, $i, $coord, $score);
            }
        }
    }
    undef (@r::r);
    return EXIT_SUCCESS;
}

# ----------------------- get_prof   --------------------------------
# Run blast and get a sequence profile
# Return EXIT_FAILURE / EXIT_SUCCESS
sub get_prof ($ $)
{
    my ($seqfile, $chk3) = @_;

    use File::Temp qw/ tempdir tempfile /;
    use lib "$ENV{HOME}/../torda/bin";
    use Blastwrap qw (blastwrap);
    if ( ! defined ($chk3)) {
        print STDERR "Bug get_prof: chk3 not defined\n"; return EXIT_FAILURE;}
    my (undef, $chk1)           = tempfile (OPEN => 0);
    my (undef, $chk2)           = tempfile (OPEN => 0);
    my (undef, $errfile)        = tempfile (OPEN => 0);
    my %opts = (
                chkpt1    => $chk1,
                chkpt2    => $chk2,
                chkpt3    => $chk3,
                e_val1    => 1e-10,
                e_val2    => 1e-8,
                e_val3    => 1e-5,
                iter1     => 3,
                iter2     => 2,
                iter3     => 2,
                keep_hits => '500',
                num_aln   => 1,
                num_brf   => 1,
                num_proc  => 2,
                out1      => '/dev/null',
                out2      => '/dev/null',
                out3      => '/dev/null',
                no_pdb    => 1,
                seg_filt  => 'm S'
                );
    if (blastwrap ($seqfile, $errfile, %opts) == EXIT_FAILURE) {
        return EXIT_FAILURE;}
    foreach my $f ($errfile, $chk1, $chk2) {
        if ( -f $f) {
            if ( ! (unlink ($f))) {
                warn "Fail unlinking errfile $errfile\n";} }  }

    return EXIT_SUCCESS;
}

# ----------------------- get_seq_vec -------------------------------
# Calculate the appropropriate probability vector for a specified
# sequence.
sub get_seq_vec ( $ $ $)
{
    my ($seq, $classfile, $abs_error) = @_;
    my $aa_clssfcn;
    if ( ! ($aa_clssfcn = aa_strct_clssfcn_read ( $classfile, $abs_error ))){
        print STDERR "Failed reading classification from $classfile\n";
        return undef;
    }
    return (aa_2_prob_vec ($seq, $aa_clssfcn));
}

# ----------------------- mymain  -----------------------------------
# Arg 1 is a sequence file. Arg 2 is a structure file.
sub mymain ()
{
    use Getopt::Std;
    my (%opts);
    my (@struct_list, @struct_dirs);
    my ($seqfile, $libfile);
    my $fatalflag = undef;
    @struct_dirs = @DFLT_STRUCT_DIRS;

    my $tmp_prf_flag;
    if ( !getopts ('a:d:h:i:m:pt:', \%opts)) {
        usage(); }
    if ( defined ( $opts { a })) {$N_ALIGNMENTS = $opts { a }}
    if ( defined ( $opts { d })) {$modeldir     = $opts { d }}
    if ( defined ( $opts { h })) {$N_BRIEF_LIST = $opts { h }}
    if ( defined ( $opts { i })) {$N_PIR        = $opts { i }}
    if ( defined ( $opts { m })) {$N_MODELS     = $opts { m }}
    if ( defined ( $opts { p })) {$tmp_prf_flag = 1}
    if ( defined ( $opts { t })) {push (@struct_dirs,split( ',', $opts { t }))}
    undef %opts;

    set_params ();

    if ( $#ARGV < 0) {
        print STDERR "Must have at least a sequence file\n";
        usage(); }
    if ( $#ARGV < 1) {
        print STDERR "Please give me a structure library / file\n";
        usage();
    }

    $seqfile   = $ARGV[0];
    $libfile   = $ARGV[1];

    check_dirs (@struct_dirs);
    if (@struct_dirs == 0) {
        die "\nNo valid structure directory. Stopping.\n"; }

    (@struct_list  = get_prot_list($libfile)) || $fatalflag++;
    check_files (@struct_dirs, @struct_list, $bin_suffix) && $fatalflag++;
    check_files (@DFLT_VEC_DIRS, @struct_list, $vec_suffix) && $fatalflag++;
    if ($fatalflag) {
        print STDERR "struct dirs were @struct_dirs\n"; }
    if ( ! ( -f "$MATRIX_DIR/$MATRIX_FILE")) {
        print STDERR "Cannot find matrix file $MATRIX_DIR/$MATRIX_FILE\n";
        $fatalflag++;
    }
    if ( ! ( -f "$PARAM_DIR/$RS_PARAM_FILE")) {
        print STDERR "Cannot find param file $PARAM_DIR/$RS_PARAM_FILE\n";
        $fatalflag++;
    }

    if ( ! ( -f "$CLASS_DIR/$CLASS_FILE")) {
        print STDERR "Cannot find classification $CLASS_DIR/$CLASS_FILE\n";
        $fatalflag++;
    }

    if ( $fatalflag) {
        print STDERR "Fatal problems\n";
        return EXIT_FAILURE;
    }

    print
        "Library from $libfile with ", $#struct_list + 1,
        " template structures\n",
        "Brief results printed for $N_BRIEF_LIST templates\n",
        "Long alignments printed for $N_ALIGNMENTS templates\n",
        "Models made for the best $N_MODELS models and put in $modeldir\n",
        "Modeller/pir format files ($N_PIR thereof) put in $pir_dir\n";
    my $top_temp = '.';
    my $tmp_dir  = tempdir (DIR => $top_temp,  CLEANUP => 1);

    my $profile_file;
    if ( ! defined ($tmp_prf_flag)) {
        my $tmpfh;
        ($tmpfh, $profile_file) =
            tempfile (DIR => $tmp_dir, OPEN => 0, SUFFIX => '.prof');
    } else {
        use File::Basename;
        my ($name, $path, $suffix) = fileparse ($seqfile, qr{\..*});
        $profile_file = "${name}.prof";
    }


    if (get_prof ($seqfile, $profile_file) == EXIT_FAILURE) {
        bad_exit ('Failed to build profile'); }
    my $t;
    if ( ! (defined ( $t = blst_chk_read ($profile_file)))) {
        print STDERR "Fail reading seq from $profile_file\n";
        return EXIT_FAILURE;
    }
    my $seq = seqprof_get_seq ($t);
    undef ($t);
    my $seqvec = get_seq_vec ($seq, "$CLASS_DIR/$CLASS_FILE", $ABS_ERROR);
    if (! defined ($seqvec)) {
        print STDERR "Failed making sequence probability vector\n";
        return EXIT_FAILURE;
    }

    print "Sequence read from \'$seqfile\' with length ",
           seq_size ($seq), " residues\n";

    my $title = $seqfile;
    $title =~ s/\..*$//;
    my $formatflag = 'not set';

    my $r = do_lib (@struct_list , @struct_dirs, $profile_file,
                    $title, $seqvec, $formatflag);
    if ($r == EXIT_FAILURE) {
        bad_exit ('calculation broke') }

        print
"__________________________________________________________________________\n",
    "Wurst gegessen at ", scalar (localtime()), "\n";
    my ($user, $system, $crap, $crap2, $host);
    ($user, $system, $crap, $crap2) = times();
    printf "I took %d:%d min user and %.0f:%.0f min sys time\n",
    $user / 60, $user % 60, $system / 60, $system % 60;
    use Sys::Hostname;
    $host = hostname() || {$host = 'no_host'};
    print "Run on $host\n";
    return EXIT_SUCCESS;
}
# ----------------------- main    -----------------------------------
exit (mymain());
