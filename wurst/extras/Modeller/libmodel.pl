#!/usr/local/bin/perl
# 18 March 2002
# Most of the other scripts here are for testing with
# deliberately crazy values.
# This one uses the installed, working binary and serious values
# for all parameters.
# rcsid = $Id: libmodel.pl,v 1.1 2007/09/28 16:57:17 mmundry Exp $

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

If you want to run without secondary structure, it is not enough
to omit the filename. Instead, use the B<-s> option described
below.

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

=item B<-m> I<N>

After alignments, I<N> models will be built, otherwise, a small
default number will appear. Set I<N> to zero if you do not want
any models.

=item B<-s>

Do not use secondary structure predictions. This will cause a
different set of parameters to be used in the calculation.

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

The selection of which scores to print out is a bit arbitrary.

=item *

the coverage picture is very ugly. It could be
beautified.

=item *

The coverage picture corresponds to the Smith and
Waterman. Perhaps it should be the Needleman and
Wunsch. Obviously, both are possible, but just a bit ugly.

=cut

use vars qw ($MATRIX_DIR $MATRIX_FILE $PARAM_DIR
             $RS_PARAM_FILE $FX9_PARAM_FILE );
use FindBin;
do "$FindBin::Bin/paths.inc" || die $@;
if ($@) {
    die "broke reading paths.inc:\n$@"; }

use strict;

#use lib "$ENV{HOME}/../torda/pl/lib";  # Where wurst lives after installation
use lib "$FindBin::Bin/../../src/Wurst/blib/lib";
use lib "$FindBin::Bin/../../src/Wurst/blib/arch";
#local paths
use lib "$FindBin::Bin/Wurst/Modeller/blib/lib";
use lib "$FindBin::Bin/Wurst/Modeller/blib/arch";

use Wurst;
use Wurst::Modeller qw(Modeller_model);
use POSIX qw(EXIT_SUCCESS EXIT_FAILURE);


# ----------------------- Defaults  ---------------------------------
# These are numbers you might reasonably want to change.
# They should (will) be changeable by options.
use vars qw ($N_BRIEF_LIST $N_MODELS $N_ALIGNMENTS);
$N_BRIEF_LIST   = 100;
$N_MODELS       =  10;
$N_ALIGNMENTS   =  20;
use vars qw ($modeldir $DFLT_MODELDIR);
*DFLT_MODELDIR = \ 'modeldir';
$modeldir = $DFLT_MODELDIR;


# ----------------------- Constants ---------------------------------
# These are declared globally, and set by set_params().
use vars qw (
             $align_type
             $sw1_pgap_open
             $sw1_qgap_open
             $sim_scale
             $sec_s_scale
             $sw1_pgap_widen
             $sw1_qgap_widen
             $fx_mat_shift
             $sim_mat_bottom
             $sw1_sec_pnlty
             );
use vars qw (
             $sw2_pgap_open
             $sw2_qgap_open
             $sw2_pgap_widen
             $sw2_qgap_widen
             $sw2_sec_pnlty
);


# These parameters will be used for extending alignments via a
# Needleman and Wunsch

use vars qw (
             $nw_pgap_open
             $nw_qgap_open
             $nw_pgap_widen
             $nw_qgap_widen
             $nw_sec_pnlty
);

# Now, some magic constants for rescoring
use vars qw (
             $NW_RS_SCR2 $NW_GAP_GEO $NW_SEQ_GAP $NW_STR_GAP $NW_STR_WDN
             $SW_RS_SCR2 $SW_GAP_GEO $SW_SEQ_GAP $SW_STR_GAP $SW_STR_WDN
             );
( $NW_RS_SCR2, $NW_GAP_GEO, $NW_SEQ_GAP, $NW_STR_GAP, $NW_STR_WDN) =
( 0.41587,     -3.27610,    -1.02167 ,   -11.52223,   -3.08488);
( $SW_RS_SCR2, $SW_GAP_GEO, $SW_SEQ_GAP, $SW_STR_GAP, $SW_STR_WDN) =
( 0.47417,     -5.16621,    -1.35113,    -16.20915,  -6.57821);

use vars qw( @DFLT_STRUCT_DIRS $phd_suffix $bin_suffix);
#*DFLT_STRUCT_DIRS = ['/rsc/appenzeller/data1/pdblib6apr00',
#                     '/rsc/appenzeller/data1/pdb_bin_pool'];
#*DFLT_STRUCT_DIRS = ['/data/procter/Gorgonzola/Clustering/CulPdbSet',
#                     '/data/procter/wurstworld/36Pci686/pdbinlib',
#                     '/bm/hdb5/data1/pdb_bin_pool']; #looks at old files last
*DFLT_STRUCT_DIRS = [ '/bm/wurst_server/foldlibs/EvaNRAug2003',
                      # Eva lib is built specially (some non-ss again)
                      '/bm/pdb_bin_pool' ]; # New recent repository

*bin_suffix       = \'.bin';
*phd_suffix       = \'.phd';

# ----------------------- set_params    -----------------------------
# We put the parameter setting in its own function since we might
# be using one set most of the time, but another set if we have no
# secondary structure info.

sub set_params ($)
{
    my $do_sec_flag = shift;
    if ( $do_sec_flag) {     # Normal parameters
           *sw1_pgap_open  =  \ 8.77616;
           *sw1_qgap_open  = \ 14.7824;
           *sim_scale      =  \ 1.07201;
           *sec_s_scale    =  \ 1.92943;
           *sw1_pgap_widen =  \ 1.12158;
           *sw1_qgap_widen =  \ 2.02476;
           *fx_mat_shift   =  \ 1.21217;
           *sim_mat_bottom = \ -4.33745;
           *sw1_sec_pnlty  =  \ 2.75143;

           *sw2_pgap_open  =  \ 9.78565;
           *sw2_qgap_open  = \ 16.0966;
           *sw2_pgap_widen =  \ 1.11091;
           *sw2_qgap_widen =  \ 1.84835;
           *sw2_sec_pnlty  =  \ 2.51036;

           *nw_pgap_open  =  \ 9.0596;
           *nw_qgap_open  = \ 14.5434;
           *nw_pgap_widen =  \ 1.00474;
           *nw_qgap_widen =  \ 2.01584;
           *nw_sec_pnlty  =  \ 2.56828;

    } else {                # Parameters when we have no sec struct guesses

        *sw1_pgap_open  =  \  9.67385398757248;
        *sw1_qgap_open  =  \  18.9017440741075;
        *sim_scale      =  \  1.00407642907837;
        *sec_s_scale    =  \  0.0;
        *sw1_pgap_widen =  \  0.743565089047655;
        *sw1_qgap_widen =  \  1.6537099328616;
        *fx_mat_shift   =  \  0.20049489028982;
        *sim_mat_bottom =  \ -2.94389560724283;
        *sw1_sec_pnlty  =  \  4.90746866278873;

        *sw2_pgap_open  =  \  8.34897;
        *sw2_qgap_open  =  \ 18.9417;
        *sw2_pgap_widen =  \  0.846406;
        *sw2_qgap_widen =  \  1.75334;
        *sw2_sec_pnlty  =  \  2.58267;

        *nw_pgap_open  =   \  9.85439;
        *nw_qgap_open  =  \  19.7665;
        *nw_pgap_widen =   \  0.80877;
        *nw_qgap_widen =   \  1.48625;
        *nw_sec_pnlty  =   \  2.65248;
    }
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
    print STDERR "Usage: $0 seq_file struct_file phd_file\n";
    exit (EXIT_FAILURE);
}

# ----------------------- get_scores --------------------------------
sub get_scores ($ $ $ $ $)
{
    my ($pair_set, $coord, $seq, $rescore_params, $to_use) = @_;

    my ($scr_tot, $coverage, $score1, $score2, $geo_gap);
    my ($str1, $crap);
    my ($open_cost, $widen_cost, $nseq_gap);
    ($crap, $score1) = pair_set_score( $pair_set );
    ($str1, $crap) =
        pair_set_coverage ($pair_set, seq_size ($seq), coord_size ($coord));
    $coverage = ($str1 =~ tr/1//);       # This is coverage as an integer
    $coverage = $coverage / seq_size ($seq); #and as fraction of sequence

    my ($k_scr2, $k_gap_geo, $k_seq_gap, $k_str_gap, $k_str_wdn);
    if ($to_use eq 's_and_w') {
        ( $k_scr2,     $k_gap_geo,  $k_seq_gap,  $k_str_gap,  $k_str_wdn) =
        ( $SW_RS_SCR2, $SW_GAP_GEO, $SW_SEQ_GAP, $SW_STR_GAP, $SW_STR_WDN);
    } else {
        ( $k_scr2,     $k_gap_geo,  $k_seq_gap,  $k_str_gap,  $k_str_wdn) =
        ( $NW_RS_SCR2, $NW_GAP_GEO, $NW_SEQ_GAP, $NW_STR_GAP, $NW_STR_WDN);
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
               $k_scr2    * $score2 +
               $k_gap_geo * $geo_gap +
               $k_seq_gap * $nseq_gap +
               $k_str_gap * $open_cost +
               $k_str_wdn * $widen_cost;
    return ( $scr_tot, $coverage, $score1, $score2, $nseq_gap, $open_cost);
}

# ----------------------- do_align ----------------------------------
sub do_align ($ $ $ $ $ $)
{
    my ($seq, $template, $sec_s_data, $fx_params, $submat, $rescore_params)
        = @_;

    my $coord  = coord_read ($template) || die "Fail on $template";

#   Now we have all files we need. Start scoring.
#   Begin by giving us four empty score matrices.
    my $sim_scr_mat =                                 # For similarity matrix
        score_mat_new (seq_size ($seq), seq_size (coord_get_seq ($coord)));
    my $sec_scr_mat =                                 # For sec struct matrix
	score_mat_new (seq_size ($seq), seq_size (coord_get_seq ($coord)));
    my $fx_scr_mat =                                  # For fx9 func scores
        score_mat_new (seq_size ($seq), seq_size (coord_get_seq ($coord)));
    my $total_scr_mat =                               # Where we put totals
        score_mat_new (seq_size ($seq), seq_size (coord_get_seq ($coord)));
    score_smat ($sim_scr_mat, $seq, coord_get_seq ($coord), $submat);
    if ($sec_s_scale != 0.0) {
	score_sec  ($sec_scr_mat, $sec_s_data, $coord); }
    score_fx   ($fx_scr_mat, $seq, $coord, $fx_params);
    score_mat_shift ($fx_scr_mat, $fx_mat_shift);
#   We have three score, matrices. Now add them together.
    $total_scr_mat = score_mat_add ($fx_scr_mat, $sim_scr_mat, $sim_scale);
    if ($sec_s_scale != 0.0) {
        $total_scr_mat =
            score_mat_add ($total_scr_mat, $sec_scr_mat,$sec_s_scale); }
#   This actually does the alignment.
    my ($sw_scr_tot, $sw_coverage, $sw_score1, $sw_score2,
        $sw_seq_gap, $sw_strct_gap);
    my ($nw_scr_tot, $nw_coverage, $nw_score1, $nw_score2,
        $nw_seq_gap, $nw_strct_gap);

    my $sw_pair_set =
        score_mat_sum_sec (my $result_mat, $total_scr_mat,
                           $coord, $sw1_sec_pnlty,
                           $sw1_pgap_open, $sw1_pgap_widen,
                           $sw1_qgap_open, $sw1_qgap_widen,
                           $S_AND_W);
    $sw_pair_set =
        score_mat_sum_sec (my $result_mat2, $total_scr_mat,
                           $coord, $sw2_sec_pnlty,
                           $sw2_pgap_open, $sw2_pgap_widen,
                           $sw2_qgap_open, $sw2_qgap_widen,
                           $S_AND_W, $sw_pair_set);

    ($sw_scr_tot, $sw_coverage, $sw_score1, $sw_score2,
        $sw_seq_gap, $sw_strct_gap) =
    get_scores ($sw_pair_set, $coord, $seq, $rescore_params, 's_and_w');

    my $nw_pair_set =
        score_mat_sum_sec (my $result_mat3, $total_scr_mat,
                           $coord, $nw_sec_pnlty,
                           $nw_pgap_open, $nw_pgap_widen,
                           $nw_qgap_open, $nw_qgap_widen,
                           $N_AND_W, $sw_pair_set);
    ($nw_scr_tot, $nw_coverage, $nw_score1, $nw_score2,
     $nw_seq_gap, $nw_strct_gap) =
    get_scores ($nw_pair_set, $coord, $seq, $rescore_params, 'n_and_w');
    my @r =
        ($nw_scr_tot, $nw_pair_set, $sw_scr_tot,
         $sw_coverage, $nw_coverage, $sw_score1, $sw_score2,
         $nw_score1, $nw_score2, $sw_pair_set);
    return (\@r);
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
sub do_lib (\@ \@ $ $)
{
    my ($structlist, $struct_dirs, $seqfile, $phdfile) = @_;
    my (@pair_sets);

    my $pfile = "$PARAM_DIR/$RS_PARAM_FILE";
    my $rescore_params = param_rs_read ($pfile) || die "Rescore params";

    my $matname = "$MATRIX_DIR" . '/' . "$MATRIX_FILE";
    my $submat = sub_mat_read ($matname) || die "Fail on $matname";
    $submat = sub_mat_shift ($submat, $sim_mat_bottom);


    my ($tmp_p, $fx_params);
    $tmp_p = "$PARAM_DIR/$FX9_PARAM_FILE";
    $fx_params = param_fx_read ($tmp_p) || die "Fail on $tmp_p";
    my $sec_s_data;
    if ($sec_s_scale != 0.0) {
        $sec_s_data = sec_s_data_read ($phdfile) || die "Fail on $phdfile";
    } else {
	$sec_s_data = "nothing here"; }

    my $seq    = seq_read ($seqfile) || die "Fail on $seqfile";
    print "Sequence from $seqfile has ", seq_size( $seq), " residues.\n";
    for (my $i = 0; $i < @$structlist; $i++) {
        my $structfile;
        $structfile = get_path (@$struct_dirs, $$structlist[$i] . $bin_suffix);
        $r::r[$i] = do_align ( $seq, $structfile, $sec_s_data,
                               $fx_params, $submat, $rescore_params);
    }
    print "\n",
"____________ Summary of best templates   _________________________________\n";
    my @indices;
    for ( my $i = 0; $i < @$structlist; $i++) {
        $indices[$i] = $i; }


    @indices = sort {
        $r::r[$b][0] <=> $r::r[$a][0];
    } @indices;

    my $todo = (@$structlist > $N_BRIEF_LIST ? $N_BRIEF_LIST : @$structlist);
    printf "%8s %8s %8s %6s %6s %8s %8s %8s %8s\n",
    'struct', 'nw scr', 'sw scr', 'sw cvr', 'nw cvr', 'sw1','sw2','nw1', 'nw2';
    for (my $i = 0; $i < $todo; $i++) {
        my $idx = $indices [$i];
        my $a = $r::r[$idx];
        printf "%8s %8.2g %8.2g %6.2f %6.2f %8.2g %8.2g %8.2g %8.2g\n",
               $$structlist[$idx],
               $$a[0], $$a[2], $$a[3], $$a[4],
               $$a[5], $$a[6], $$a[7], $$a[8];
    }
    print "\n";


#   Now, the somewhat experimental coverage printout
    print
"____________ Summary of coverage of query sequence _______________________\n";
    for (my $i = 0; $i < $todo; $i++) {
        my $idx = $indices [$i];
        my $a = $r::r[$idx];
        my $csize = coord_size (coord_read (
                  get_path (@$struct_dirs, $$structlist[$idx] . $bin_suffix)));
        my ($str, $crap) = pair_set_coverage ($$a[9], seq_size ($seq), $csize);
        print "S & W coverage with $$structlist[$idx]\n";
        $str =~ s/1/X/g;
        $str =~ s/0/\-/g;
        print $str, "\n";
    }


    $todo = (@$structlist > $N_ALIGNMENTS ? $N_ALIGNMENTS : @$structlist);
    print "\n",
"____________ Best detailed alignments       ______________________________\n";
    for (my $i = 0; $i < $todo; $i++) {
        my $idx = $indices [$i];
        my $f = $$structlist[$idx];
        my $structfile = get_path (@$struct_dirs, $f . $bin_suffix);
        my $pair_set = $r::r[$idx][1];
        my $coord = coord_read ($structfile) || die;
        my ($score, $sw_cvr, $nw_cvr);
        my $a = $r::r[$idx];

        print
"__________________________________________________________________________\n",
              "Alignment to $$structlist[$idx]\n";
        ($score, $sw_cvr, $nw_cvr) = ($$a[0], $$a[3], $$a[4]);
        printf "Score is %.4g sw cover: %.2f nw cover %.2f\n",
                $score, $sw_cvr, $nw_cvr;


        if ($sec_s_scale != 0.0) {
            print pair_set_pretty_string ($pair_set, $seq,
                                          coord_get_seq ($coord),
                                          $sec_s_data, $coord);
        } else {
            print pair_set_pretty_string ($pair_set, $seq,
                                          coord_get_seq ($coord));
        }
    }


    $todo = (@$structlist > $N_MODELS ? $N_MODELS : @$structlist);
    if ( ! (-d "$modeldir")) {
        mkdir("$modeldir",0777) ||die "Fail create modeldir ($modeldir): $!"; }

    print "Writing models for $todo structures in $modeldir\n";
    my $made_models = 0;
    for (my $i = 0; $i < $todo; $i++) {
        my $idx = $indices [$i];
        my $f = $$structlist[$idx];
        my $structfile = get_path (@$struct_dirs, $f . $bin_suffix);
        my $pair_set = $r::r[$idx][1];
        my $coord = coord_read ($structfile) || die;
        my $new_model = make_model ($pair_set, $seq, $coord); # Build model
        my $pdbfile = "${f}";
        $pdbfile =~ s/\.bin//;
        $pdbfile = "${modeldir}/${pdbfile}.pdb";
        coord_2_pdb ($pdbfile, $new_model, $seq);       # and write coordinates
        if ($made_models<5) {
            my @gm = Modeller_model("seq",$seq, "$f", $pair_set, $new_model, $coord, "$f",$pdbfile);
            print "Made these models : \n".(join ",",@gm)."\n";
            $made_models++;
        }
    }


    undef (@r::r);
}

# ----------------------- mymain  -----------------------------------
# Arg 1 is a sequence file. Arg 2 is a structure file.
sub mymain ()
{
    use Getopt::Std;
    my (%opts);
    my (@struct_list, @struct_dirs);
    my ($seqfile, $libfile, $phdfile);
    my ($fatalflag, $do_sec_flag) = (0, 1);
    @struct_dirs = @DFLT_STRUCT_DIRS;

    if ( !getopts ('a:d:h:m:st:', \%opts)) {
        usage(); }
    if ( defined ( $opts { a })) {$N_ALIGNMENTS = $opts { a }}
    if ( defined ( $opts { d })) {$modeldir     = $opts { d }}
    if ( defined ( $opts { h })) {$N_BRIEF_LIST = $opts { h }}
    if ( defined ( $opts { m })) {$N_MODELS     = $opts { m }}
    if ( defined ( $opts { s })) {$do_sec_flag = 0}
    if ( defined ( $opts { t })) {push (@struct_dirs,split( ',', $opts { t }))}
    undef %opts;

    set_params ($do_sec_flag);

    if ( $#ARGV < 0) {
        print STDERR "Must have at least a sequence file\n";
        usage(); }
    if ( $#ARGV < 1) {
        print STDERR "Please give me a structure library / file\n";
        usage();
    }

    $seqfile   = $ARGV[0];
    $libfile   = $ARGV[1];


    if ( $#ARGV < 2) {
        use File::Basename;
        my ($n, $d, $t) = fileparse ($seqfile, '\..*');
        $phdfile = $d . $n . $phd_suffix;
    } else {
        $phdfile   = $ARGV[2];
    }
    check_dirs (@struct_dirs);
    if (@struct_dirs == 0) {
        die "\nNo valid structure directory. Stopping.\n"; }

    (@struct_list  = get_prot_list($libfile)) || $fatalflag++;
    check_files (@struct_dirs, @struct_list, $bin_suffix) && $fatalflag++;

    if (( ! -f $phdfile ) && ($do_sec_flag)) {
        print STDERR "No phdfile found in \"$phdfile\"\n";
        $fatalflag++;
    }
    if ( $fatalflag) {
        print STDERR "Fatal problems\n";
        return EXIT_FAILURE;
    }

    print
        "Sequence in $seqfile\n",
        "Library from $libfile with ", $#struct_list + 1,
        " template structures\n",
        "Brief results printed for $N_BRIEF_LIST templates\n",
        "Long alignments printed for $N_ALIGNMENTS templates\n",
        "Models made for the best $N_MODELS models and put in $modeldir\n";

    if ($do_sec_flag) {
        print "Secondary structure info in $phdfile\n"; }
    else {
        print "*Not* using secondary structure predictions\n"; }

    do_lib (@struct_list , @struct_dirs, $seqfile, $phdfile);

        print
"__________________________________________________________________________\n",
    "$0 Ganz fertig at ", scalar (localtime()), "\n";
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
