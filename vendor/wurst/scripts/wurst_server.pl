#!/usr/bin/perl
# This is a relatively minimal alignment script for the wurst server.
# This version will attempt to generate a sequence profile by running blast.
# rcsid = $Id: wurst_server.pl,v 1.1 2007/09/28 16:57:00 mmundry Exp $;

=pod

=head1 NAME

wurst_server.pl

=head1 SYNOPSIS

wurst_server.pl [-a n]
S<[-h n_prot_hit_list]>
S<[-d n_detailed_alignments]>
S<[-f fmt]>
S<[-m n_models]> title email_address sequence

=head1 DESCRIPTION

This is a simple script for the wurst web server. Many values are
hardcoded. The mission is to take

=over

=item *

a sequence

=item *

an email address

=item *

a title for the sequence

=back

and mail back

=over

=item *

a results file

=item *

files with built models

=back

Alignments use

=over

=item S<profile to structure alignments>

=item S<profile to profile alignments>


=back

=head1 OPTIONS

=over

=item -a n

If I<n> is zero, then all results will be sent back as embedded
text.  If I<n> is a number than, of the models built, I<n> will
be sent per mail message.

=item -h n_prot_hit_list

In the terse output, print out this many folds. These are not
detailed alignments.

=item -d n_detail_alignments

In the next section print out I<n_detail_alignments> alignments.

=item -f format

if I<format> is 'meta', then we can process the coordinate files
so they are suitable for the Leszek's meta server.

=item -m n_models

Create I<n_models> models for the sequence and mail them back.

=back

=head1 IMPLEMENTATION NOTES

=over

=item S<temporary directories>

Each script gets a temporary directory to play in. They should be
cleared up by perl's automatic directory cleaning.

=item output

We do print to stderr and stdin, but since this is going to be in
a server, we join the two streams near the start of the
script. When results are mailed back, the whole lot gets mailed
back to the return address.

=item S<hacks for server>

For the livebench server, we have to force the number of models
to 10, otherwise it uses the default of 5. This is done by
testing for the 'meta' flag amongst the options.

=item S<execution dir>

We run from the top level temporary directory.

=item Signals

We trap common signals and send back our stdin/stdout. This is
may well happen on reboots.

=back

=head1 BUGS AND CHANGES

There is a terrible bug on which can be provoked by some
errors. If you bad_exit() tries to call the close_and_mail()
function, this can end up with infinite recursion.

=cut

use vars qw ($LIB_FILE $MATRIX_DIR  $PARAM_DIR $CLASS_DIR $CLASS_FILE
             $RS_PARAM_FILE $FX9_PARAM_FILE );
use FindBin;
do "$FindBin::Bin/paths.inc" || die $@;
if ($@) {
    die "broke reading paths.inc:\n$@"; }

use strict;
use warnings;

BEGIN {
    use vars qw ($wurst_home);
    $wurst_home = '/home/other/wurst';
};

use lib "$FindBin::Bin../blib/arch";
use lib "$FindBin::Bin../blib/lib";
#use lib "$FindBin::Bin/../src/Wurst/blib/lib";
#use lib "$FindBin::Bin/../src/Wurst/blib/arch";

use Wurst;
use POSIX qw(EXIT_SUCCESS EXIT_FAILURE);

# ----------------------- Defaults  ---------------------------------
# These are numbers you might reasonably want to change.
# They should (will) be changeable by options.
use vars qw ($N_BRIEF_LIST $N_MODELS $N_ALIGNMENTS $DFLT_MAX_ATTACH);
$N_BRIEF_LIST   = 50;
$N_MODELS       =  5;
$N_ALIGNMENTS   = 20;
$DFLT_MAX_ATTACH = 5;

# This is a top level temporary directory. We can have a cron job
# run around and clean it up at night. Each job creates its own
# temporary directory under this one.
use vars qw ($log_base $top_temp);
*top_temp   = \"$wurst_home/wurst_delete_able_temp";
# Where we will write logs to
*log_base   = \"$wurst_home/wurst_logs/log";

# Define our mail program, reply-to address and from address.
use vars qw ($mail_from_addr $mail_prog);
*mail_from_addr = \'"Wurst results" <noreply@zbh.uni-hamburg.de>';
*mail_prog      = \'/usr/bin/mailx';

# Switches..
# During testing, we do not want to be able to switch off things
# like the calculation, mailing... These are turned on and off
# here.

use vars qw ($redirect_io $really_mail
             $do_the_calculation $fake_args $verbose);
my $nothing = undef;
*redirect_io  =       \1;
*really_mail  =       \1;
*do_the_calculation = \1;
*fake_args          = \$nothing;
*verbose            = \$nothing;
undef $nothing;

# ----------------------- Global variables  -------------------------
# Unfortunately, we need to store some things here, mainly in
# case we have to quickly die. The bad_exit() routine can mail
# back something informative if it knows the address and job
# title.
use vars qw ($email_address $tmp_dir $title);
# The directory where we put models will live under our $tmp_dir
use vars qw ($modeldir);

# ----------------------- Fixed paths   -----------------------------
use vars qw ($LIB_FILE);
#*LIB_FILE = \"$wurst_home/wurst_server/FoldLibs/debug_list.folds";
*LIB_FILE = \"$wurst_home/wurst_server/FoldLibs/pdb90.list";

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

# Now, some magic constants for rescoring
#                        sw_geog   sw_nsg sw_open sw_widen
# From optimisation run, -1.021    -0.517    -2.809     0.272
use vars qw ( $SW_RS_GEOG $SW_RS_NSG $SW_RS_OPEN $SW_RS_WIDEN $SW_SCR2);
( $SW_RS_GEOG, $SW_RS_NSG, $SW_RS_OPEN, $SW_RS_WIDEN, $SW_SCR2) =
( -0.678,        -0.021,    -2.04,       -0.188,      0.0654);
use vars qw( @DFLT_STRUCT_DIRS  @PROFILE_DIRS @DFLT_VEC_DIRS
             $CLASS_DIR $CLASS_FILE $ABS_ERROR
             $phd_suffix $bin_suffix $prof_suffix $vec_suffix);
*DFLT_STRUCT_DIRS = ["$wurst_home/wurst_server/FoldLibs/pdb90_bin",
                     "$wurst_home/wurst_server/FoldLibs/pdb90_bin.old" ];
*PROFILE_DIRS     = ["$wurst_home/wurst_server/FoldLibs/pdb90_prof",
                     "$wurst_home/wurst_server/FoldLibs/pdb90_prof.old",
                     ];
*DFLT_VEC_DIRS    = ["$wurst_home/wurst_server/FoldLibs/pdb90_vec_5mer8",
                     "$wurst_home/wurst_server/FoldLibs/pdb90_vec_5mer8.old"];
*bin_suffix       = \'.bin';
*prof_suffix      = \'.prof';
*phd_suffix       = \'.phd';
*vec_suffix       = \'.vec';
*MATRIX_FILE      = \'741.out';
*CLASS_DIR        = *PARAM_DIR;
*CLASS_FILE       = \'5mer_seq_angle4/r8.influ-o-data-1';
*ABS_ERROR        = \ 0.4;
# ----------------------- set_params    -----------------------------
# This gets its own function because it can be more complicated
# if, in the future, we have a version depending on various
# options like whether or not we have secondary structure
# information.

sub set_params ()
{
    *sw1_pgap_open  =  \  1.681;
    *sw1_pgap_widen =  \  0.311;
    *sw1_qgap_open  =  \ $sw1_pgap_open;
    *sw1_qgap_widen =  \ $sw1_qgap_widen;
    *wt1            =  \  0.198;
    *wt2            =  \  0.619;
    *wt3            =  \ (1 - ($wt1 + $wt2));
    *off1           =  \  0.573;
    *nw_pgap_open  =   \  $sw1_pgap_open;
    *nw_qgap_open  =   \  $sw1_qgap_open;
    *nw_pgap_widen =   \  $sw1_pgap_widen;
    *nw_qgap_widen =   \  $sw1_qgap_widen;
    *nw_sec_pnlty  =   \  0.0;
}

# ----------------------- log_job   ---------------------------------
# Minimal logging of job.
# We just save the first few characters of the title. It is useful
# for checking jobs from eva/livebench.
# The title gets single quotes, since it is the only thing with a
# totally unpredictable amount of white space.
# The arguments should be
# email_address, 'start' or 'end', title.
sub log_job ($ $ $)
{
    my ($addr, $text, $title) = @_;

    $title =~ s/^ +//;                 # Remove leading white space
    $title = substr ($title, 0, 15);
    $title = "'$title'";
    my ($sec, $min, $hour, $mday, $mon, $year, $wday, $yday, $iddst) =
        localtime (time());
    $mon +=1;
    $year += 1900;
    my $name = "${log_base}_${mday}_${mon}_${year}";  # Our file name
    my $hostname = `uname -n`;
    chomp $hostname;
    my $logtext = "$addr $text ". localtime (time()). " $hostname $title\n";
    if ( ! (open (LOGFILE, ">>$name"))) {
        print STDERR "Failed logging to $name\n"; return; }
    print LOGFILE $logtext;
    close (LOGFILE);
    return 1;
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
# Return zero if all OK. Return a number of missing files if some
# are not found.
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

# ----------------------- mail_file   -------------------------------
# I do not yet know the best way to send mail.
# We really have to use mailx since it is posix standard and mail
# is deprecated. The problems include
# 1. -a for attachments is not standard. Posix does not specify the
#    mechanism
# 2. Currently, our mail uses -q to read from a file. Posix does
#    not support this.
# Because it is so difficult and messy, all mail should be routed
# through here.
# Return 1 on success, undef otherwise.
sub mail_file ($ $ $ \@)
{
    my ($textfile, $subject, $address, $attachments) = @_;
    $ENV{sendwait} = '1';    # Might persade mailer to finish before returning
    $ENV{encoding} = '8bit'; # Otherwise mailer thinks about quoted-printable
                             # This worked with mail/nail, but is not in
                             # the posix description of mailx.

    if ( ! defined ($address)) {
        print STDERR "No address found ! $subject $address\n"; return undef;}

    if ( ! -f ($textfile)) {
        print STDERR ('Cannot find mail file ! $subject $address'); return undef;}

    my @cmdline = $mail_prog;

    if ( defined ($subject)) {
        push (@cmdline, '-s', $subject); }

    foreach my $a (@$attachments) {
        if ( ! -f $a ) {
            print STDERR "model $a disappeared $subject $address\n"; }
        push (@cmdline, '-a', $a);
    }

    if ($textfile) {
        push (@cmdline, '-q', $textfile); }

#   We used to use a -r or -S like
#       push (@cmdline, '-S', "from=$mail_from_addr");
#   but -r works on old mail and -S on newer mailers.
#   Setting the environment works on more. Reading the posix page
#   suggests that one cannot do this portably.
    $ENV{from} = "$mail_from_addr";
    push (@cmdline, $address);
    my $mail_ret = 0;

    if ($really_mail) {
        if (! (open (MAIL, '|-', @cmdline))) {
            print STDERR "Open fail on mail pipe: $!\n";}
        close (MAIL);
        $mail_ret = $?;
    } else {
        open (STDOUT, ">del_me_out");
        open (STDERR, ">del_me_err");
        open (CRAP, "<$textfile");
        while (<CRAP>) {print };
        print " I would invoke @cmdline\n";
        close (STDERR); close (STDOUT);
        $mail_ret = 0;
    }
    wait;
    if ($mail_ret != 0) {
        return undef;}
    else {
        return 1;}
}

# ----------------------- close_up_and_mail -------------------------
# This could be a happy or unhappy exit.
# Grab any files we can find and send them away.
sub close_up_and_mail ($ $ $)
{
    my ($subject, $address, $max_attach) = @_;

#   Send back the model files...
    my @flist = ();
    if ( ! opendir (MODELS, $modeldir)) {
        $max_attach = 0;
    } else {
        @flist = grep (!/^\.$|^\.\.$/, readdir(MODELS));
        closedir (MODELS);
    }
#   Now, we send back our output.
#   From here on, we run the risk of losing error messages.
#   Anything from this point will go to a general error file.
    close (STDERR);
    close (STDOUT);
    open (STDERR, ">>emergency_wurst_error");
    open (STDOUT, ">>emergency_wurst_stdout");
    restore_handlers();
    my @nothing = ();
    if ($max_attach == 0) {
        mail_file ("$tmp_dir/tempout", $subject, $address, @nothing);
        for (my $i = 0; $i < @flist; $i++) {
            my $s = 'Model '. ($i+1) . ' of ' . @flist .
                ' on ' . $flist[$i] . " from $subject";
            mail_file ("$modeldir/$flist[$i]", $s, $address, @nothing);
        }
    } else {
        my $count = 1;
        my $total = int (@flist / $max_attach);
        if (($total * $max_attach) < @flist) {$total++;}
        while (@flist) {
            my @x = splice (@flist, 0, $max_attach);
            map ($_ = "$modeldir/$_", @x);
            my $s = "$subject, part $count of $total";
            if ( !(mail_file ("$tmp_dir/tempout", $s, $address, @x))) {
                my $host = hostname();
                print STDERR "Mail bust. $host";
            }
            $count++;
        }
    }
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
    restore_handlers();  # otherwise an ugly loop is possible
    print STDERR "Error: \"$msg\"\n";
    if (! defined ($title)) {
        $title = 'unknown'; }
    my $subject = "Failed calculating on $title";
    my @attach = undef;
    if (defined ($email_address)) {
        close_up_and_mail ($subject, $email_address, @attach) ;}
    else {
        ;        # Here, we should write to syslog or something.
    }
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

    score_fx_prof ($fx_mat, $profile, $coord, $fx_params);
    $fx_mat = score_mat_scale ($fx_mat, $wt1);
    $fx_mat = score_mat_shift ($fx_mat, $off1);

    score_prof_prof ($sim_scr_mat, $profile, $tmplt_prof, $submat);
    $sim_scr_mat = score_mat_scale ($sim_scr_mat, $wt2);

    score_pvec ($seq_strct_mat, $seqvec, $tmplt_vec);
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

    my $num_scrs = 1000;
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
    my $mod_z = $sw_coverage * $z_scr; # delete this line

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

# ----------------------- parent_name_casp --------------------------
# Given a coordinate name as returned by wurst, convert it to the
# form liked by the livebench server / casp. This is the version that
# uses the information in the file.
sub old_parent_name_casp ($)
{
    my $in = shift;
    if (length ($in) <= 4) { # Do not mess with the name if it is too short
        return $in; };
    my $pdbname = lc(substr ($in, 0, 4));
    if (length ($in) > 4) {
        $pdbname = "$pdbname" . '_' . substr($in, 4, 1); }
    return $pdbname;
}
# Now is the (temporary) version which will use the file name.
# Our .bin files have temporarily got broken chain identifiers.
sub parent_name_casp ($)
{
    my $in = shift;
    use File::Basename;
    my ($name, $path, $suffix) = fileparse ($in, '.bin');
    if (length($in) > 4) {
        $name = substr($name, 0, 4) . '_' . substr($name, 4, 1);}
    return $name;
}

# ----------------------- polish_fmt --------------------------------
# For the livebench server, we have to add some information.
# We get the filename to fix, the title/identifier, model number and
# parent coordinates.
sub polish_fmt ($ $ $ $ $) {
    my ($fname, $title, $m_num, $prnt_coord, $score) = @_;
    my $header =
'PFRMAT TS
AUTHOR Huber-Torda-server
';
    my $trailer =
'TER
END
';

#   Build the header

    $header = "${header}TARGET $title\n";
    my $tmp2 = sprintf ('%6g', $score);
    $header = "${header}SCORE $tmp2\n";

    $header = "${header}MODEL " . ($m_num + 1) . "\n";

#   Version based on information in .bin file:
#   my $tmp = coord_name ($prnt_coord);
#   if (defined $tmp) {
#       if (length ($tmp) > 3) {
#           $header = "${header}PARENT " .  parent_name_casp ($tmp). "\n"; } }
    $header = "${header}PARENT " .  parent_name_casp ($fname). "\n";

#   Read up the bare coordinates in one go.
    if (! ( open (ORIG, "<$fname"))) {
        warn "open fail for read on $fname\n";
        return EXIT_FAILURE;
    }
    local $/;        # localise slurp mode
    my $content = <ORIG>;
    close (ORIG);
#   For casp, we should also remove ORIG and SCALE lines.
    $content =~ s/^ORIGX.+\n//gm;
    $content =~ s/^SCALE.+\n//gm;

    if (! (open (MANGLE, ">$fname"))) {
        warn "open fail for write on $fname\n";
        return EXIT_FAILURE;
    }

    print MANGLE "${header}${content}${trailer}";
    close (MANGLE);
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
    my (@pair_sets);

    my $pfile = "$PARAM_DIR/$RS_PARAM_FILE";
    my $rescore_params = param_rs_read ($pfile) || bad_exit( "Rescore params");

    my $matname = "$MATRIX_DIR" . '/' . "$MATRIX_FILE";
    my $submat = sub_mat_read ($matname) || bad_exit( "Fail on $matname");

    my ($tmp_p, $fx_params);
    $tmp_p = "$PARAM_DIR/$FX9_PARAM_FILE";
    $fx_params = param_fx_read ($tmp_p) || bad_exit ( "Fail on $tmp_p");

    my $profile = blst_chk_read ($profile_file);
    my $seq     = seqprof_get_seq ($profile);
    if ( ! $profile) {
        warn "Failed reading profile $profile_file\n"; return EXIT_FAILURE; }
    for (my $i = 0; $i < @$structlist; $i++) {
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
    if ( seq_size ($seq) > 250) {
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

    $todo = (@$structlist > $N_MODELS ? $N_MODELS : @$structlist);
    if ( ! mkdir("$modeldir",0777)) {
            bad_exit( "Fail create modeldir ($modeldir): $!"); }

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
sub get_prof ($ $)
{
    my ($seq, $chk3) = @_;
    use File::Temp qw/ tempdir tempfile /;
    use lib "$wurst_home/bin";
    use Blastwrap qw (blastwrap);
    my (undef, $errfile)        = tempfile (DIR => $tmp_dir, OPEN => 0);
    my (undef, $seqfile)        = tempfile (DIR => $tmp_dir, OPEN => 0);
    if ( ! (open (TSEQ, ">$seqfile"))) {
        print STDERR "Seq writing failure\n"; return EXIT_FAILURE; }
    print TSEQ $seq;
    close (TSEQ);
    my %opts = (
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
                seg_filt  => 'm S',
                verbosity => 0
                );

    if (blastwrap ($seqfile, $errfile, %opts) == EXIT_FAILURE) {
        return EXIT_FAILURE;}
    foreach my $f ($errfile, $seqfile) {
        if ( -f $f) {
            if ( ! (unlink ($f))) {
                warn "Fail unlinking errfile $errfile\n";}
        }
    }
    return EXIT_SUCCESS;
}

# ----------------------- break_string ------------------------------
# Given a long input string, like a sequence, return a new string
# which is like a readable sequence format.
sub break_string ($)
{
    my $in = shift;
    my $out;

    my @words = split (/\n/, $in);
    if ($words[0] =~ m/^\>/) {
        $out = "$words[0]\n";
        $in =~ s/$words[0]//;
        $in =~ s/\s//g;
    }

    my $n = 0;
    while (my $s = substr ($in, 0, 10)) {
        $out .= $s;
        substr ($in, 0, 10) = '';
        if ($n++ == 5) {
            $out .= "\n";
            $n = 0;
         } else {
            $out .= " ";
        }
    }

    return $out;
}
# ----------------------- catch_kill     ----------------------------
# The main thing is, if we get a KILL or TERM, to call exit and get
# out of here. This means there is a better chance of closing files
# wherever we were up to.
sub catch_kill
{
    my ($sig) = @_;
    bad_exit ("signal $sig received");
}

# ----------------------- kill_handlers  ----------------------------
# set up signal catchers so we can call exit() and die gracefully.
sub kill_handlers ()
{
    $SIG{INT } = \&catch_kill;
    $SIG{QUIT} = \&catch_kill;
    $SIG{TERM} = \&catch_kill;
}

# ----------------------- restore_handlers --------------------------
# If we are at the stage of mailing, we no longer want to trap
# interrupts. Otherwise, they will call the bad_exit routine again.
sub restore_handlers ()
{
    $SIG{INT } = 'DEFAULT';
    $SIG{QUIT} = 'DEFAULT';
    $SIG{TERM} = 'DEFAULT';
}

# ----------------------- x_to_m  -----------------------------------
# Sometimes we get sequences with an 'x' or many of them.
# Strictly, this means unknown. In practice, they are most often
# methionines, changed to selenomethionine.
# Mainly for the sake of the competition, lets change them here to
# met (m) residues.

sub x_to_m ($)
{
    my $in = shift;
    my $out = '';
    my @words = split (/\n/, $in);
    if ($words[0] =~ m/^\>/) {     # If we have a comment line,
        $out = "$words[0]";        # do not change it
        shift (@words);            # but remove it
    }
    foreach my $t (@words) {
        if (length ($t) == 0) {
            next;}
        $t =~ s/x/m/gi;
        $out = "${out}\n$t";
    }
    return $out;
}

# ----------------------- reduce_priority ---------------------------
# We can reduce our own priority to be sociable.
sub reduce_priority ()
{
    my $PRIO_PGRP = 1;
    my $low_priority = 2;      # 0 is normal. 10 is very low.
    setpriority ($PRIO_PGRP, getpgrp (0), $low_priority);
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
    if ($fake_args) {
        @ARGV = ('-h', '60', '-a', '2', '-d', '5', '-m', '1', '-f', 'meta',
                 'testing script only. Fake args',
                 'tordatest@zbh.uni-hamburg.de',
                 'a c d e l k w c e d e d w l v i l v w l w e d a c d e
                  l k w c e d e d w l v i l v');
        print "using fake args\n";
    }

    my $horror_debug = undef();
    if ($horror_debug) {
        open (DBG, ">/tmp/debugme") || die "splat";
        for (my $i = 0; $i < @ARGV; $i++) {
            print DBG "$ARGV[$i]\n"; }
        close (DBG);
        exit 0;
    }
    reduce_priority();
    my (%opts);
    my (@struct_list);
    my $seqfile = undef;
    my $fatalflag =  0;
    my @struct_dirs = @DFLT_STRUCT_DIRS;

#   Set up directory for output and models. It is the only place
#   we will be allowed to write to, so cleanup is easy after disaster

    if ( ! -d $top_temp ) {
        mkdir ($top_temp) || bad_exit ("Fail creating $top_temp: $!"); }
    if ( ! chdir ($top_temp)) {
        bad_exit ("Failed to cd to $top_temp: $!"); }
    use File::Temp qw (tempdir);
    $tmp_dir  = tempdir (DIR => $top_temp,  CLEANUP => 1);
    chmod 0777, $tmp_dir; # So normal people can look for disasters.
    $modeldir = "$tmp_dir/modeldir"; # Where models go. Automatically cleaned.
#   If the machines get rebooted or jobs are killed, try to mail
#   back some information. Trap the errors.
    kill_handlers();
    if ( $redirect_io ) {
#       Before doing anything, fix up stderr, stdout.
        open my $oldout, ">&STDOUT" || bad_exit ("dup stdout");
        open my $olderr, ">&STDERR" || bad_exit ("dup stderr");
        open (STDOUT, ">$tmp_dir/tempout")|| bad_exit ("redirect stdout fail");
        open (STDERR, ">&STDOUT")         || bad_exit ("redirect stderr fail");
        close $oldout || bad_exit ("close fail on stdout");
        close $olderr || bad_exit ("close fail on stderr");
#       From here on, any print statements, to stderr or stdout will
#       go to '$tmp_dir/tempout' in the temporary directory.
    }

    my $max_attach = $DFLT_MAX_ATTACH;
    my $formatflag;
    if ( !getopts ('a:d:f:h:m:', \%opts)) {
        bad_exit ('script failure - getopts broke'); }
    if ( defined ( $opts { a })) {$max_attach   = $opts { a }};
    if ( defined ( $opts { d })) {$N_ALIGNMENTS = $opts { d }}
    if ( defined ( $opts { f })) {$formatflag   = $opts { f }}
    if ( defined ( $opts { h })) {$N_BRIEF_LIST = $opts { h }}
    if ( defined ( $opts { m })) {$N_MODELS     = $opts { m }}
    undef %opts;

    if ($#ARGV != 2) {
        bad_exit ("scripting bug. \$\#ARGV was $#ARGV. Should be 2"); }

    $title           = $ARGV[0];
    $email_address   = $ARGV[1];
    my $sequence     = $ARGV[2];
    log_job ($email_address, 'start', $title);

    if ( lc ($formatflag) =~ m/meta/ ) {
        $formatflag = 'meta'; }   # This means we are running for the polish
    else {                        # meta server
        undef ($formatflag) ; }

    set_params ();
    check_dirs (@struct_dirs);
    check_dirs (@DFLT_VEC_DIRS);

    if (@struct_dirs == 0) {
        bad_exit( "\nNo valid structure directory. Stopping."); }
    @struct_list  = get_prot_list($LIB_FILE);
    if ($#struct_list <= 0) {
        $fatalflag++ ; }
    if (check_files (@struct_dirs, @struct_list, $bin_suffix) != 0) {
        $fatalflag++; }
    if (check_files (@DFLT_VEC_DIRS, @struct_list, $vec_suffix) != 0) {
        $fatalflag++; }

    if ($fatalflag) {
        print STDERR "struct dirs were @struct_dirs\n"; }

    if ( $fatalflag) {
        bad_exit ('Fatal problems'); }

    $sequence =~ s/\015/\n/;
    $sequence = x_to_m ($sequence);
    my $seq = seq_from_string ($sequence);
    if ( ! $seq ) {
        print STDERR "Problem reading sequence:\n$sequence\n";
        bad_exit ('Seq reading failure');
    }

#   Hack for bioinfo.pl server
    if ( defined ($formatflag)) {
        if ($formatflag eq 'meta') {
            if ($N_MODELS < 10) {
                $N_MODELS = 10; }
            if ($max_attach < 10) {
                $max_attach = 10; }
        }
    }

    use File::Basename;
    my $temp = fileparse ($LIB_FILE);
    print
        "Title: \"$title\"\n",
        "Results sent to \"$email_address\"\n",
        "Library from $temp with ", $#struct_list + 1,
        " template structures.\n",
        "Brief results printed for $N_BRIEF_LIST templates.\n",
        "Long alignments printed for $N_ALIGNMENTS templates.\n",
        "Models made for the best $N_MODELS models.\n";
    if ($max_attach) {
        print "Models will be sent as $max_attach attachments per file\n";}
    else {
        print "The model coordinates will be embedded in the mail text\n"; }
    print
        "Sequence length ", seq_size ($seq), " is\n",
        break_string (seq_print ($seq)), "\n";
#   Get the profile in a temporary file.
    my ($tmpfh, $profile_file) =
        tempfile (DIR => $tmp_dir, OPEN => 0, SUFFIX => '.prof');

    if (get_prof ($sequence, $profile_file) == EXIT_FAILURE) {
        bad_exit ('Failed to build profile'); }

    my $seqvec = get_seq_vec ($seq, "$CLASS_DIR/$CLASS_FILE", $ABS_ERROR);
    if (! defined ($seqvec)) {
        print STDERR "Failed making sequence probability vector\n";
        return EXIT_FAILURE;
    }

    if ($do_the_calculation) {
        my $r = do_lib (@struct_list , @struct_dirs, $profile_file,
                        $title, $seqvec, $formatflag);
        if ($r == EXIT_FAILURE) {
            bad_exit ('calculation broke') }
    }

    unlink ($profile_file);
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
    STDOUT->autoflush(1);
    my $subject = "Wurst calculation on $title";
    close_up_and_mail ($subject, $email_address, $max_attach);
    log_job ($email_address, 'stop', $title);
    sleep 25; # Yuk, but I am so worried about the mailer.
    if ( -f "$tmp_dir/tempout" ) {
        unlink ("$tmp_dir/tempout"); }
    return EXIT_SUCCESS;
}

# ----------------------- main    -----------------------------------
exit (mymain());
