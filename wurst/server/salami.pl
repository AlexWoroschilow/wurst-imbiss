#!/usr/bin/perl
##$ -clear
#$ -w e
##$ -l arch=glinux -l short=0
#$ -p -50
#$ -S /home/torda/bin/perl
#$ -cwd
#$ -j y
#$ -q stud.q 

=pod

=head1 NAME

salami_2.pl - Given a structure, align it to a library of templates

=head1 SYNOPSIS

libsrch.pl [options]


=head1 DESCRIPTION

Given a structure, align it to every member of a library of
templates. Rank the matches.

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

=item B<-s> I<N>

Number of models to be built. Check if this switch is actually used.
It does set a variable, but we have to check if that line of code is
ever used.

=item B<-e> I<address>

An email address to which results will be mailed. If I<address> has
the value, I<Anonymous>, no mail will be sent. This is for the web
page interface, which must function without requiring an email
address.

=item B<-n> I<title>

This title will be shown in the results page and in the mail which
contains the link to the results web page.

=item B<-q> I<structfile>

This is the path to the query structure file. It is not really an
option  since the script will stop if it is not present.

=item B<-l> I<libfile>

This contains a list of structures (the library) which will be
searched. They are PDB identifiers.

=item B<-t> I<additional_dirs>

A list of additional directories where one may look for
structures. For the web page, this is unlikely to be used. For running
from the command line, it might be useful.

=item B<-r> I<rmsd>

When results are superimposed, pairs of atoms will be removed from the
alignment until the calculated rmsd is less than I<rmsd>.

=item B<-f> I<minFracDME>

When structures are aligned, the frac DME is calculated. Only Values
better than I<minFracDME> will be printed out.

=back

=head1 Running

When called from the web page, the script is typically faced with a list like this

 perl ./salami_2.pl -o 5ptiA_t0Hn0c -n SalamiSearjunkch -e junk@mailinator.com -q /home/other/wurst/salami_jobs/5ptiA.pdb -l /smallfiles/public/no_backup/bm/pdb90.list -s 0 -a 0 -r 3 -f 0.75 -i 100 -x no -p str

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
F<modeldir>.

=head1 THREADING

This script has to read 10000s of coordinate files which can be very
slow, so it uses threads to fetch them in the background and copy them
to a directory called C<$LOCAL_STORAGE> which is currently set to
/dev/shm. This is a ram file system. After working on each file, it is
deleted. The number of file to prefetch is controlled by C<$n_to_buf>
which seems to work fine with values like 5 or 15. One could set it
higher, but we do not want to risk filling the temporary storage.

=head2 Coding style of threads.

The rules are

=over

=item * Threads should be started before any objects are created at
the C level. 

=item * Threads are cleaned up just before exiting.

=back

These are perl threads, not linux/posix threads and this leads to the
problem we have often seen. If a thread exits and there are objects
from the C level, the perl interpreter sees their reference count go
to zero and calls their destructors which calls C<free()>, even when
another thread is using them. This is why it is essential that threads
be created before wurst code is called.

At the moment, there are several running for most of the time.  There
is a check for .bin files at the start which also runs in the
background, but for most of the time one sees

=over

=item * the initial thread

=item * the main reader thread

=item * typically two helper reader threads

=back

The main reader thread pulls protein names from a queue of files to
read ($prot_to_read). For each protein, there is at least one
coordinate file and one or two .vec files to read. For one protein,
these names are given to helper threads (low_reader) which copy them
to local storage. When the .vec and .bin files are copied, the name of
the protein is put in the prot_ready_q queue.

At the same time, the main loop pulls names from the prot_ready_q and
calculates o the protein. When it is finished, it removes the file
from the temporary local storage.

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

=head1 TODO

When we have finished the loop over the library, we could tell the
queues that they are not needed any more. In perl 5.16, this does not
seem to be implemented, but in 5.20, we could do the following:

    $prot_to_read_q->end();  # These function exist in newer perls,
    $prot_ready_q->end();    # but not the current one.
    for (my $i = 0; $i <= $#low_reader_q; $i++) {
        $low_reader_q[$i]->end();
        $low_reader_done_q[$i]->end();
    }

=cut

use FindBin;

use lib "/home/work/other/salamiServer/v02";
use Salamisrvini;
use lib $LIB_LIB;  #initialize in local Salamisrvini.pm;
use lib $LIB_ARCH; #initialize in local Salamisrvini.pm;
use lib "$FindBin::Bin/../pca";
use Similarity::Builder;

use Wurst;

use vars qw ($MATRIX_DIR $PARAM_DIR
  $RS_PARAM_FILE $FX9_PARAM_FILE );

if ($@) {
    die "broke reading paths.inc:\n$@";
}
if ( defined( $ENV{SGE_ROOT} ) ) {
    $MATRIX_DIR = "$ENV{HOME}/../../torda/c/wurst/matrix";
    $PARAM_DIR  = "$ENV{HOME}/../../torda/c/wurst/params";
}

use strict;
use warnings;

use POSIX qw(EXIT_SUCCESS EXIT_FAILURE);

# global variable for using the combined or the single ca/strct functions
our $function;

# ----------------------- Defaults  ---------------------------------
# These are numbers you might reasonably want to change.
# They should (will) be changeable by options.
use vars qw ($N_BRIEF_LIST $N_MODELS $N_ALIGNMENTS $MAG_NUM_ASREST);
$N_BRIEF_LIST = 100;
$N_MODELS     = 5;
$N_ALIGNMENTS = 5;
$MAG_NUM_ASREST = 50;

use vars qw ($modeldir $DFLT_MODELDIR);
*DFLT_MODELDIR = \'modeldir';
$modeldir      = $DFLT_MODELDIR;

# Define our mail program, reply-to address and from address.
use vars qw ($mail_from_addr $mail_reply_to $mail_prog);
*mail_from_addr = \'"Wurst results" <nobody@zbh.uni-hamburg.de>';
*mail_reply_to  = \'nobody@zbh.uni-hamburg.de';
*mail_prog      = \'/usr/bin/mailx';

# Switches..
# During testing, we do not want to be able to switch off things
# like the calculation, mailing... These are turned on and off
# here.

use vars qw ($really_mail $fake_args );
*really_mail  = \1;
*fake_args    = \undef;

# ----------------------- Global variables  -------------------------
# Unfortunately, we need to store some things here, mainly in
# case we have to quickly die. The bad_exit() routine can mail
# back something informative if it knows the address and job
# title.
use vars qw ($email_address $tmp_dir $title);

# ----------------------- Sequence Alignment Constants -----------------------
# These are declared globally, and set by set_params().
# changed by Iryna  
use vars qw (
  $align_type
  $sw1_pgap_open
  $sw1_qgap_open
  $sw1_pgap_widen
  $sw1_qgap_widen
  $s_const
  $s_factor
  $weight
  $gauss_err
  $tau_error
  $ca_dist_error
  $corr_num
  $m_s_scale
  $m_shift
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

# ----------------------- set_params    -----------------------------
# This gets its own function because it can be more complicated
# if, in the future, we have a version depending on various
# options like whether or not we have secondary structure
# information.
# pgap controls penalties in the sequence, qgap in the structure.

sub set_params () {
#    changed by Iryna
#     *sw1_pgap_open  = \3.25;
#     *sw1_pgap_widen = \0.8942;
    *sw1_pgap_open  = \1.915;
    *sw1_pgap_widen = \1.85;
    *sw1_qgap_open  = \$sw1_pgap_open;
    *sw1_qgap_widen = \$sw1_pgap_widen;
    *m_shift        = \-0.1;
    *s_factor       = \-0.003186;
    *s_const        = \-0.03642;
    *weight         = \0.7642;
    *gauss_err      = \0.4;

    *tau_error      = \0.15;
    *ca_dist_error  = \0.385;
    *corr_num       = \4;
    *m_s_scale      = \0;
    *m_shift        = \0;
    *nw_pgap_open  = \$sw1_pgap_open;
    *nw_qgap_open  = \$sw1_qgap_open;
    *nw_pgap_widen = \$sw1_pgap_widen;
    *nw_qgap_widen = \$sw1_qgap_widen;
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
    my $name = "${LOG_BASE}_${mday}_${mon}_${year}";  # Our file name
    my $hostname = hostname();
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
sub get_prot_list ($) {
    my $f = shift;
    my @a;
    if ( !open( F, "<$f" ) ) {
        print STDERR "Open fail on $f: $!\n";
        return undef;
    }
    while ( my $line = <F> ) {
        chomp($line);
        my @words = split( ' ', $line );
        if ( !defined $words[0] ) { next; }
        $line = $words[0];
        $line =~ s/#.*//;     # Toss comments away
        $line =~ s/\..*//;    # Toss filetypes away
        $line =~ s/^ +//;     # Leading and
        $line =~ s/ +$//;     # trailing spaces.
        if ( $line eq '' ) {
            next;
        }
        substr( $line, 0, 4 ) = lc( substr( $line, 0, 4 ) );    # 1AGC2 to 1agc2
        if ( length($line) == 4 ) {    # Convert 1abc to 1abc_
            $line .= '_';
        }
        push( @a, $line );
    }
    close(F);
    return (@a);
}

# ----------------------- get_path  ---------------------------------
# We have a filename and a list of directories where it could
# be. Return the path if we can find it, otherwise return undef.
sub get_path (\@ $) {
    my ( $dirs, $fname ) = @_;
    foreach my $d (@$dirs) {
        my $p = "$d/$fname";
        if ( -f $p ) {
            return $p;
        }
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
sub check_dirs (\@) {
    my $a    = shift;
    my $last = @$a;
    for ( my $i = 0 ; $i < $last ; $i++ ) {
        if ( !-d $$a[$i] ) {
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
sub check_files (\@ \@ $) {
    my ( $dirs, $fnames, $ext ) = @_;
    my $errors = 0;
    my @missing_files;
    foreach my $f (@$fnames) {
        my $name = "$f$ext";
        if ( !get_path( @$dirs, $name ) ) {
            $errors++;
            push (@missing_files, $name);
        }
    }
    return $errors;
}

# ----------------------- usage   -----------------------------------
sub usage () {
    print STDERR "Arguments wrong. Not starting\n";
    exit(EXIT_FAILURE);
}

# ----------------------- get_scores --------------------------------
sub get_scores ($ $ $ $) {
    my ( $pair_set, $coord1, $coord2, $to_use ) = @_;

    my ( $scr_tot, $coverage, $score1, $geo_gap, $score1_gap );
    my ( $str1, $crap );
    my ( $open_cost, $widen_cost, $nseq_gap );
    ( $score1_gap, $score1 ) = pair_set_score($pair_set);
    ( $str1,       $crap )   =
      pair_set_coverage( $pair_set, coord_size($coord1), coord_size($coord2) );
    $coverage = ( $str1 =~ tr/1// );    # This is coverage as an integer
    $coverage =
      $coverage / seq_size( coord_get_seq($coord1) )
      ;                                 #and as fraction of query structure

    my ( $k_scr2, $k_gap_geo, $k_seq_gap, $k_str_gap, $k_str_wdn );

    if ( $coverage < .05 ) {
        $geo_gap    = 0;
        $nseq_gap   = 0;
        $open_cost  = 0;
        $widen_cost = 0;
    }
    else {
        ( $open_cost, $widen_cost ) = pair_set_gap( $pair_set, 1, 1 );
    }

    $k_str_gap = 1 * $sw1_pgap_open;
    $k_str_wdn = 1 * $sw1_pgap_widen;

    $scr_tot = $score1 + $k_str_gap * $open_cost + $k_str_wdn * $widen_cost;
    return ( $scr_tot, $coverage, $score1, $score1_gap, $nseq_gap, $open_cost );
}

# ----------------------- close_up_and_mail -------------------------
# This could be a happy or unhappy exit.

sub close_up_and_mail ($ $ $)
{
    my ($subject, $address, $text_to_mail) = @_;

    $ENV{sendwait} = '1';    # Might persade mailer to finish before returning
    $ENV{encoding} = '8bit'; # Otherwise mailer thinks about quoted-printable

    if ( ! defined ($address)) {
        print STDERR "No address found ! $subject $address\n"; return undef;}

    my @cmdline = $mail_prog;

    if ( defined ($subject)) {
        push (@cmdline, '-s', $subject); }

    push (@cmdline, '-r', $mail_from_addr);
    push (@cmdline, $address);
    my $mail_ret = 0;
    if ($really_mail) {
        open (MAIL, '|-', @cmdline);
        print MAIL $text_to_mail;
        close (MAIL);
        $mail_ret = $?;
    } else {
        print "Wanted to say, $text_to_mail\n";
        print " I would invoke @cmdline\n";
        $mail_ret = 0;
    }

    if ($mail_ret != EXIT_SUCCESS) {
        return undef;}
    else {
        return 1;}

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
    thread_cleanup();
    restore_handlers();  # otherwise an ugly loop is possible
    print STDERR "Error: \"$msg\"\n";
    if (! defined ($title)) {
        $title = 'unknown'; }
    my $subject = "Failed calculating on $title";
    if (defined ($email_address)) {
        close_up_and_mail ($subject, $email_address, $msg) ;}
    exit (EXIT_FAILURE);
}


# ----------------------- get_alt_scores ---------------------------------
# calculates scores on random paths through the scoring matrix
# parameters: number_of_paths/scores, scoring_matrix, pair_set_of_optimal_path
# return: the scores
# This function takes a significant amount of cpu time. We should see
# if there is any scope optimisation.
sub get_alt_scores($ $ $) {
    my ( $num_scrs, $scr_mat, $pair_set ) = @_;
    my @scr_fin;

    for ( my $i = 0 ; $i < $num_scrs ; $i++ ) {
        $scr_fin[$i] = find_alt_path_score_simple( $scr_mat, $pair_set );}

    return \@scr_fin;
}

# ----------------------- normalize_alt_scores ------------------------------
#
sub normalize_alt_scores($) {
    my ($scrs_ref) = @_;
    my $mean = 0.0;

    foreach my $scr ( @{$scrs_ref} ) {
        $mean += $scr; }
    $mean /= @$scrs_ref;

    my $deviation = 0.0;

    foreach my $scr (@$scrs_ref) {
        my $tmp = $scr - $mean;
        $deviation += ( $tmp * $tmp );
    }
    $deviation /= ( @$scrs_ref - 1 );
    $deviation = sqrt($deviation);

    return ( $mean, $deviation );
}

# ----------------------- get_dme_thresh ----------------------------
# Given an alignment and two proteins, return the thresholded DME
# On error, return 0.0.
sub get_dme_thresh ($ $ $) {
    my ( $pair_set, $c1, $c2 ) = @_;
    my $model = make_model( $pair_set, coord_get_seq($c1), $c2 );
    if (! $model) {   # There is no error message here, since the C code
        return 0.0; } # already prints a warning.
    if ( coord_size($model) < 10 ) {
        return 0.0;
    }
    my $frac;
    if ( !dme_thresh( $frac, $c1, $model, 3.0 ) ) {
        print STDERR "dme_thresh broke.\n"; }
    return ($frac);
}

# ----------------------- do_align ----------------------------------
# This does the alignment.

sub do_align ($ $ $ $ $ $) {
    my ( $coord1,$coord2,$pvec1_strct,$pvec2_strct,$pvec1_ca,$pvec2_ca)= @_;
    my ($matrix_ca,$matrix_strct,$tot_matrix) = 0;
    my $seq_ptr1 = coord_get_seq($coord1);
    my $seq_ptr2 = coord_get_seq($coord2);

#   changed by Iryna
#   build score matrix
#   my $matrix = score_mat_new( seq_size($seq_ptr1), seq_size($seq_ptr2) );
# 
#   score_pvec( $matrix, $pvec1, $pvec2);
#   $matrix = score_mat_shift( $matrix, $m_shift );

    if ($function == 1 || $function == 0) {
        $matrix_ca =
            score_mat_new( seq_size($seq_ptr1), seq_size($seq_ptr2) );
        score_pvec( $matrix_ca,    $pvec1_ca,    $pvec2_ca );
    }
    if ($function == 2 || $function == 0) {
        $matrix_strct =
            score_mat_new(seq_size($seq_ptr1),seq_size($seq_ptr2) );
        score_pvec( $matrix_strct, $pvec1_strct, $pvec2_strct );
    }

#   shift and score matrix
    my $smallsize = (coord_size($coord1) > coord_size($coord2)
                      ? coord_size($coord2) : coord_size($coord1));
    my $shift = ($s_factor * $smallsize) + $s_const;
    if ($function == 0) {
        $tot_matrix =
            score_mat_add ($matrix_ca, $matrix_strct, $weight, $shift);
    } elsif ($function == 1) {
        $tot_matrix = score_mat_shift($matrix_ca, $shift);
    } else {
        # Andrew had another routine for the shift:
        my $scaler = $m_s_scale * $smallsize * 0.001;
        $tot_matrix = score_mat_shift($matrix_strct, $m_shift + $scaler);
    }
    #end change
    my (
        $sw_scr_tot,    $sw_coverage, $sw_score1,
        $sw_score1_gap, $sw_seq_gap,  $sw_strct_gap
    );
    my (
        $nw_scr_tot,    $nw_coverage, $nw_score1,
        $nw_score1_gap, $nw_seq_gap,  $nw_strct_gap
    );

    my $sw_pair_set = score_mat_sum_smpl(
        my $crap_mat,   $tot_matrix,   $sw1_pgap_open, $sw1_pgap_widen,
        $sw1_qgap_open, $sw1_qgap_widen, $S_AND_W
    );

    (
        $sw_scr_tot,    $sw_coverage, $sw_score1,
        $sw_score1_gap, $sw_seq_gap,  $sw_strct_gap
    )
        = get_scores( $sw_pair_set, $coord1, $coord2, 's_and_w' );
    my $frac_dme = get_dme_thresh( $sw_pair_set, $coord1, $coord2 );

    my $num_scrs     = 1000;
    my $alt_scrs_ref = get_alt_scores( $num_scrs, $tot_matrix, $sw_pair_set );

    my ( $mean, $deviation ) = normalize_alt_scores($alt_scrs_ref);

    undef($alt_scrs_ref);
    my $z_scr;
    if ( $deviation != 0 ) {
        $z_scr = ( $sw_scr_tot - $mean ) / $deviation;
    } else {
        $z_scr = 0.0; }    # Should not really happen

    #   If the alignment is tiny, one can get a ridiculous z-score
    if ( $sw_coverage < 0.03 ) {    # silently wipe these guys out
        $z_scr = 0.0; }

    my $sec_strct_pnlty = 0.0;
    my $newcoord;
    # ====== added for multiple alignment
    my $name;
    my $vectors;
    my $seq;
    my $rmsd;
    # ====================================
    
    my $patscor;
    my @r = (
        $sw_scr_tot, $rmsd, $nw_scr_tot,  $sw_coverage, $nw_coverage,
        $sw_score1,  $nw_score1,   $sw_pair_set, $z_scr, $frac_dme,
        $newcoord, $name, $vectors, $seq, $patscor
    );
    return ( \@r );
}

# ----------------------- zero_shift_mat ----------------------------
# This shifts a substitution matrix, not an alignment score matrix.
sub zero_shift_mat ($ $) {
    my ( $sub_mat, $shift ) = @_;
    for ( my $i = 0 ; $i < 20 ; $i++ ) {
        for ( my $j = $i ; $j < 20 ; $j++ ) {
            my $t = sub_mat_get_by_i( $sub_mat, $i, $j ) + $shift;
            sub_mat_set_by_i( $sub_mat, $i, $j, $t );
        }
    }
}

#added by Iryna
#------------------------ read_struct_pl -----------------------------
sub read_struct_pl ($ $) {
   my ($vecname, $coord1) = @_;
   my ($pvec_strct, $pvec_ca);
   if ($function == 1 || $function == 0) {
        my $vecfile_ca = "$PVEC_CA_DIR/$vecname.vec";


        if ( -e $vecfile_ca){
           $pvec_ca = prob_vec_read($vecfile_ca);
        } else {
            #my $CA_CLASSFILE =
            #    '/cluster/ploeffler/wurst_server/clssfcns/F7_ca';
            my $ca_classfcn = ac_read_calpha ($CA_CLASSFILE, $tau_error,
                                                   $ca_dist_error, $corr_num);
            $pvec_ca = calpha_strct_2_prob_vec($coord1, $ca_classfcn);
        }
    }
    if ($function == 2 || $function == 0) {
        my $vecfile_strct = "$PVEC_STRCT_DIR/$vecname.vec";
        if ( -e $vecfile_strct){
           $pvec_strct = prob_vec_read($vecfile_strct);
        } else {
            #my $strct_classfile =
            #    '/cluster/ploeffler/wurst_server/clssfcns/F6_struct';
            my $strct_classfcn =
                                aa_strct_clssfcn_read($CLASSFILE,$gauss_err);
           $pvec_strct = strct_2_prob_vec($coord1, $strct_classfcn);
        }
   }
   
   return ($pvec_strct, $pvec_ca);
}

# These two could be declared within main() and passed as arguments, but
# the argument list is already quite daunting.
# Declare them here and set them up in do_lib() just before the main loop.
# close them at the very end of the script.

use vars qw ( @low_reader_thr
              @low_reader_q
              @low_reader_done_q
              $prot_reader
              $prot_to_read_q
              $prot_ready_q
              $check_file_thr
              $prot_reader_first
              $thread_error
);

# ----------------------- ini_file_check ----------------------------
# Historically, we have had problems with files disappearing. This
# should run before we do any work, but it takes too long. We still want
# to check if files exist, so we have a compromise. This starts in the
# background, checks the files and reports any errors. If the main
# loop sees that we had a problem, it will stop immediately.
sub ini_file_check (\@ \@ $) {
    my ( $dirs, $fnames, $ext ) = @_;
    my $errors = undef;
    my @missing_files;
    my $count;
    foreach my $f (@$fnames) {
        my $name = "$f$ext";
        if ( !get_path( @$dirs, $name ) ) {
            $errors++;
            push (@missing_files, $name);
            if ($errors >= 20) {
                last;}
        }
    }
    if ($errors) {
        my $s;
        for (my $i = 0; $i < $errors; $i++) {
            $s = "$missing_files[$i]"; }
        $thread_error = "$thread_error Missing files: $s\n";
    }
    threads->exit();
}

# ----------------------- queue_measurer ----------------------------
# This is an interesting function for instrumenting, but not for
# routine use. We simply look at the length of the queue with
# proteins that are ready to be worked on. This is the only way we
# can survey if we are buffering enough, or more than we need.
{
    my $x = 1;
    sub queue_measurer ()
    {
        if ($x) {
            undef $x;
            setpriority (0, 0, 20);
            sleep 2;
        }
        my $i = 0;
        my @p;
        my @to_r;
        my @low_0;
        my @low_1;
        while (1) {
            $p[$i] = $prot_ready_q->pending;
            $to_r[$i] = $prot_to_read_q->pending;
            if ($i == 10) {
                $i = 0;
                print STDERR "ready   @p\nto_read @to_r\n";
            }
            sleep 1;
            $i++;
        }
    }
}

# ----------------------- start_threads  ----------------------------
# These are the file readers and corresponding queues.
sub start_threads($)
{
    my $a = shift;
    my $struct_dirs = $$a[0];
    my $struct_list = $$a[1];

    $check_file_thr = threads->create('ini_file_check', $struct_dirs, $struct_list, $BIN_SUFFIX);
    $prot_to_read_q = Thread::Queue->new();
    $prot_ready_q   = Thread::Queue->new();
    $prot_reader    = threads->create('prot_reader');
}


# ----------------------- get_local_name ----------------------------
# Given a name like /a/b/c/blah, returns a string like
# $local_storage/a_b_c_blah
sub get_local_name ($)
{
    my $in = shift;
    $in =~ s/\//_/g;
    $in = "$LOCAL_STORAGE/$in";
    return $in;
}

use File::Copy qw/ copy /;
# ----------------------- copy_to_local -----------------------------
# Given a source and destination, copy the file to local storage.
sub copy_to_local ($)
{
    my $q_num = shift;
    while (my $source = $low_reader_q[$q_num]->dequeue()) {
        if ($source eq 'STOP') {
            undef ($low_reader_q[$q_num]);
            threads->exit();}
        my $dest = get_local_name ($source);
        if (! File::Copy::copy ($source, $dest)) {
            my $err = "File copy failed from \n  $source to\n  $dest\n$!";
            $low_reader_done_q[$q_num]->enqueue('ERROR');
            $thread_error = "$err";
            warn "$err";
        } else {
            $low_reader_done_q[$q_num]->enqueue($dest);
        }
    }
}

$prot_reader_first = 1;
# ----------------------- prot_reader--------------------------------
# From the queue, we get
#    name
#    number of files to read
#    an array reference with the file names to copy to local storage
sub prot_reader ()
{
    my $a;
    if ($prot_reader_first) {             # At the start, we do not know
        undef ($prot_reader_first);       # how many threads we will need
        $a = $prot_to_read_q->dequeue();  # so just once, we start them.
        my $name = $$a[$#$a];
        for( my $i = 0; $i < $#$a; $i++) {
            $low_reader_q[$i] = Thread::Queue->new();
            $low_reader_done_q[$i] = Thread::Queue->new();
            $low_reader_thr[$i] = threads->create ('copy_to_local', $i);
        }
    }
    do {
        my $name = pop(@$a);
        if ($name eq 'STOP') {
            for (my $i = 0; $i <= $#low_reader_thr; $i++) {
                $low_reader_q[$i]->enqueue('STOP');}
            for (my $i = 0; $i <= $#low_reader_thr; $i++) {
                $low_reader_thr[$i]->join; }
            for (my $i = 0; $i <= $#low_reader_thr; $i++) {
                undef ($low_reader_q[$i]) }
            undef ($prot_to_read_q);
            threads->exit(); }
        if (! $name) {
            print STDERR "prot_reader() given empty name\n";
            return;
        }
        for (my $i = 0; $i < @$a; $i++) {
            $low_reader_q[$i]->enqueue($$a[$i]); }
        for (my $i = 0; $i < @$a; $i++) {
            my $r = $low_reader_done_q[$i]->dequeue();
            if ($r eq 'ERROR') {
                $name = 'ERROR' }
        }
        $prot_ready_q->enqueue ($name);
    } while ($a = $prot_to_read_q->dequeue())
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

sub do_lib (\@ \@ $ $ $ $ $ $ $ $ $ $ $ \@) {
    my ( $structlist, $struct_dirs, $query_struct, $coord1, $title,
        $rmsd_thresh, $minFracDME, $maxItNum, $jobfolder, $text_to_mail, $vis_dis, $filter, $filter_thresh, $focf) = @_;
    my (@pair_sets);
    my ($vecname, $coord2, $pvec1_strct, $pvec1_ca, $pvec2_strct, $pvec2_ca);
    use File::Basename;
    $vecname = basename($query_struct, ".bin");
    my $queryseq = coord_get_seq($coord1);

    ($pvec1_strct, $pvec1_ca) = read_struct_pl($vecname, $coord1);

    if (@$structlist == 0){
        log_time("NO STRUCTURES TO ALIGN");
        return EXIT_SUCCESS;
    }
    if (@$structlist != @$focf){
        log_time("LENGTH STRUCTLIST != FOCF");
        return EXIT_FAILURE; 
    }
    # Fill buffer
    my $minsize = seq_size( coord_get_seq($coord1) );
    my $n_to_buf = 10;
    if ($n_to_buf >= ($#$structlist + 1)) {
        $n_to_buf = $#$structlist;} # Very unlikely, but does happen when debugging
    for ( my $i = 0; $i < $n_to_buf; $i++) {
        my @a;
        my $name = $$structlist[$i];
        if ($function == 1 || $function == 0) {
            push (@a, "$PVEC_CA_DIR/$name.vec");
        }
        if ($function == 2 || $function == 0) {
            push (@a, "$PVEC_STRCT_DIR/$name.vec");
        }
        push (@a, "$STRUCT_DIR/$name$BIN_SUFFIX");
        push (@a, $name);
        log_time("enqueued $name");
        $prot_to_read_q->enqueue (\@a);
    }

    log_time("Before main loop");
    # LOOP
    # Go over all structures
    for ( my $i = 0; $i < @$structlist; $i++) {
        if ( $thread_error) {
            print STDERR "$thread_error";
            return EXIT_FAILURE;
        }
        my $libprotname = $prot_ready_q->dequeue();
        if ($libprotname eq 'ERROR') {
            print STDERR "Error from reader thread on $$structlist[$i]\n";
            return EXIT_FAILURE;
        }
        if ($libprotname ne $$structlist[$i]) {
            print STDERR "Reading from queue out of sync with file list\n",
            "Probably a programming bug\nStopping\n";
            return EXIT_FAILURE;
        }
        my @for_reader;
        my $for_reader_n = 0;
        if ($function == 1 || $function == 0) {
            my $tname = "$PVEC_CA_DIR/$$structlist[$i].vec";
            my $cached = get_local_name ($tname);
            $pvec2_ca = prob_vec_read($cached);
            if (!$pvec2_ca) {
                print STDERR "$tname (cached) not found\n";
                if (! -e $tname) {
                    print STDERR "$tname does not exist\n";}
                return EXIT_FAILURE;
            }
            unlink ($cached);
            if ( $i < (@$structlist - $n_to_buf)) {
                my $next = "$PVEC_CA_DIR/$$structlist[$i+$n_to_buf].vec";
                push (@for_reader, $next);
            }
        }
        if ($function == 2 || $function == 0) {
            my $tname = "$PVEC_STRCT_DIR/$$structlist[$i].vec";
            $pvec2_strct = prob_vec_read ($tname);
            my $cached = get_local_name($tname);
            if (!$pvec2_strct) {
                print STDERR "$tname (cached) not found\n";
                if ( ! -e $tname) {
                    print STDERR "$tname does not exist\n";}
                return EXIT_FAILURE;
            }
            unlink ($cached);
            if ( $i < (@$structlist - $n_to_buf)) {
                my $next = "$PVEC_STRCT_DIR/$$structlist[$i+$n_to_buf].vec";
                push (@for_reader, $next);
            }
        }
        {
            my $cached = get_local_name("$STRUCT_DIR/$$structlist[$i]$BIN_SUFFIX");
            $coord2 = coord_read($cached);
            if (!$coord2) {
                warn "coord_read $$structlist[$i] failed";
                return EXIT_FAILURE;
            }
            unlink ($cached);
        }
        if ( $i < (@$structlist - $n_to_buf)) {
            my $tname = $$structlist[$i+$n_to_buf];
            push (@for_reader, "$STRUCT_DIR/$tname$BIN_SUFFIX");
            push (@for_reader, $tname);
            $prot_to_read_q->enqueue (\@for_reader);
        }

        #added by Iryna
        if ($function == 1) { # undefine the unused pvec
            $pvec1_strct = undef;
            $pvec2_strct = undef;
        } elsif ($function == 2) {
            $pvec1_ca = undef;
            $pvec2_ca = undef;
        }

        $r::r[$i] = do_align( $coord2, $coord1, $pvec2_strct, $pvec1_strct,
                                                $pvec2_ca,    $pvec1_ca );

        $r::r[$i][4] =
          get_seq_id_simple( $r::r[$i][7], coord_get_seq($coord2),
            coord_get_seq($coord1) );
        $r::r[$i][6] = $r::r[$i][3] * seq_size( coord_get_seq($coord2) );
        if ( $r::r[$i][6] < 25 ) {
            $r::r[$i][8] = 0;
        }
        $r::r[$i][2] = $r::r[$i][9] * ( ( $r::r[$i][6] ) / $minsize );
        $r::r[$i][11] = $r::r[$i][6] ** $r::r[$i][9];

    }
    log_time("Main loop done");

    my @indices;
    for (my $i = 0 ; $i < @$structlist ; $i++ ) {
        $indices[$i] = $i;
    }
#   changed by Iryna
#   @indices = sort { $r::r[$b][2] <=> $r::r[$a][2]; } @indices;
    @indices = sort { $r::r[$b][11] <=> $r::r[$a][11]; } @indices;

    #DIRS
    my $htmldir = $RESDIR.$jobfolder;
    my ($htmlbase, undef, undef) = fileparse ($htmldir);
    my $url = "$RESURL/jobs/$htmlbase/index.html";
    $$text_to_mail .= "SALAMI has finished searching for proteins similar to $vecname.\n";
    $$text_to_mail .= " Your results can be viewed in our interactive results browser at:\n";
    $$text_to_mail .= " $url\n";
    coord_2_pdb("$htmldir/modeldir/$vecname.pdb", $coord1);

    my $todo = ( @$structlist > $maxItNum ? $maxItNum : @$structlist );

    if ( ! -d $modeldir ) {
        if ( !mkdir( "$modeldir", 0777 ) ) {
            bad_exit("Fail create modeldir ($modeldir): $!"); } }
    my $p1 = coord_name($coord1);

    print STDERR " TODO: $todo merge_alignments............\n";
    my $veclist = initveclist($todo+1);
    my $seqlist = initseqlist($todo+1);
    print STDERR "initlist done \n";

    my $pvec1;
    if ($function == 2 || $function == 0) {
        #XXX todo: use both vectors
        $pvec1 = $pvec1_strct;
        addvec($veclist, prob_vec_read("$PVEC_STRCT_DIR/$$structlist[$indices[0]].vec"));
    }
    if ($function == 1){
        $pvec1 = $pvec1_ca;
        addvec($veclist, prob_vec_read("$PVEC_CA_DIR/$$structlist[$indices[0]].vec"));
    }
#    addvec($veclist, prob_vec_read("$PVEC_STRCT_DIR/$$structlist[$indices[0]].vec"));
    my $seq = coord_get_seq(
                          coord_read(
                                  get_path( @$struct_dirs, $$structlist[$indices[0]] . $BIN_SUFFIX)));
    addseq($seqlist, $seq);
    addvec($veclist, $pvec1);             #add the query sequence's vector
    addseq($seqlist, coord_get_seq($coord1));
    print STDERR "addvec done todo = $todo \n";
    my $malignm = $r::r[$indices[0]][7];
    for (my $i = 1 ; $i < $todo ; $i++ ) {
      if ($function == 2 || $function == 0) {
          addvec($veclist, prob_vec_read("$PVEC_STRCT_DIR/$$structlist[$indices[$i]].vec"));
      }
      else {
          addvec($veclist, prob_vec_read("$PVEC_CA_DIR/$$structlist[$indices[$i]].vec"));
      }
       $seq = coord_get_seq(
                            coord_read(
                                       get_path( @$struct_dirs, $$structlist[$indices[$i]] . $BIN_SUFFIX )));
       addseq($seqlist, $seq);


       $malignm = merge_localigns($malignm, $r::r[$indices[$i]][7], 1, 1);
       $malignm = remove_seq($malignm, get_pair_set_m($malignm)-1);
    }
    my $test = get_pair_set_m($malignm);  #test

    print STDERR "m = $test \n";
    print STDERR "str X = ", pair_set_string($malignm, coord_get_seq($coord1), coord_get_seq($coord1));
    print STDERR "pair_set_stringI \n";
    print STDERR  pair_set_stringI($malignm, 1, $queryseq), "\n";
    #conservierung computing
    print STDERR "getconservvec \n";
    my $conservs = getconservvec($malignm, $veclist);
    print STDERR "get_seq_conserv \n";
    my $seqconserv = get_seq_conserv($malignm, $seqlist);
    # ====================================
    my $j = 0;
    my $shift = 0;
    my @proteins = ();

MINFRAGDME: { #for the dme thresh
    for (my $i = 0 ; $i < ($todo + $shift); $i++ ) {
        my $idx      = $indices[$i];
        my $pair_set = $r::r[$idx][7];
        my $a        = $r::r[$idx];

        if ($$a[9] < $minFracDME) {
           if (@$structlist < ($todo + $shift)) {
              last;
           } else {
              $shift++;
              next; #MINFRAGDME;
           }
        }

        my $coord2   =
        coord_read(
              get_path( @$struct_dirs, $$structlist[$idx] . $BIN_SUFFIX ) );

        ($$a[1], $$a[10], my $crap_b ) =
          coord_rmsd( $$a[7], $coord2, $coord1, 0 );
        my $p2 = $$structlist[$idx];
        my $sid;
        if ( $$a[6] != 0 ) {
            $sid = $$a[4] / $$a[6];
            $sid *= 100;
        }
        else {
            $sid = 0.0;
        }
        my $pdbid=$$structlist[$idx];
        push(@proteins, $pdbid);

        my $alistr;
        my $seq1 = coord_get_seq($coord1);
        my $seq2 = coord_get_seq($coord2);
        my $testbitset = pair_set_sel_geti($pair_set);
        my $set_alg_tmp = selected_pair_set_get($pair_set, $testbitset);
        ($$a[1], $$a[10], my $coord11) = coord_rmsd( $set_alg_tmp, $$a[10], $coord1, 0);
        my $n = get_pair_set_n($set_alg_tmp);
        while ( ($$a[1] > $rmsd_thresh) && ($n > $MAG_NUM_ASREST)) {
            ($$a[1], $$a[10], my $coord11) = coord_rmsd( $set_alg_tmp, $$a[10], $coord11, 0);
            $testbitset = pair_set_sel_delmaxdistance($$a[10], $coord11, $pair_set, $testbitset);
            $set_alg_tmp = selected_pair_set_get($pair_set, $testbitset);
            $n = get_pair_set_n($set_alg_tmp);
        }
        my $pdbi = $pdbid;
        $pdbi =~ s/(.*?)(.)$/$1/;

        printf RESULTS "
  <result id='$pdbid'>
    <pdbid>$pdbid</pdbid>
    <z_scr>%8.3g</z_scr>
    <f_dme>%8.3g</f_dme>
    <sw_cvr>%6.2f</sw_cvr>
    <seq_id>%.0f</seq_id>
    <asize>%8.4g</asize>
    <rmsd>%6.2f</rmsd>
    <q_scr>%4.4f</q_scr>
    <focf>%1.3f</focf>
  </result>
" ,$$a[8], $$a[9], $$a[6]/seq_size(coord_get_seq($coord1)), $sid, $$a[6], $$a[1], $$a[2], $$focf[$i];

        my $alistr1 = pair_set_stringI($pair_set, 1, $seq1);
        my $alistr2 = pair_set_stringI($pair_set, 0, $seq2);
        my $num1 = pair_set_get_strNum($pair_set, 1, $coord1);
        my $num2 = pair_set_get_strNum($pair_set, 0, $coord2);
        my $used = pair_set_sel_print($testbitset, get_pair_set_n($pair_set));
        my $sp1 = pair_set_get_startpos($pair_set, 1);
        my $sp2 = pair_set_get_startpos($pair_set, 0);
        my $conserv;
        my $strc;
        if ($i < 1) {
          $conserv = printconserv($conservs, get_pair_set_n($malignm), $malignm, $sp1, get_pair_set_n($pair_set), $j);
          $strc = printconserv($seqconserv, get_pair_set_n($malignm), $malignm, $sp1, get_pair_set_n($pair_set), $j);
        } else {
          $conserv = printconserv($conservs, get_pair_set_n($malignm), $malignm, $sp1, get_pair_set_n($pair_set), $j+1);
          $strc = printconserv($seqconserv, get_pair_set_n($malignm), $malignm, $sp1, get_pair_set_n($pair_set), $j+1);
        }
        $j++;

        open ALI, ">$htmldir/json/$$structlist[$idx].json";
        print ALI "{
    \"chains\" : {
        \"$$structlist[$idx]\" : \"$alistr2\",
        \"$vecname\" : \"$alistr1\"
    },
    \"$$structlist[$idx]_positions\" :
    [$num2],
    \"".$vecname."_positions\" :
    [$num1],
    \"$$structlist[$idx]_startpos\" :
    $sp1,
    \"".$vecname."\_startpos\" :
    $sp2,
    \"conservation\" :
    [$conserv],
    \"seq_conservation\" :
    [$strc],
    \"used\" :
    \"$used\"
}";
        coord_2_pdb("$htmldir/modeldir/$$structlist[$idx].pdb", $$a[10]);
    close ALI;    
    }
    } #minFracDME:

    # Visulize distances code, produces pca.xml
    if ($vis_dis eq "Y"){
       log_time("Before similarity calculation");
       open PCA, ">$htmldir/pca.xml";
       print STDERR "vis_dis YES\n",
       my $builder = new Similarity::Builder({
          ref1 => \@proteins,
          ref2 => \@proteins,
          pcc => 3
       });
       print PCA $builder->xml($builder->all_to_all());
       close PCA;
       copy_to_flensburg("pca.xml", "", $jobfolder);
       log_time("After similarity calculation");
    } else {
        print STDERR "vis_dis NO\n",
    }

    if($N_MODELS > 0){
        coord_2_pdb("$vecname.pdb", $coord1);
        for (my $i = 0 ; $i < $N_MODELS && $i < $todo ; $i++ ) {
            my $idx = $indices[$i];
            my $a = $r::r[$idx];
            coord_2_pdb("$modeldir/$$structlist[$idx].pdb", $$a[10]);
        }
    }
    undef(@r::r);
    print STDERR "***************\nhtmldir=$htmldir \nOUTURL=$OUTURL \nhtmlbase=$htmlbase \nRESDIR=$RESDIR \nurl=$url \njobfolder=$jobfolder";

    copy_to_flensburg("json", "", $jobfolder);
    copy_to_flensburg("modeldir", "", $jobfolder);
    log_time("Finished copying to flensburg");

    return EXIT_SUCCESS;
}

# ----------------------- copy_to_flensburg ----------------------------
# copies files to flensburg. they all land into the declared
# jobfolder

sub copy_to_flensburg
{
   my ($what, $new_name, $jobfolder) = @_;

   my $ret = `rsync -a $RESDIR$jobfolder/$what wurst\@flensburg:$OUTURL/$jobfolder/$new_name 2>&1`;
   if ($? != EXIT_SUCCESS) {
       print STDERR "Fail copying results back to flensburg: $!\n",
       "scp command said\n$ret\n";
       print STDERR "Please die and put a message on the web page\n";
   }
   print STDERR $ret;
}


# ----------------------- log_time ----------------------------
# writes time into STDERR
#
sub log_time
{
   my ($message) = @_;
   print STDERR " ### ". localtime." ### ". $message. "\n";
}



# ----------------------- catch_kill     ----------------------------
# The main thing is, if we get a KILL or TERM, to call exit and get
# out of here. This means there is a better chance of closing files
# wherever we were up to.
sub catch_kill
{
    my ($sig) = @_;
    thread_cleanup();
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

# ----------------------- thread_cleanup ----------------------------
# This used to be a general cleanup, but now the threads are closed
# as soon as possible, so this is only called when there is a problem.
sub thread_cleanup ()
{
    my @t;
    push (@t, 'STOP');
    $prot_to_read_q->enqueue(\@t);
    $prot_reader->join;
    $check_file_thr->join;
    undef ($prot_to_read_q);
    undef ($prot_to_read_q);
    undef ($prot_ready_q);
}

# ----------------------- mymain  -----------------------------------
# Arg 1 is a structure file. Arg 2 is a structure list file.
sub mymain () {
    $DB::single = 1;
    use threads;
    use Thread::Queue;
    use threads::shared;
    use File::Temp qw/ tempdir /;

    if ( !($LOCAL_STORAGE = tempdir ( DIR => $LOCAL_STORAGE, CLEANUP => 1))) {
        print STDERR "Failed to make tempdir in $LOCAL_STORAGE\n$!\n";
        return EXIT_FAILURE;
    }
    print "start of mymain\n";
    use Getopt::Std;
    my (%opts);
    my ( @struct_list, @struct_dirs );
    my ( $structfile,  $libfile );
    my $fatalflag = undef;
    @struct_dirs = @DFLT_STRUCT_DIRS;

#--------------------------------server stuff-------------------------

#   Set up directory for output and models. It is the only place
#   we will be allowed to write to, so cleanup is easy after disaster

    if ( ! -d $TOP_TEMP ) {
        mkdir ($TOP_TEMP) || bad_exit ("Fail creating $TOP_TEMP: $!");     }
    if ( ! chdir ($TOP_TEMP)) {
        bad_exit ("Failed to cd to $TOP_TEMP: $!"); }

    # *+-* TEMPDIR AT $TOP_TEMP ON CLUSTER AT  home/other/wurst/wurst_delete....

    $tmp_dir  = tempdir (DIR => $TOP_TEMP,  CLEANUP => 1);

    chmod 0777, $tmp_dir;      # So normal people can look for disasters.
    $modeldir = "$tmp_dir/modeldir"; # Where models go. Automatically cleaned.
#   If the machines get rebooted or jobs are killed, try to mail
#   back some information. Trap the errors.
    kill_handlers();

    my $email_address = "Anonymous";
    my $minFracDME;
    my $maxItNum;
    my $rmsd_t;
    my $jobfolder;
    my $vis_dis;
    my $filter;

    if ($fake_args) {    # Set $fake_args (above) to a non zero value
        my $junkdir = "$ENV{HOME}/junk";  # will turn on these values
        $opts{o} = 'fake_arg_folder';     # for testing and
        $opts{n} = 'fake_search';         # debugging.
        $opts{e} = 'rubbish@mailinator.com';
        $opts{q} = "$junkdir/4fyuA.pdb";
        $opts{l} = '/smallfiles/public/no_backup/bm/pdb90.list';
        $opts{l} = "$junkdir/short.list";
        $opts{l} = "$junkdir/8.list";
        $opts{s} = '0';
        $opts{r} = '3';
        $opts{f} = '0.75';
        $opts{i} = '100';
        $opts{p} = 'str';
        $opts{c} = 'N';
        $opts{x} = 'Y';
    }

    if ( !getopts( 'a:d:e:n:h:s:t:q:l:r:f:i:o:x:p:b:c:', \%opts ) ) {
        usage(); }
    if ( defined( $opts{a} ) ) { $N_ALIGNMENTS = $opts{a} }
    if ( defined( $opts{d} ) ) { $modeldir     = $opts{d} }
    if ( defined( $opts{o} ) ) { $jobfolder    = $opts{o} }
    if ( defined( $opts{h} ) ) { $N_BRIEF_LIST = $opts{h} }
    if ( defined( $opts{s} ) ) { $N_MODELS     = $opts{s} }
    if ( defined( $opts{e} ) ) { $email_address= $opts{e} }
    if ( defined( $opts{n} ) ) { $title	       = $opts{n} }
    if ( defined( $opts{q} ) ) { $structfile   = $opts{q} }
    if ( defined( $opts{c} ) ) { $vis_dis      = $opts{c} }
    if ( defined( $opts{x} ) ) { $filter       = $opts{x} }
    else{
        print STDERR "Must have at least a query structure file\n";
        usage();
    }
    if ( defined( $opts{l} ) ) { $libfile      = $opts{l} }
    else{
        print STDERR "Please give me a structure library / file\n";
        usage();
    }
    if ( defined( $opts{t} ) ) {
        push( @struct_dirs, split( ',', $opts{t} ) );
    }
    if ( defined( $opts{r} ) ) { $rmsd_t       = $opts{r} }
    if ( defined( $opts{f} ) ) { $minFracDME   = $opts{f} }
    if ( defined( $opts{i} ) ) { $maxItNum     = $opts{i} }
    # new and old vectors (smallfiles))
    if ( defined( $opts{p} ) ) {
      if ($opts{p} eq 'ca' ) {
        $function = 1;
        # change parameters to best ca parameters
        # 2902 set: -0.698265545936722 with..
        $sw1_pgap_open = $sw1_qgap_open = 2.421;
        $sw1_pgap_widen = $sw1_qgap_widen = 1.133;
        $s_factor = -0.001941;
        $s_const = 0.02312;
      } elsif ($opts{p} eq 'str') {
        $function  = 2;
        # I used Andrews o.336, best for monsterset
        $sw1_pgap_open = $sw1_qgap_open = 1.07;
        $sw1_pgap_widen = $sw1_qgap_widen = 0.855;
        $m_shift = -0.885;
        $m_s_scale = -0.4501;
      } else {
        print STDERR "I dont know this function, try ca/str\n\n";
        usage();
      }
    }
    undef %opts;
    if (! $fake_args) {
        log_job ($email_address, 'start', $title); }

    set_params();
    my $query_struct;
    my $query_struct_name = $structfile;

    check_dirs(@struct_dirs);
    if ( @struct_dirs == 0 ) {
        die "\nNo valid structure directory. Stopping.\n"; }
    $DB::single = 1;
    ( @struct_list = get_prot_list($libfile) ) || $fatalflag++;
    share ($thread_error);
    {
        my @checker_list=  (\@struct_dirs, \@struct_list, $BIN_SUFFIX);
        start_threads(\@checker_list);
    }

    if ($filter eq "Y"){
        print STDERR "filter YES\n";
    } else {
        print STDERR "filter NO\n";
    }

    if ($fatalflag) {
        print STDERR "struct dirs were @struct_dirs\n";
        print STDERR "Fatal problems\n";
        return EXIT_FAILURE;
    }

    $query_struct = pdb_read( $query_struct_name, '', '' );
    if ($query_struct == 0) {
        print STDERR "\npdb_read of query failed. ($query_struct_name)\n";
        return EXIT_FAILURE;
    }

    log_time("before creating dirs");
    my $htmldir = $RESDIR.$jobfolder;
    my $mkdret = mkdir $htmldir;
    if ($mkdret != 1) {
        if (-e $htmldir) {
            print STDERR "Output folder \"$htmldir\" exists. Possible problem.\n";
        } else {
            print STDERR "makedir failure $jobfolder in  $RESDIR. Must die. $!";
            return EXIT_FAILURE;
        }
    }
    $htmldir =~ s/\n$//;
    mkdir("$htmldir/modeldir");
    mkdir("$htmldir/json");
    log_time("after creating dirs");

    ########### FILTER #########################
    my $prob_vec = strct_2_prob_vec($query_struct, paa_strct_clss6());

    my @sorted_structs = sort @struct_list;
    my $cat_list = join(",", @sorted_structs);

    log_time("Before get_focf");
    my $target_file = "/home/other/wurst/salamiServer/v02/fpc.bin";
    #write_fpc_bin($PVEC_STRCT_DIR, $new_target); #holy writing!
    my $focf = get_focf($cat_list, @sorted_structs, $prob_vec, $target_file);
    log_time("After get_focf, before get_threshs");
    my $num_threshs = 4;
    my $c_threshs = get_threshs($focf, @sorted_structs, 5, $num_threshs);
    my @thresholds = (1.0);
    my $tt;
    if ($vis_dis eq "N"){
       for (my $i=0; $i<$num_threshs; $i=$i+1){
          $tt = get_float_value($c_threshs, $i)*2;
          $tt = sprintf("%1.1f", $tt);
          $tt = $tt/2;
          if ($thresholds[@thresholds-1] != $tt){
             push(@thresholds, $tt);
          }
       }
    }
    push(@thresholds, 0.0);
    log_time("After get_threshs and arrangment of perl cutoff array");

    log_time("Before sorting structures");
    my @struct_lists;
    my @focf_lists;
    for (my $i=0; $i<@thresholds; $i++){
       push(@{$struct_lists[$i]}, ());  
       push(@{$focf_lists[$i]}, ());  
    }
    push(@struct_lists, ());  
    push(@focf_lists, ());  

    open FOCF, ">$RESDIR$jobfolder/focf.table";
    print FOCF "ID, focf, bin\n";
    my $rounded;
    my @focf_filtered;
    my @rest_structs; my @rest_focf; # used for next iteration

    # binning
    my $v;
    for ( my $i = 0; $i < @sorted_structs; $i++) {
       $v = get_float_value($focf, $i);
       $rounded = sprintf("%1.3f", $v);
       my $f;
       for ($f=1; $f<@thresholds; $f++){
          if ($rounded >= $thresholds[$f] && $rounded < $thresholds[$f-1]) {
             push(@{$struct_lists[$f-1]}, $sorted_structs[$i]);
             push(@{$focf_lists[$f-1]}, $rounded);
             print FOCF "$sorted_structs[$i], $rounded, ". ($f-1) ."\n";
             last;
          }
        }
    }
    close FOCF;
    log_time("Done with writing to FOCF and structure sorting");
    copy_to_flensburg("focf.table", "", $jobfolder);
    $query_struct_name =~ s/\..+$//;
    my $text_to_mail = "";

    # used for server to recognize results
    open READY, ">$RESDIR$jobfolder/results.ready";
    print READY "the results are ready now\n";
    close READY;

    # writing parameters to results.xml
    open RESULTS, ">$RESDIR$jobfolder/results.xml";
    print RESULTS "
<query>
 <params>
   <structure>".basename($query_struct_name)."</structure>
   <title>$title</title>
   <rmsd_thresh>$rmsd_t</rmsd_thresh>
   <min_f_dme>$minFracDME</min_f_dme>
   <max_n>$maxItNum</max_n>
   <jobfolder>$jobfolder</jobfolder>
   <vis_dis>$vis_dis</vis_dis>
   <filter>$filter</filter>
 </params>
 <results>
";

    log_time("*** Starting do_lib loop: ***");
    for (my $f=1; $f<@thresholds; $f++){
       my $r = do_lib( @{$struct_lists[$f-1]}, @struct_dirs,
                    $query_struct_name, $query_struct,
                    $title, $rmsd_t, $minFracDME, $maxItNum,
                    $jobfolder, \$text_to_mail, 
                    $vis_dis, $filter, $thresholds[$f], @{$focf_lists[$f-1]});
       if ( $r == EXIT_FAILURE ) { bad_exit('calculation broke'); }
       log_time(" ------------ After do_lib ". ($f-1). " ----------------------------");
       
       if ($vis_dis eq "Y"){
          last;
       } else {
          # copy pre-results of specific cutoff
          my $rescut = "results_". ($thresholds[$f]*100) . ".xml";
          close RESULTS;
          copy("$RESDIR$jobfolder/results.xml", "$RESDIR$jobfolder/$rescut") or die "Copy failed on $rescut";
          open RESULTS, ">>$RESDIR$jobfolder/results.xml";
          open RESCUT, ">>$RESDIR$jobfolder/$rescut";
          print RESCUT " </results>\n</query>";
          close RESCUT;
          copy_to_flensburg($rescut, $rescut, $jobfolder);
          copy_to_flensburg("results.ready", "", $jobfolder);
       }
    }

    # copy results of specific cutoff
    print RESULTS " </results>\n</query>";
    close RESULTS;
    copy_to_flensburg("results.xml", "", $jobfolder);
    copy_to_flensburg("results.ready", "", $jobfolder);

    {
        my @t;
        push (@t, 'STOP');
        $prot_to_read_q->enqueue(\@t);
        $prot_reader->join();
        $check_file_thr->join();
    }

    
    use Sys::Hostname;
    my $fmt = "Finished on %s (Hamburg time) on %s\n";
    my $s = sprintf ($fmt, scalar(localtime()), hostname());
    $text_to_mail = "$text_to_mail$s\n";
    
    if($email_address ne "Anonymous"){
        close_up_and_mail ($title, $email_address, $text_to_mail); }
#   We could call thread_cleanup() now, but at the moment, there
#   are no thread leaks and it is not necessary.
#   thread_cleanup();
    if ( ! $fake_args) {
        log_job ($email_address, 'stop', $title);}

    log_time("Finished with all");
    return EXIT_SUCCESS;
}

# ----------------------- main    -----------------------------------
exit( mymain() );
