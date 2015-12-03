# 22 May 2002
# This is the monster script to do a big calculation, taking a
# pile of sequences and aligning each to everyone in a library.
# For each alignment a model is calculated and its score and gap
# penalties calculated. The closeness of the model is also
# compared to the native structure.
# rscid = $Id: ff2_param.pl,v 1.1 2007/09/28 16:57:02 mmundry Exp $

=pod
=head1 NAME

ff2_param.pl - align a set of sequences to a set of structures
    and write lots of scores to a file for later optimisation

=head1 SYNOPSIS

ff2_param.pl [ options ] seq_file_name struct_file_name output_sw output_nw

=head1 DESCRIPTION

A nice program.

=head2 COMMANDLINE

The commandline must give the following:

=over

=item seq_filename

This is a file of protein whose sequences will be aligned across
    the library. Do not give names corresponding to sequence
    files. Give names corresponding to structure files with an
    extension like F<.bin>. The exact extension is hardwired near
    the top of the script.  The format is lists of protein codes,
    with or without filetypes appended. If a name has a filetype,
    it will be thrown away and changed to the default
    extension. Anything after a hash, C<#>, will be treated as a
    comment. Blank lines are allowed. The file looks like

     # Comment with date
     1abcA      # protein 1abc, chain A
     1qrs_      # protein 1qrs, no chain specified.
     1xyz       # If no chain is specified, it will be converted to 1xyz_
     1def_.bin  # File type (.bin) will be tossed out
     2abcA.phd  # This file type will also be thrown away.

For every file, the code will check to see if a F<.bin> file and
    a F<.phd> file exists.

=item struct_filename

This is identical to F<seq_filename> above, but is the list of
    structures to be aligned to. It defines our library.

=item output_sw

This is where lots of output will go to. It is the output from
all the Smith and Waterman models and comparisons.

=item output_nw

This is the same as the output from output_sw, but will contain
numbers calculated from Needleman and Wunsch alignments.

=back

=head2 OPTIONS

=over

=item -p PHD_DIRECTORY[,PHD_DIR]

Add (not replace) PHD_DIRECTORY to the list of places where we
will look for phd files. This is a comma separated list of alternative
directories. No spaces are allowed in the list.

=item -s SEQ_DIRECTORY[,SEQ_DIR]

Add (not replace) the directory SEQ_DIRECTORY to the list of
places where we will look for *.bin files in order to get
sequences.

=item -t STRUCT_DIRECTORY[,STRUCT_DIR]

Add (not replace) the directory STRUCT_DIRECTORY to the list of
places where we will look for F<*.bin> files in order to find
structures.

=item -v VERBOSITY

Set verbosity to VERBOSITY. A big number means noisy. A zero
means be very quiet.

=back

=cut
use FindBin;
use lib "$FindBin::Bin../blib/arch";
use lib "$FindBin::Bin../blib/lib";

use strict;
use POSIX qw(EXIT_SUCCESS EXIT_FAILURE);

use Wurst;

use vars qw ( $MATRIX_DIR $MATRIX_FILE $PARAM_DIR
              $FX9_PARAM_FILE $RS_PARAM_FILE );
do "$FindBin::Bin/paths.inc" ||
    die "Error reading paths";
if ($@) {
    die "broke reading paths.inc:\n$@"; }


# ----------------------- Various constants    ----------------------
use vars qw($phd_suffix $bin_suffix
            @DFLT_SEQ_DIRS @DFLT_STRUCT_DIRS @DFLT_PHD_DIRS);
*phd_suffix       = \'.phd';
*bin_suffix       = \'.bin';
*DFLT_SEQ_DIRS    = ['/rsc/appenzeller/data1/pdblib6apr00',
                     '/rsc/appenzeller/data1/pdb_bin_pool'];
*DFLT_STRUCT_DIRS = ['/rsc/appenzeller/data1/pdblib6apr00',
                     '/rsc/appenzeller/data1/pdb_bin_pool'];
*DFLT_PHD_DIRS    = ["$ENV{HOME}/phd/results"];

use vars qw ($DFLT_VERBOSITY);
*DFLT_VERBOSITY = \1;

# Threshold for calculating thresholded-DME-based similarity
use vars qw ($dmethresh $viol_max $scale);
*dmethresh = \4.0;
*viol_max  = \10.0;  # Geometric gaps bigger than this are set to this
*scale     = \1.0;   # The scale with which gap penalties will be multiplied


# ----------------------- alignment parameters ----------------------
# These guys come from pasting from an optimisation run.
# Don't make fun of all the decimal places.
# The first set of parameters are for the initial score matrix
# and Smith and Waterman. The second (smaller) set of penalties
# will be used for extension of alignments via a Needleman and
# Wunsch.
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
            print "$$a[$i] not a valid directory. Removing from list\n";
            splice @$a, $i, 1;
            $last--;
            $i--;
        }
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
# We have a filename, an extension and a list of directories
# where it could be. Return the path if we can find it, otherwise
# return undef.
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

# ----------------------- write_info --------------------------------
# We have one sequence, coordinates and an alignment. Write
# anything we can think of to a file.
sub write_info ($ $ $ $ $ $ $ $)
{
    my ($pair_set, $seq, $seqname, $coord, $structname, $rs_prams,
        $ofile, $c_seq) = @_;
    my ($big_score, $smpl_score);
    my ($geo_quad, $geo_linear, $geo_logi, $n_seq_gap);
    my ($open_cost, $widen_cost);
    my $model;
    my ($score2);
    $model  = make_model ($pair_set, $seq, $coord);
    ($big_score, $smpl_score) = pair_set_score ($pair_set);
    ($geo_quad, $geo_linear, $geo_logi, $n_seq_gap) =
                                     coord_geo_gap ($model, $scale, $viol_max);
    ($open_cost, $widen_cost) = pair_set_gap ($pair_set, $scale, $scale);
    
    my ($tdme, $coverage);
    dme_thresh ($tdme, $model, $c_seq, $dmethresh) ||
        print STDERR "Horrible DME error on seq ", coord_name($c_seq),
        ' aligned to ', coord_name ($coord);

    my ($str1, $str2) = 
        pair_set_coverage ($pair_set, seq_size ($seq), coord_size ($coord));
    $coverage = ($str1 =~ tr/1//);  # This is coverage as an integer
    $coverage = $coverage / seq_size ($seq); # and as fraction of sequence
    $score2   = score_rs ($model, $rs_prams);  # This is the second force field

    print $ofile 
        "$tdme ",
        coord_size($coord),
        " $coverage $smpl_score $score2 $geo_quad  $geo_linear  $geo_logi ",
        "$n_seq_gap  $open_cost  $widen_cost ", coord_name($coord),
        ' ', coord_name ($c_seq), "\n";
}
# This next line is here so it can be compared to the one above
use vars qw($output_header);
*output_header = \
'#Tdme coord coverage score1 score2 geo_gap   geo_gap geo_gap n_seq n_strct n_strct coord sequ
#     size                         quadratic linear  logi    gap   gap     widen   name  name
';

# ----------------------- do_seq_lib --------------------------------
# Take one sequence and align to a library (given by a list).
# After doing each alignment, call write_info which has a look at the
# quality of the model and writes out scores and gap information.
sub do_seq_lib ($ \@ \@ \@ $ \@ $ $ $ $)
{
    my ($seqname, $struct_list, $seq_dirs, $struct_dirs,
        $submat, $phd_dirs, $fx_params, $rs_prams,
        $ofile_sw, $ofile_nw) = @_;
    my ($seqfile, $c_seq, $seq, $phdfile, $sec_s_data);

    if ( ! ($seqfile = get_path (@$seq_dirs, "$seqname$bin_suffix"))) {
        print STDERR "$seqname$bin_suffix not found\n";
        return undef;
    }
    if ( ! ($c_seq = coord_read ($seqfile))) {
        print STDERR "Failed reading coordinates from $seqfile\n";
        return undef;
    }
    $seq = coord_get_seq ($c_seq);

    if ( ! ($phdfile = get_path (@$phd_dirs, "$seqname$phd_suffix"))) {
        print STDERR "$seqname$phd_suffix not found\n";
        return undef;
    }
    if ( ! ($sec_s_data = sec_s_data_read ($phdfile))) {
        print STDERR "Failed reading phd data from $phdfile\n";
        return undef;
    }

    my $tmp = "# " . coord_name ($c_seq) . "\n" .
              "# " . scalar @$struct_list .
              " lines follow after sequence size\n" .
              coord_size ($c_seq) . " residues \n".
              $output_header;
    print $ofile_sw $tmp;
    print $ofile_nw $tmp;
    undef ($tmp);

    for (my $i = 0; $i < @$struct_list; $i++) {
        my $structname = $$struct_list[$i];
        my $fullname   = "$structname$bin_suffix";
        my($structpath, $coord);
        my($sec_src_mat, $sim_scr_mat, $fx_scr_mat, $tot_scr_mat,$sec_scr_mat);
        if ( ! ($structpath = get_path (@$struct_dirs, $fullname))) {
            print STDERR "$fullname structure not found\n";
            return undef;
        }
        if (! ($coord = coord_read ( $structpath ))) {
            print STDERR "coord_open fail on $structpath\n";
            return undef;
        }
        my $c_size = coord_size ($coord);
        my $s_size = seq_size ($seq);
#       Create three types of score matrix and one for the total.
        $sec_scr_mat = score_mat_new( $s_size, $c_size);
        $sim_scr_mat = score_mat_new( $s_size, $c_size);
        $fx_scr_mat  = score_mat_new( $s_size, $c_size);
        
#       Now do each score and then add them together.
        score_sec  ($sec_scr_mat, $sec_s_data, $coord);
        score_smat ($sim_scr_mat, $seq, coord_get_seq ($coord), $submat);
        score_fx   ($fx_scr_mat, $seq, $coord, $fx_params);

        $tot_scr_mat = score_mat_add ($fx_scr_mat, $sim_scr_mat, $sim_scale);
        $tot_scr_mat = score_mat_add ($tot_scr_mat, $sec_scr_mat,$sec_s_scale);

        my $pair_set =
            score_mat_sum_sec (my $result_mat, $tot_scr_mat,
                               $coord, $sw1_sec_pnlty,
                               $sw1_pgap_open, $sw1_pgap_widen,
                               $sw1_qgap_open, $sw1_qgap_widen,
                               $S_AND_W);
        $pair_set =
            score_mat_sum_sec (my $result_mat2, $tot_scr_mat,
                               $coord, $sw2_sec_pnlty,
                               $sw2_pgap_open, $sw2_pgap_widen,
                               $sw2_qgap_open, $sw2_qgap_widen,
                               $S_AND_W, $pair_set);

        write_info( $pair_set, $seq, $seqname, $coord, $structname, $rs_prams,
                    $ofile_sw, $c_seq);
        $pair_set =
            score_mat_sum_sec (my $result_mat3, $tot_scr_mat,
                               $coord, $nw_sec_pnlty,
                               $nw_pgap_open, $nw_pgap_widen,
                               $nw_qgap_open, $nw_qgap_widen,
                               $N_AND_W, $pair_set );
        write_info( $pair_set, $seq, $seqname, $coord, $structname, $rs_prams,
                    $ofile_nw, $c_seq);
    }
    return 1;
}


# ----------------------- usage   -----------------------------------
sub usage ()
{
    print STDERR "$0 [options] seq_namefile struct_namefile outputfilename\n";
    exit (EXIT_FAILURE);
}

# ----------------------- mymain  -----------------------------------
sub mymain ()
{
#   Get the list of files with sequences.
    my $verbosity   = $DFLT_VERBOSITY;
    my @seq_dirs    = @DFLT_SEQ_DIRS;
    my @struct_dirs = @DFLT_STRUCT_DIRS;
    my @phd_dirs    = @DFLT_PHD_DIRS;
    my (@seq_list, @struct_list);
    my $fatalflag = 0;  # Sometimes we will set a flag, but die a
                        # bit later

    use Getopt::Std;
    use vars qw (%opts);
    if ( !getopts ('p:s:t:v:', \%opts)) {
        usage(); }
    if ( defined ( $opts { h }) || defined ( $opts { '?' })) {
        usage(); }
    if ( defined ( $opts { p })) {push (@phd_dirs,   split( ',', $opts { p }))}
    if ( defined ( $opts { s })) {push (@seq_dirs,   split( ',', $opts { s }))}
    if ( defined ( $opts { t })) {push (@struct_dirs,split( ',', $opts { t }))}
    if ( defined ( $opts { v })) { $verbosity = $opts { v }; }
    undef (%opts);

    check_dirs (@seq_dirs);
    check_dirs (@struct_dirs);
    check_dirs (@phd_dirs);

    if (@seq_dirs == 0) {
        die "\nNo valid sequence directories. Stopping.\n"; }
    if (@struct_dirs == 0) {
        die "\nNo valid structure directories. Stopping.\n"; }
    if (@phd_dirs == 0) {
        die "\nNo valid phd directories. Stopping.\n"; }

    if ( @ARGV != 4) {
        print STDERR "Wrong number of arguments\n";
        usage();
    }

    my $seq_filename    = $ARGV[0];
    my $struct_filename = $ARGV[1];
    my $output_1        = $ARGV[2];
    my $output_2        = $ARGV[3];

    print "$0 running with verbosity at $verbosity.\n",
          "Starting at ", scalar (localtime()), "\n";
    foreach my $i (@seq_dirs) {
        print "Looking in $i for sequences (from structures)\n";}
    foreach my $i (@struct_dirs) {
        print "Looking in $i for structures\n";}
    foreach my $i (@phd_dirs) {
        print "Looking in $i for phd files\n";}


    (@seq_list     = get_prot_list($seq_filename))    || $fatalflag++;
    (@struct_list  = get_prot_list($struct_filename)) || $fatalflag++;

#   We just got the (potentially) large of sequences and
#   structures to be aligned to. Now, before doing serious
#   calculations, just check if all the files exist. If there
#   is an error, we keep going for a while and stop after all
#   file checking is done. This is to get out as many errors
#   as possible in each pass.

    check_files (@seq_dirs, @seq_list, $bin_suffix)       && $fatalflag++;
    check_files (@struct_dirs, @struct_list, $bin_suffix) && $fatalflag++;
    check_files (@phd_dirs, @seq_list, $phd_suffix)       && $fatalflag++;
    if ($fatalflag) {
        die " Fatal errors"; }


    my $matname = "$MATRIX_DIR" . '/' . "$MATRIX_FILE";
    my $submat = sub_mat_read ($matname) || die "Fail on $matname";
    $submat = sub_mat_shift ($submat, $sim_mat_bottom);

    my ($tmp_p, $fx_params, $rs_prams);
    $tmp_p = "$PARAM_DIR/$FX9_PARAM_FILE";

    if (! -f $tmp_p) {
       die "Cannot find param file $tmp_p\n";}
    $fx_params = param_fx_read ($tmp_p) || die "Fail on fx9 param file $tmp_p";

    $tmp_p = "$PARAM_DIR/$RS_PARAM_FILE";
    $rs_prams = param_rs_read ($tmp_p) || die "Fail on rs param file $tmp_p";
    if ($verbosity > 1) {
        print "I have got\n",
        @seq_list, "sequences.\n",
        @struct_list, "structures in my library.\n"
    }

    my ($ofile_sw, $ofile_nw);
    open (OUT1, ">$output_1") ||
        die "Failed opening $output_1 for output: $!\n";
    $ofile_sw = \*OUT1;
    open (OUT2, ">$output_2") ||
        die "Failed opening $output_2 for output: $!\n";
    $ofile_nw = \*OUT2;


    my ($tmp1, $tmp);
    $tmp1 = $#seq_list + 1;
    $tmp = "# $tmp1 sequences\n";
    print $ofile_sw $tmp;
    print $ofile_nw $tmp;
    foreach my $s (@seq_list) {
        my $r = do_seq_lib ($s, @struct_list, @seq_dirs, @struct_dirs,
                            $submat, @phd_dirs, $fx_params, $rs_prams,
                            $ofile_sw, $ofile_nw);
        if ( ! $r ) {
            die "Error processing sequence $s\n"; }
    }
    close ($ofile_sw);
    close ($ofile_nw);
    print "\nGanz fertig at ", scalar (localtime()), "\n";
    my ($user, $system, $crap, $crap2) = times();
    printf "I took %d:%d min user and %.0f:%.0f min sys time\n",
    $user / 60, $user % 60, $system / 60, $system % 60;
    return EXIT_SUCCESS;
}
exit (mymain());
