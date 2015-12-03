package Wurst::Modeller;

use 5.006;
use strict;
use warnings;

require Exporter;
use AutoLoader qw(AUTOLOAD);

our @ISA = qw(Exporter);

# Items to export into callers namespace by default. Note: do not export
# names by default without a very good reason. Use EXPORT_OK instead.
# Do not simply export all your public functions/methods/constants.

# This allows declaration        use Wurst::Modeller ':all';
# If you do not need this, moving things directly into @EXPORT or @EXPORT_OK
# will save memory.
our %EXPORT_TAGS = ( 'all' => [ qw(Modeller_path Local_disk_path Get_model_parameters
                                   Write_Mshell
                                   model_top_script
                                   Pir_from_raw_model
                                   Modeller_model
                                   Modeller_m_model
                                   pir_mulpairs
                                   ) ] );

our @EXPORT_OK = ( @{ $EXPORT_TAGS{'all'} } );

our @EXPORT = qw( Modeller_model
                  Modeller_m_model
                  );

our $VERSION = '0.1';

# behaviours
our $DEBUG = 0;
my $TIDY_UP=1;
use File::Temp qw (tempdir tempfile);
use FileHandle;
use Wurst;

# -------------- MODELLER SHELL SCRIPTS ---------------------------------------
# Modeller is run via system(), in a temporary directory.
# MODINSTALL6v2 should point to your MODELLER installation.

# Shell script to start a modeller job
use vars qw ($MODEL_SHELL 
             $MODEL_DEFTOP
             $TOP_MAKEPDBPIR
             $PDB_FLAT_DIR_LOC
             );

*PDB_FLAT_DIR_LOC = \ "/projects/biodata/pdb/data/structures/all/pdb";
#use Env qw (MODINSTALL6v2);

# initialise the modeller paths
our $MODINSTALL6v2;
$SIG{TRAP} = 'IGNORE'; kill 'TRAP', $$;
# verify the installation *really* exists
my @modvars = ( 'MODINSTALL6v2',
                'EXECUTABLE_TYPE6v2',
                'LIBS_LIB6v2',
                'KEY_MODELLER6v2');
my $LOCAL_MODINSTALL = 1;
foreach my $m_var (@modvars) {
    if (exists($ENV{$m_var})) {
        $LOCAL_MODINSTALL = $LOCAL_MODINSTALL & defined($ENV{$m_var});
    } else { 
        $LOCAL_MODINSTALL = 0;
    }
}
if ($LOCAL_MODINSTALL) {
    # assume everything is set up correctly
    # so only make sure we use the right paths
    $MODINSTALL6v2 = $ENV{MODINSTALL6v2};
    $MODEL_SHELL = '"#!/usr/bin/tcsh
setenv MODINSTALL6v2      ".$MODINSTALL6v2."
setenv EXECUTABLE_TYPE6v2 i386-absoft
setenv LIBS_LIB6v2        ".$MODINSTALL6v2."/modlib/libs.lib
setenv KEY_MODELLER6v2    MODELIRANJE
#alias  mod                mod6v2
set    path=(\$path ".$MODINSTALL6v2."/bin)
limit  stacksize unlimited
cd $TEMP_DIR
mod6v2 $TEMP_MODSCRIPT
"';
    ;
} else {
    # find our 'own' modeller installation
    $MODINSTALL6v2 = "/home/procter/bin/modeller6v2/";
    $MODEL_SHELL = '"#!/usr/bin/tcsh
setenv MODINSTALL6v2      ".$MODINSTALL6v2."
setenv EXECUTABLE_TYPE6v2 i386-absoft
setenv LIBS_LIB6v2        ".$MODINSTALL6v2."/modlib/libs.lib
setenv KEY_MODELLER6v2    MODELIRANJE
alias  mod                mod6v2
set    path=(\$path ".$MODINSTALL6v2."/bin)
limit  stacksize unlimited
cd $TEMP_DIR
mod6v2 $TEMP_MODSCRIPT
"';
}

sub Write_Mshell( $ @; % ) {
    my ($TEMP_DIR, $modscripts, $parm) = @_;
    my $TEMP_MODSCRIPT=shift @{$modscripts};
    my ($tf) = ($TEMP_MODSCRIPT=~m/(.+).top/);
    my ($fh, $tempfile) = tempfile(TEMPLATE=>$tf."XXXXXXX", SUFFIX=>'.sh', DIR=>$TEMP_DIR);

    $fh->print(eval($MODEL_SHELL));
    while (scalar @{$modscripts}) {
        $fh->print("mod6v2 ".(shift @{$modscripts})."\n");
    }
    close $fh;
    return $tempfile;
}

# ----------- TOP SCRIPT GENERATION ----------------------------------------
# .top scripts are run by MODELLER
# they tell it what to do, which alignments to use, and which coordinates
# correspond to which sequence.
# We also put in some helpful hints for further refinement
my $MODEL_TOPREMARKS = '"
# These can all be varied to make modeller work harder
# Defaults are given, but commented out. These definitions
# come from __defs.top under ".$MODINSTALL6v2."/bin

# SET LIBRARY_SCHEDULE = 4
# 1 ... thorough var target func 
# 4 ... faster var target func schedule

# maximal numb of iterations for the cycles
# the variable target function method
#SET MAX_VAR_ITERATIONS = 200  

# MD_LEVEL and LOOP_MD_LEVEL define
# what kind of optimization is done 
# after the variable target function method:
# \'nothing\'            ... nothing;
# \'refine_1\'            ... very fast MD annealing;
# \'refine_2\'            ... fast MD annealing;
# \'refine_3\'            ... slow MD annealing;
# \'refine_4\'            ... very slow MD annealing;
# \'refine_5\'            ... very slow/large dt MD annealing;


SET RSTRS_REFINED = 1         # the types of restraints used to define
# hot spots when MD_LEVEL <> \'nothing\':
# 0 ... stereochemistry only;
# 1 ... stereochemistry and dihedral;
# 2 ... all restraints;

SET EXIT_STAGE    = 0         # 0 ... no effect;
# 1 ... exit without any optimization after
#       restraints and an initial model are 
#       calculated (more efficient than
#       REPEAT_OPTIMIZATION=0);
# 2 ... exit after the initial model is calculated
#       (restraints are not calculated)

SET REPEAT_OPTIMIZATION = 1   # how many times the whole optimization
# schedule (variable target function
# method and refinement) is repeated
# for each initial model;

SET TRACE_OUTPUT =  10        # every which CG or MD cycle is reported;

SET MAX_MOLPDF   = 100E3      # abort optimization of the current model if 
# the molecular pdf is larger than this and
# continue with the next model;
SET DEVIATION     = 4.0  # the amount of randomization of the initial model
                         # must be > 0 if different final models are wanted;

"';
$MODEL_DEFTOP = '"# Default model generation for $SEQ_NAME from $Template
INCLUDE
# directory with input atom files
SET ATOM_FILES_DIRECTORY = \'./:$PDB_FLAT_DIR_LOC\' 
# input file with templates and alignment to model sequence
SET ALNFILE     = \'$algfile\'             
# Template as Full pdb or just backbone and Cbeta chain
SET KNOWNS      = \'$Template_Full\' \'$Template\' 
# tag for sequence within ALNFILE
SET SEQUENCE    = \'$MDL_SEQ_NAME\' 
# you can probably set some of these 1\'s to 0 to make the logfile smaller
SET OUTPUT_CONTROL = 1 1 1 1 0 

# this is the range of basic models, they will be called
# $MDL_SEQ_NAME.B9999(0+)($Starting_Model .. $Ending_Model)
SET STARTING_MODEL = $Starting_Model
SET ENDING_MODEL = $Ending_Model


# If this flag is set to 1 then loop modelling will also be done
SET DO_LOOPS = $Do_Loops

# if we do loops then additional models will be made :
# $MDL_SEQ_NAME.BL(0+)($Starting_Model .. $Ending_Model)($LM_Start .. $LM_End)
# $MDL_SEQ_NAME (base model)(loop model num from this base model)
SET LOOP_STARTING_MODEL = $LM_Start
SET LOOP_ENDING_MODEL = $LM_End

# if you want to do really careful loop modelling then you
# need to pick the regions yourself
# by default, modeller just picks anywhere where there are 
# unaligned bits of sequence
#  SUBROUTINE ROUTINE = \'select_loop_atoms\'
#   # Make sure your atom selection is NOT a function of coordinates!!!!
#   # (because atom selections for restraints construction and optimization
#   #  will be different!) This means SELECTION_SEGMENT = \'SPHERE\' is bad!
# 
#   # 8 residue insertion:
#    PICK_ATOMS SELECTION_SEGMENT = \'28:A\' \'35:A\', SELECTION_STATUS = \'initialize\'
#   # Add 4 residue insertion:
#    PICK_ATOMS SELECTION_SEGMENT = \'173:A\' \'176:A\', SELECTION_STATUS = \'add\'
#
#
#   RETURN
# END_SUBROUTINE
# There are much more things you can do with the loop modelling procedure
# look at ".$MODINSTALL6v2."/bin/__loop.top for more info

# see below for the types of optimization schemes that can be set here
SET LOOP_MD_LEVEL = \'$Loop_MD\'
SET MD_LEVEL = \'$Model_MD\'
SET REFINE_HOT_ONLY = $Hot_refine_only
# 1 ... select and optimize only HOT atoms in refine;
# 0 ... select and optimize all atoms in refine;
# usually about half of the atoms are hot; in such cases,
# 0 is faster for sequences longer than about 100 aa
# because a faster non-bonded pairs algorithm can be used.

CALL ROUTINE    = \'model\'                # get a model

# Now, you might want to play with setting the parameters below
# to make better models, but remember to change the STARTING_MODEL
# and ENDING_MODEL range, if you want to keep your current models.
"';
sub model_top_script( $ $ $ $ $ % ) {
    # no defaults here
    my $SEQ_NAME = shift @_; # sequence name
    my $Template = shift @_; # Name of the template used to model - also names the template 
    #               alignment to the sequence
    my $algfile = shift @_;;
    $algfile =~ s/.+[^\\]\///; # strips
    my $Template_Full = shift @_; # name for pir entry for full pdb file
    my $MDL_SEQ_NAME = shift @_; # safe name of sequence for a tag in the PIR file
    my $ph = shift @_;
    # essentially defaults
    sub xod( $ $ % ) {
        my ($dval, $k, $hsh) = @_;
        if (defined($hsh)
            and (exists($hsh->{$k})) and (defined($hsh->{$k}))) {
            return $hsh->{$k};
        }
        return $dval;
    }
    my ($Starting_Model,$Ending_Model) = 
        (xod(1,"STARTING_MODEL", $ph),
        xod(1,"ENDING_MODEL", $ph));
    my $Model_MD = xod("refine_1",
                       "MD_LEVEL", $ph); # default MD level
    my $Loop_MD = xod("refine_3",
                       "LOOP_MD_LEVEL", $ph); # default loop MD level
    my ($Do_Loops, $LM_Start, $LM_End) = 
        (xod(0,"DO_LOOPS", $ph),
         xod(1,"LOOP_MODEL_START",$ph),
         xod(1,"LOOP_MODEL_END",$ph)); # default loop settings
    my ($Hot_refine_only) = xod(0,"REFINE_HOT_ONLY", $ph);
    my $script = eval($MODEL_DEFTOP).
        "\n".eval($MODEL_TOPREMARKS); # put all the variables in place;
    # sometimes we may not have a defined Full_Template model
    $script =~ s/''//g; 
    return $script;
}


$TOP_MAKEPDBPIR = '"
# expand an alignment via Modeller
# this one just adds the original PDB file 
# entry to a template/sequence alignment

INCLUDE

# directory with input atom files
SET ATOM_FILES_DIRECTORY = \'./:$PDB_FLAT_DIR_LOC\' 

# input file with templates and alignment to model sequence

READ_ALIGNMENT FILE=\'$algfile\'             

READ_MODEL FILE = \'pdb$Template_PDBID\' , MODEL_SEGMENT=\'FIRST:$Template_CHAIN\' \'LAST:$Template_CHAIN\'

SEQUENCE_TO_ALI ADD_SEQUENCE=on, ALIGN_CODES = ALIGN_CODES \'pdb$Template_PDBID\' , ATOM_FILES = ALIGN_CODES
ALIGN ALIGN_BLOCK = $Align_block_last , ALIGN_WHAT=\'LAST\', GAP_PENALTIES_1D=-0.1 -0.1 

WRITE_ALIGNMENT FILE=\'$pdbalgfile\'
"';

sub Make_pdb2seqtop ( $ $ $; $ ) {
    my ($Pdb_file, $algfile, $pdbalgfile, $Align_block_last) = @_;
    defined($Align_block_last) or $Align_block_last = 2;
    (defined($algfile) and (length $algfile > 4)) 
        or die ("Script error. Make_pdb2seqtop needs an alignment to work on\n");
    my ($Template_PDBID, $Template_CHAIN) = $Pdb_file=~/([1-9][A-Za-z0-9]{3})([A-Z0-9_])/;
    ($Template_CHAIN eq '_') and ($Template_CHAIN=' ');
    # remove any absolute filename refs
    $algfile =~ s/.+[^\\]\///; # strips
    $pdbalgfile =~ s/.+[^\\]\///; # strips
    my $script = eval($TOP_MAKEPDBPIR);
    ($@) and warn ("Possible problems with TOP_MAKEPDBPIR translation:", $@);
    return ("pdb".$Template_PDBID, $script);
}

# another script would be :
# pick a couple of good models, read them in, make any obvious restraints
# minimise/simulate and writeout.

# ------------------------ PIR FILE GENERATION ----------------------------------
# The simplest PIR file generator
# returns model.pdb coord and PIR as a string
sub Pir_from_raw_model ( $ $ $ $ $ $ ; $ ) {
    my ($snam, $seq, $tnam, $tstrn, $templ, $pairset, $temp_model) = @_;
    my @p;
    my $t_model_seq;
    my ($s_aligned, $t_aligned) = pair_set_coverage($pairset, seq_size($seq), coord_size($templ));

    my $t_model = (defined($temp_model)) ? $temp_model : (make_model($pairset, $seq, $templ));
    my @txt_algnment = split "\n",(pair_set_string($pairset, $seq, coord_get_seq($templ)));
    my ($al_seq, $al_tem, $sseq, @sq);
    $sseq = seq_print($seq);
    $sseq =~ s/>.+\s+//;
    $sseq =~ s/\s//g;
    chomp $sseq;
    @sq = split '',$sseq;
    my $tseq = seq_print(coord_get_seq($templ));
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
            if (defined($tins)) { $al_seq.=("-"x(length $tins));
                                  $al_tem.=$tins; }
            if (defined($sins)) { $al_seq.=$sins;
                                  $al_tem.=("-"x(length $sins));
                              };
        };
    }

    # this case is really easy but we do have to use the original template model
    my ($m_start, $m_end, $a_start, $a_end) = 
        (model_pdb_num($templ, 0),model_pdb_num($templ, -1+coord_size($templ)),
         1, seq_size($seq));
    # Sequence comes first
    push @p, ">P1; $snam";
    push @p, "sequence:$a_start:1: :$a_end:\@:.:.:.:.";
    push @p, "".(uc $al_seq)."*","";
    
    # template sequence last so we can use ALIGN_WHAT=last in top script
    # to expand alignment to whole PDB file
    push @p, ">P1; $tnam";
    push @p, "structure:$tstrn:$m_start:\@:$m_end:\@:.:.:.:.";
    push @p, "".(uc $al_tem)."*","";
    return ($templ, (join "\n", @p,""));
}

# Analyse a modeller log file in several different ways
# return codes
# -1 : total failure
# 0 : Warnings only
# remains : count of error lines in file...
#
my @non_serious_e_list = ('Numbers of expected, actual arguments', 'VIOL_REPORT_CUT');

sub inspect_modeller_log( $; % ) {
    my ($logf, $parms) = @_;
    my $line;
    (-f $logf) or return (-1, {E=>["Modeller process failed. No logfile.\n"]});
    my $fh = new FileHandle;
    my %msgs;
    $fh->open($logf) or die $@."\nProblems opening $logf";
    do {
        $line=<$fh>;
        my ($e,$m) = ($line=~/[^>](.)>(.+)/);
        if ($e) {
#            if ($line=~/  need to check for non-seriousity of the error
            (exists ($msgs{$e})) or ($msgs{$e}=[]);
            push @{$msgs{$e}}, $m;
        }
    } while (not $fh->eof);
    $fh->close;
    return (scalar @{(exists($msgs{E}) ? $msgs{E} : [] )}, {(%msgs)});
}


# utilities

sub make_local_dir ( % ) {
    my $pset = shift @_;
    if (not (exists ($$pset{localdisk}) and (defined $$pset{localdisk}) and (-d $$pset{localdisk}))) {
        if (-d "/local") { $$pset{localdisk} = "/local"; }
        elsif (-d "/tmp") { $$pset{localdisk} = "/tmp"; };
    }
    my $tempdir = tempdir(TEMPLATE=>"modellerXXXXXXXX", DIR=>$$pset{localdisk}, CLEANUP=>$TIDY_UP);
    if (not -d $tempdir) { die ("Broke trying to make a tempdir in $$pset{localdisk}"); };
    return $tempdir;
};


# The main model function

sub pdbfile_exists ( $ ) {
    my $pd = shift @_;
    if (-f $PDB_FLAT_DIR_LOC."/pdb".(substr $pd, 0, 4).".ent") { return 1; };
    if (-f $PDB_FLAT_DIR_LOC."/pdb".(substr $pd, 0, 4).".ent.Z") { 
        return 1; };
    if (-f $PDB_FLAT_DIR_LOC."/pdb".(substr $pd, 0, 4).".ent.gz") {
        print "Warning. Can only find $pd as a .gz file (Modeller doesn't know about these!\n";
        return 1; };
    return undef;
}


sub Modeller_model( $ $ $ $ $ $ $; %) {
    my ($s_nam, $seq, $t_nam, $alignment, $t_wmdl, $t_chain, $t_pdb, $xtra_params) = @_;
    my $t_seq = coord_get_seq($t_chain);
    my %model_params; # = m_params_merge($xtra_params); # how !
    my $pnam;
    foreach $pnam (keys %{$xtra_params}) {
        $model_params{$pnam} = $$xtra_params{$pnam};
    }
    my $wk_dir = make_local_dir(\%model_params);
    ($DEBUG) and print "wd:\n$wk_dir\n";
    my $root_fname = "mdl_".$s_nam."_f_".$t_nam;
    my $mini_pirfile = "$wk_dir/".$root_fname."_simple.pir";
    my $mini_xlattop = "$wk_dir/pdb_add.top";
    my $logf_xlattop = "$wk_dir/pdb_add.log";
    my $full_pirfile = "$wk_dir/$root_fname.pir";
    my $full_topfile = "$wk_dir/$root_fname.top";
    my $logf_topfile = "$wk_dir/$root_fname.log";
    my @gen_models;
    my @w;
    # Do we try to use original pdb file ?
    if (defined($t_pdb) and (pdbfile_exists($t_pdb))) {
        # make the mini_pirfile

        (open PIRFILE, ">$mini_pirfile") or
            die ("Modeller.pm script error: Couldn't open $mini_pirfile\n");
        
        @w = Pir_from_raw_model($s_nam, $seq, 
                                $t_nam, $t_pdb, $t_chain, 
                                $alignment, $t_wmdl);
        print PIRFILE $w[1];
        close PIRFILE;
    } else {
        $mini_pirfile = "";
        # make the full_pirfile and skip the translation

        (open PIRFILE, ">$full_pirfile") or
            die ("Modeller.pm script error: Couldn't open $full_pirfile\n");
        
        @w = Pir_from_raw_model($s_nam, $seq, $t_nam, $t_pdb, 
                                $t_chain, $alignment, $t_wmdl);
        print PIRFILE $w[1];
        close PIRFILE;

    }

    coord_2_pdb("$wk_dir/$t_pdb.atm",$w[0]); # this is the template atom file mentioned in the pir
    # generate the TOP file
    my @topfs;
    (open TOPFILE, ">$full_topfile") or 
        die ("Modeller.pm script error: Couldn't open $full_topfile\n");

    if ($mini_pirfile) {

        open MTOPFILE, ">$mini_xlattop" or
            die ("Modeller.pm script error: Couldn't open $mini_xlattop\n");
        my ($pdb_tag, $script) = 
            Make_pdb2seqtop( $t_pdb, 
                             $mini_pirfile, $full_pirfile);
        print MTOPFILE $script;
        close MTOPFILE;
        print TOPFILE 
            model_top_script($s_nam, $pdb_tag, $full_pirfile, 
                             '', $s_nam, \%model_params);
        push @topfs, $mini_xlattop, $full_topfile;

    } else {
        print TOPFILE 
            model_top_script($s_nam, $t_nam, $full_pirfile, 
                             '', $s_nam, \%model_params);
        push @topfs, $full_topfile;
    }
    close TOPFILE;
    
    # run the script as a new thread ?
    my $shell_scr = Write_Mshell( $wk_dir, \@topfs, %model_params);
    my $cwd = $ENV{PWD};
    system("tcsh",$shell_scr);
    chdir($cwd);
    # wait for the script to finish and check on the output
    my @model_status;
    if ($mini_pirfile) {
        @model_status = inspect_modeller_log($logf_xlattop);
        if ($model_status[0]==0) {
            @model_status = inspect_modeller_log($logf_topfile);
        } else {
            print "Problems with initial alignment of pdb chain to template chain\n";
        }
    } else { @model_status = inspect_modeller_log($logf_topfile); };

    if ($model_status[0]==0) {
        # we always copy the models if they exist.
        my @ms = (glob("$wk_dir/$s_nam.B*"));
        foreach my $i (@ms) {
            my $newnam = $i;
            my $ts = $s_nam."_".$t_nam;
            $newnam=~s!/$s_nam!/$ts!x;
            move($i, $newnam);
            use File::Copy;
            $DEBUG and print "Moving $newnam to ./\n";
            move($newnam, "./");
            push @gen_models, $newnam;
        }
        @gen_models = map { $_ =~s/$wk_dir\///g; $_ } @gen_models;
        
        
        # No terminating errors so now extract some more info
        # from the models
#        my @model_set = the filenames;
#        my @model_lset = the loop model filenames;
        # hash of info
        # the OBJECTIVE function values
        # ? nothing more ?
    } else {
        print "".(join "\n",@{$model_status[1]->{E}}, "\n");
        # we always copy the models if they exist.
        my @ms = (glob("$wk_dir/$s_nam.B*"));
        foreach my $i (@ms) {
            my $newnam = $i;
            my $ts = $s_nam."_".$t_nam;
            $newnam=~s!/$s_nam!/$ts!x;
            move($i, $newnam);
            use File::Copy;
            $DEBUG and print "Moving $newnam to ./ (despite errors)\n";
            move($newnam, "./");
            push @gen_models, $newnam;
        }
        @gen_models = map { $_ =~s/$wk_dir\///g; $_ } @gen_models;
        # somehow we have to bail out here or restart with a better configuration
        # generate the summary of what broke
        # write it if necessary
        # return failure code
    }
# keep the files for the mo.
    # make a sensible name
    system("mkdir $wk_dir/".$s_nam."_$t_nam");
    system("cp $wk_dir/*.* $wk_dir/".$s_nam."_$t_nam/"); # only work for unix!
    system("tar -cf - -C $wk_dir "."./".$s_nam."_$t_nam \| gzip >./$root_fname.tgz");
    
    return @gen_models; # the model set if they exist
}

# need an expansion function
# given two alignment strings creates a new stack of alignments based on another pair
# containing at least one common sequence
    
sub modeller_pir ( $ $ $; $ $ ) {
    my ($pair_set, $sid, $seq, $snam);
}

sub pick_sequence_insertions ( $ ) {
    my ($pair_set) = @_; 
};
# given a set of pick regions, rationalise them to contiguous segments
# and generate the TOP code
sub create_loop_picks( @ ) {
    
}

# Multiple alignment of many templates to one sequence
# pdb2wurst_template magic can only be carried out on a pir
# where the template structure entry is the last in the file.
# This bit of perl reorders a given PIR tag entry to make it the last one.
my $Reorder_pir_pl = 'perl -e \'my (@p,@m,$t, $w); $t = shift @ARGV; $w = 0; while (<>) { ($_=~/>P1/) and (($_=~/$t/) ? ($w=1) : ($w=2)); ($w==1) and push @mtag, $_; ($w==2) and push @pir, $_; ($_ =~ /\*\s*/) and $w=0;} print join "",@pir, @mtag;\'';
# expands to a perl inline to rearrange a pir
sub reorder_pir_com ( $ $ $) {
    my ($pir_filen, $tag, $new_f) =@_;
    return ($Reorder_pir_pl." $tag $pir_filen > $new_f");
};

# Generation of a multiple alignment from several pairwise alignments

# pair_q
# returns true if any of the arrays in an array of array refs 
# have any entries.
sub pair_q {
    my @pairs = @_;
    
    while (scalar @pairs) {
        if (scalar @{(shift @pairs)}) {
            return 1;
        }
    }
    return 0;
}

# given sequence, array of structures, corresponding
# array of pair_set objects.

sub pir_mulpairs ( $ @ @ ) {
    # Compute coverage strings for all pairsets.
    # These are all registered in terms of sequence coverage string.
    # create block alignment from walking down coverage strings.

    my ($seq, $crds, $pairs) = @_;
    my $n_tplt = scalar @{$crds};
    my (@msa, $bp);
    my (@seqs, @covset);
    my $p;

    for ($p=0; $p<scalar @{$crds}; $p++) {
        my ($a,$b) = pair_set_coverage( $$pairs[$p], seq_size($seq), 
                                        seq_size(coord_get_seq($$crds[$p])));
        $covset[$p*2] = [ split "", $a ];
        $covset[$p*2+1] = [ split "",$b ];
        $msa[$p] = [];
        $seqs[$p] = [ split "", seq_print(coord_get_seq($$crds[$p])) ];
        pop @{$seqs[$p]};
    }
    $msa[$n_tplt] = []; # sequence entry
    $seqs[$n_tplt] = [ (split "", seq_print($seq))];
    pop @{$seqs[$n_tplt]};
    my $q;
    while (scalar @{$seqs[$n_tplt]}) {
        my $ins = 0; # length of really aligned row region
        my @al_atrow = []; # no seq->str aligned in this column
        # put in all structure gaps first, then any sequence-gap,
        # and finally put in an aligned column
        for ($p=0; $p<$n_tplt; $p++) {
            my $i = 0;
            if (scalar @{$covset[$p*2+1]}) {
                while((scalar @{$covset[$p*2+1]})
                      and ($covset[$p*2+1]->[0] eq '0')) {
                    shift @{$covset[$p*2+1]};
                    push @{$msa[$p]}, (shift @{$seqs[$p]});
                    for ($q=0; $q<=$n_tplt; $q++) {
                        ($q==$p) or (push @{$msa[$q]}, '-');
                    }
                }
            }
            if (scalar @{$covset[$p*2]}) {
                if (($covset[$p*2]->[0] eq '0')
                    or (not scalar @{$seqs[$p]})) {
                    $al_atrow[$p] = '-';
                } else {
                    shift @{$covset[$p*2+1]};
                    $al_atrow[$p] = shift @{$seqs[$p]};
                }
                shift @{$covset[$p*2]};
            } else {
                $al_atrow[$p] = '-';
            }
        }
        # make the tail of the directed star
        foreach ($p=0; $p<=$n_tplt; $p++) {
            if ($p<$n_tplt) {
                push @{$msa[$p]}, shift @al_atrow;
            } else {
                push @{$msa[$p]}, shift @{$seqs[$p]};
            }
        }
    }
    # final tails added at end
    $p=0;
    while ($p<$n_tplt) {
        if (scalar @{$seqs[$p]}) {
            push @{$msa[$p]}, (@{$seqs[$p]});
            for ($q=0; $q<=$n_tplt; $q++) {
                ($q!=$p) and (push @{$msa[$q]}, ('-'x(scalar @{$seqs[$p]})));
            }
        }
        $p++;
    }
    for ($p=0; $p<=$n_tplt; $p++) {
        $seqs[$p] = join "",@{$msa[$p]};
    }
    if ($DEBUG) { 
        print "Msa (debug)\n";
        $p=0;
        while ((72*$p)<length($seqs[0])) {
            for ($q=0; $q<=$n_tplt; $q++) {
                printf "%2d : %s\n",$q+1,(substr $seqs[$q], $p*72, 72);
            }
            $p++;
        }
        print "---\n";
    }
    return (@seqs);
}
# snam, sequence, vector of names, coords,  and pairsets
sub Pir_from_manymodels ($ $ $ @ @ @; $) {
    my ($snam, $comment, $seq, $names, $coord, $pairs) = @_;
    
    # first, make the msa
    my @msa = pir_mulpairs ( $seq, $coord, $pairs);
    my $seq_al = pop @msa; # sequence line is last entry.

    my @pir;
    my @pirtags;

    # Sequence comes first
    push @pir, "# $comment", "# Wurst::Modeller generated ".localtime(time);
    push @pir, ">P1; $snam";
    push @pir, "sequence:$snam:1: :".(seq_size($seq)).":\@:.:.:.:.";
    push @pir, "".(uc $seq_al)."*","";

    # Generate PIR tag entries (in reverse)
    my $t = 0;
    do {
        my $templ = $$coord[$t];
        my $tnam = $$names[$t];
        my $tstrn = $tnam; # for wurst generated chain coords
        # pdb entry lengths
        my ($m_start, $m_end)
            = (model_pdb_num($templ, 0),model_pdb_num($templ, -1+coord_size($templ)));
        push @pir, ">P1; $tnam";
        push @pir, "structure:$tstrn:$m_start:\@:$m_end:\@:.:.:.:.";
        push @pir, "".(uc $msa[$t])."*","";
        unshift @pirtags, $tnam;
    } while (++$t<scalar @{$coord});
    
    return (\@pirtags, (join "\n", @pir,""));
}

# shell operations
# $mini_pir for first tag."\nperl -e rearrange > new pir"
# $mini_pir.new for second "\n etc


# Takes a name, a sequence, a descriptive comment,
# a vector of template chain pdb_id's,
# a vector of coord objects
# a vector of alignments from each template to sequence,
# optional extra hash of parameters to control modeller.

sub Modeller_m_model( $ $ $ @ @ @; %) {
    my ($s_nam, $seq, $comment, 
        $ts_nm, $t_chains, 
        $alignments, 
        $xtra_params) = @_;
    my @t_seq = map { coord_get_seq($_); } @{$t_chains};

    my %model_params; # Merge any extra paramters into here
    my $pnam;
    foreach $pnam (keys %{$xtra_params}) {
        $model_params{$pnam} = $$xtra_params{$pnam};
    }

    my $wk_dir = make_local_dir(\%model_params);
    ($DEBUG) and print "wd:\n$wk_dir\n";

    my $root_fname = "mdl_".$s_nam; # hope this is unique'ish
    # these are the files for the final modelling procedure
    my $full_pirfile = "$root_fname.pir";
    my $full_topfile = "$root_fname.top";
    my $logf_topfile = "$root_fname.log";

    # these are all stubs used to generate $full_pirfile
    # these get substituted each time
    my $mini_pirfile = "$root_fname";
    my $mini_xlattop = '"pdb_add_".$pdbid.".top"'; 
    my $logf_xlattop = '"pdb_add_".$pdbid.".log"';

    my @gen_models; # final array of generated models

    
    my @w;
    my ($t,$pdbid);
    my @x_pdb;

    # write basic pir msa
    (open PIRFILE, ">$wk_dir/".$mini_pirfile."_ini.pir") or
        die
        ("Modeller.pm script error: Couldn't open $mini_pirfile"."_ini.pir\n");
    
    # here, all tags correspond to the wurst chain model
    my ($tags, $pir) = 
        Pir_from_manymodels( $s_nam, $comment, $seq,
                             $ts_nm, 
                             $t_chains, $alignments);
    print PIRFILE $pir."\n";
    close PIRFILE;
    my @sh_coms; # the command pipeline
    my ($o_pir, $curp, $newp) = ($mini_pirfile."_ini.pir",
                                 0, 1); 
    # rearrange is after first transform, and final reorders the sequence
    my @r_tags = (@{$tags}, $s_nam);
    shift @r_tags; 
    my @newtags; # the tags for the pdb_chainid entries
    my @xlogfls;
    my $al_block = 1+scalar (@r_tags); # this is to define the block of alignments
    # for the transformation top.
    foreach $t (@{$tags}) { 
        if (not pdbfile_exists($t)) {
            warn "Can't find pdb $t in any of the known pdb locations.\n";
        }
        # generate top to transform the alignment for this pdb file
        $pdbid = $t;
        (open MTOPFILE, ">$wk_dir/".(eval($mini_xlattop))) or
            die ("Modeller.pm script error: Couldn't open $mini_xlattop".
                 (eval ($mini_xlattop))."\n");
        
        my ($pdb_tag, $script) = 
            Make_pdb2seqtop( $t, 
                             $o_pir, $mini_pirfile.$curp.".pir", $al_block++ ); # always goes up for next
        print MTOPFILE $script;
        close MTOPFILE;
        push @newtags, $pdb_tag;
        push @xlogfls, "$wk_dir/".eval($logf_xlattop);
        push @sh_coms, eval($mini_xlattop)."\n".
            (reorder_pir_com($mini_pirfile.$curp.".pir",
                             shift @r_tags,
                             $mini_pirfile.$curp."r.pir"));
        $o_pir = $mini_pirfile.$curp."r.pir";
        $curp++;
    }
    # add the cp to the last call to make the final pir file 
    $sh_coms[$#sh_coms] = $sh_coms[$#sh_coms]."\ncp $o_pir $full_pirfile";
    
    # and make the multi-model topfile
    # we are really lazy here and just use the original
    # model_top_script code, with a specially formed pdb_tag
    
    (open TOPFILE, ">$wk_dir/$full_topfile") or 
        die ("Modeller.pm script error: Couldn't open $full_topfile\n");
    print TOPFILE "# Wurst::Modeller::Modeller_m_model generated top file\n";
    print TOPFILE "# generated ".(localtime(time))."\n";
    print TOPFILE "# $comment\n";
    print TOPFILE 
        model_top_script($s_nam, (join "\' \'",@newtags),
                         $full_pirfile, 
                         '', $s_nam, \%model_params);
    push @sh_coms, $full_topfile; # complete the pipeline
    
    # run the pipeline
    (not exists $model_params{NO_SHELL_RUN}) and ($model_params{NO_SHELL_RUN} = 0);
    my $shell_scr = Write_Mshell( $wk_dir, \@sh_coms, %model_params);
    my $cwd = $ENV{PWD};
    if (not $model_params{NO_SHELL_RUN}) {
        # only do the checking if we really run the job
        system("tcsh",$shell_scr);
        chdir($cwd);
        # wait for the script to finish and check on the output
        my @model_status;
        # check each translation went fine
        my $mr = 1; # count of number of modeller commands
        my $e_mr = 0; # error at stage
        while (($mr) and (scalar @xlogfls)) {
            my $mf = shift @xlogfls;
            @model_status = inspect_modeller_log($mf);
            if ($model_status[0]!=0) {
                print "Problems with an alignment of pdb chain to template chain\n".
                    "Log file is $mf from run $mr\n";
                $e_mr = $mr;
                $mr = 0;
            } else {
                $mr++;
            }
        }
        if ($mr) {
            @model_status = inspect_modeller_log("$wk_dir/".$logf_topfile);
        }
        # we always copy the models if they exist.
        # we also do very little name translation
        my @ms = (glob("$wk_dir/$s_nam.B*"));
        foreach my $i (@ms) {
            my $newnam = $i;
            use File::Copy;
            $DEBUG and print "Moving $newnam to ./\n";
            move($newnam, "./");
            push @gen_models, $newnam;
        }
        @gen_models = map { $_ =~s/$wk_dir\///g; $_ } @gen_models;

        if ($model_status[0]!=0) {
            print "There were errors in ".
                (($mr>(scalar @newtags))?"Main modelbuilding":"Translation stage $e_mr").
                "\n";
            print "".(join "\n",@{$model_status[1]->{E}}, "\n");
        }
    }
    # keep the input files.
    # make a sensible name for the tarball
    
    system("mkdir $wk_dir/".$s_nam."_multi");
    system("cp $wk_dir/*.* $wk_dir/".$s_nam."_multi/"); # only work for unix!
    # move old files out the way
    if (-f "./$root_fname.tgz") {
        my $i=0;
        while (-f "./$root_fname.tgz.$i") {
            $i++;
        }
        move("./$root_fname.tgz", "./$root_fname.tgz.$i");
    }
    system("tar -cf - -C $wk_dir "."./".$s_nam."_multi \| gzip >./$root_fname.tgz");
    
    return @gen_models; # the model set if they exist
}

1; 
# -----------------excess comments area-----------------------
# Overall resources
# 1. Modeller installation and running script
#  - locate model template as a PDB file
#    . Also use -CB model as a template to ensure that
#    . altloc fuckups don't happen (perhaps)
#  - Create temp-directory (somewhere local to node)
#  - (cpoy)write a TOP file
#    . default options
#    . ways to work harder - refinement & loops
#    . multiple templates ?
#    . multiple runs - number of models
#  - spawn a modeller process
#    . ideally this should be through globus but 
#      we just assume we are a single threaded process anyhow.
#    . wait for finish (timelimits for too long a modelling process)
#    . analyse result
#        Find structures - know how many .B/.BL structures we should
#        have made.
#        .log file
#        E> - deal with a broken modeller process
#        Success - make a summary of energy/etc as a return value
#    . Tidy up
#      . decide what to be done with files - keep logs etc
#      . package up all supporting files
#      . remove tempdir
#      . return with appropriate info (where to get the models, etc)


__END__
=head1 NAME

Wurst::Modeller - MODELLER jobs with wurst Alignments

=head1 SYNOPSIS

  use Wurst::Modeller;
  my @models = Modeller_model ( $sequence_name, $sequence,
                                $template_name, $pair_set,
                                $template_model, $template_coord,
                                $t_pdb_name; %extra_params );
  
=head1 ABSTRACT

  Tries to make beautiful proteins with modeller. Forks a job with 
  a generic set of parameters, based on the C-B model of the 
  sequence aligned onto the template.
  

=head1 DESCRIPTION

There is very little to add to from the abstract. The module,
at use Wurst::Modeller time, tries to find a real installation
of modeller, but is easily fooled. Either unset all variables
if modeller doesn't run, so the module will use an installation
that it knows about, or make sure that the variables below
are properly set in the shell environment :

  MODINSTALL6v2      /path/to/my/dir/that/contains/modeller6v2
  EXECUTABLE_TYPE6v2 i386-absoft
  LIBS_LIB6v2        $(MODINSTALL6v2)/modlib/libs.lib
  KEY_MODELLER6v2    xxxmodelxxxiranxxxjexxx (for instance)

The Modeller_model function is best used for turning a wurst
alignment into an initial set of modeller input files, and
making a quick model.
Errors are captured from the modeltopscript.log file, should
you wish to inspect them. If a model file really is produced, then
its filename is returned by the function. 

See the modellermodel.pl script for a working use of the function.
This script may be *all* that you need.


=head1 MODELLER HELP

The best place to go is always the documentation :

  http://salilab.org/modeller/manual/manual.html
  file:///data/procter/modeller6v2/doc/manual/manual.html

There are also example scripts for nearly all of the commands :
/data/procter/modeller6v2/examples/commands/

Finally, there are some hints within the mdl_yourmodel_template.top 
file which is buried in the .tgz file that contains all the info 
from the modeller run. Make a new directory in the place where
the script ran from, change to it, and extract the files by 
using this command :

   tar -zxf ../mdl_yourmodel_template.tgz


=head2 HOMO-OLIGOMERS

These can be modelled in modeller, but its a bit tricky (not as bad as
hetero-oligomers, but thats just another story). 

The basic script is described here :

   http://salilab.org/modeller/manual/node101.html#28724

Or locally,
   /data/procter/modeller6v2/examples/commands/define_symmetry.top
There is also a helpful email snippet below :

=over

Hi,

You can have more than one pair of segments defined trough the
DEFINE_SYMMETRY command, so you can restrict the three molecules to be
identical by defining, for example, three pairs. If you have three
monomers A, B, and C you would define pairs A-B, B-C, and C-A, that
should make all three identical.

The 'defsym' routine in the TOP file would look something like this
(note the values for the ADD_SYMMETRY option):

SUBROUTINE ROUTINE = 'defsym'

       SET RES_TYPES = 'ALL'
       SET ATOM_TYPES = 'MNCH'
       SET SELECTION_STATUS = 'INITIALIZE'
       SET SELECTION_SEARCH = 'SEGMENT'

       SET SYMMETRY_WEIGHT = 0.5
       PICK_ATOMS PICK_ATOMS_SET = 2, SELECTION_SEGMENT = '1:' '102:'
       PICK_ATOMS PICK_ATOMS_SET = 3, SELECTION_SEGMENT = '103:' '204:'
       DEFINE_SYMMETRY ADD_SYMMETRY = on off
       PICK_ATOMS PICK_ATOMS_SET = 2, SELECTION_SEGMENT = '103:' '204:'
       PICK_ATOMS PICK_ATOMS_SET = 3, SELECTION_SEGMENT = '205:' '306:'
       DEFINE_SYMMETRY ADD_SYMMETRY = on off
       PICK_ATOMS PICK_ATOMS_SET = 2, SELECTION_SEGMENT = '205:' '306:'
       PICK_ATOMS PICK_ATOMS_SET = 3, SELECTION_SEGMENT = '1:' '102:'
       DEFINE_SYMMETRY ADD_SYMMETRY = on off
       
       RETURN

END_SUBROUTINE

=back

  
=head2 EXPORT



=head1 SEE ALSO

    Wurst
    Modeller web-page at www.salilab.org

=head1 AUTHOR
    
    Jim
    
=head1 COPYRIGHT AND LICENSE

Copyleft 2004 by James Procter

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself. 

=cut
