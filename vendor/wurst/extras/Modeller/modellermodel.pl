#!/usr/bin/perl
#$ -S /usr/bin/perl
#$ -cwd
#$ -j y
#$ 
# uses latest sequence structure parameters to make an alignment
# and generate a model using modeller

use strict;
use FindBin;
my $wurst_installpath = "$FindBin::Bin/../../";
#__wurstli__ The wurst library location


use lib "$FindBin::Bin/../../src/Wurst/blib/lib";
use lib "$FindBin::Bin/../../src/Wurst/blib/arch";
#local paths
use lib "$FindBin::Bin/Wurst/Modeller/blib/lib";
use lib "$FindBin::Bin/Wurst/Modeller/blib/arch";
use lib "$FindBin::Bin/../Guts/src/Wurst/Guts/blib/lib";
use lib "$FindBin::Bin/../Guts/src/Wurst/Guts/blib/arch";

use lib "$FindBin::Bin/../AlDb/"; # for SeqStrCmp.pm

use Wurst;
use Wurst::Guts qw(seqprof_trim pair_set_shift);
use Wurst::SeqStrCmp;
use Wurst::Modeller;

my $verb = 0;

my @bindir = ("/bm/wurst_server/FoldLibs/pdb90/","/bm/wurst_server/foldlibs/EvaNR12feb2004/", "/bm/pdb_bin_pool");
my @profdir = ("/bm/wurst_server/FoldLibs/pdb90_prof/","/bm/wurst_server/foldlibs/EvaNR12feb2004_prof/"); # , "/bm/pdb_bin_prof");



# From the wurst_server script 10th June 2004
my %prof_params = (
                   "align_type" => 'prof_fx_j10',
                   "sw1_pgap_open" => 1.140,
                   "sw1_pgap_widen" => 0.190,
                   "sw1_qgap_open"  =>  2.733,    
                   "sw1_qgap_widen" =>  0.952,
                   "sec_s_scale"    =>   0.0,
                   "sw1_sec_s_scale" =>0.0,
                   "sw2_sec_s_scale" => 0.0,
                   "nw_sec_s_scale" => 0.0,
                   "nw_sec_pnlty" => 0.,
                   "fx_mat_shift"       =>   1.987,
                   "seq_wt"         =>   0.830,
                   "zero_shift"     =>   0.163,
                   
                   "nw_pgap_open"  => 1.140,
                   "nw_qgap_open"  =>    2.733,
                   "nw_pgap_widen" =>    0.190,
                   "nw_qgap_widen" =>    0.952,
                   "nw_rs_scr2" => .13,
                   "sw_rs_scr2" => .18,
                   "sw_gap_geo" => -1.12,
                   "sw_seq_gap" => -0.57,
                   "sw_str_gap" => -2.18,
                   "sw_str_wdn" => -0.21,
                   "nw_gap_geo" => -1.27,
                   "nw_seq_gap" => -0.09,
                   "nw_str_gap" => -2.35,
                   "nw_str_wdn" => -0.93,
#                   "sim_mat_bottom" => 
                   "MATR_NAME" => $wurst_installpath."/matrix/725tmp.out",
                   "rescore_paramsfile" => $wurst_installpath."/params/cc_allat_p891+0.param",
                   
                   "PARAM_FILE" => $wurst_installpath."/params/g9.param"
                   );

# ----------------------------------------
# Check arguments
my $Usage = "Usage : $0 [-v+plgxchmX,YnX,Y] Modelname BlastProfileFile Template_pdbid[:seqX-seqY] [Template_pdbid[:seqX-seqY]]+ [ dir1 .. ] 
seqX-seqY specifies range of sequence to align against a particular template
Currently looking for pdbid.bin files in \n".(join ",",@bindir)."\n and corresponding profiles in \n".(join ",",@profdir)."\n";
my $Opts = (join "\n",
            "-v+ more v's for more verbosity",
            "-p  prepare files - do not actually run anything",
            "-l  turns on loop modelling",
            "-g  sets all refine levels to refine_3",
            "-mX,Y sets MODEL_START,MODEL_END",
            "-nX,Y sets LOOP_MODEL_START, LOOP_MODEL_END",
            "-x  model_level to refine_3 and loop_level to refine_5 (xtremely long)",
            "-c  careful - sets refine_hot_only on",
            "-t  forces the use of /tmp \n    for when you get errors like\n    'can't create file in /local : permission denied'\n",
            "-h  this message.");
sub prof_bin_valid {
    my ($chainid) = @_;
    my ($coord, $prof);
    my ($di, @d);
    @d = @bindir;
    $coord = undef;
    do {
        $di=shift @d;
        (-f $di."$chainid.bin") and ($coord = coord_read($di."$chainid.bin"));
    } while ((not defined($coord)) and scalar @d);
    @d = @profdir;
    $prof = undef;
    do {
        $di=shift @d;
        (-f $di."$chainid.prof") and ($prof = blst_chk_read($di."$chainid.prof"));
    } while ((not defined($prof)) and scalar @d);
    (defined($coord) 
     and defined($prof)) and return ($coord, $prof);
    defined($coord) or 
	warn "Can't find $chainid in ".(join ",",@bindir)."\n";
    defined($prof) or 
	warn "Can't find $chainid in ".(join ",",@profdir)."\n";
    return undef;
};

my ($modelname, $s_profn, $s_prof, $pdbid);
my ($m_low, $m_high, $lm_low, $lm_high, 
    $loop_level, $loopm, $ref_level, $hot_refine);
my $argp = shift @ARGV;
my $proc = 1;
# set up the modeller parameter hash from arguments
my %m_parms;

while ($proc and defined($argp)) {
    $proc = 0;
    if ($argp=~/^-/) {
	if ($argp=~/[Vv]+/) {
	    ($verb) = ($argp=~/([Vv]+)/);
            $verb = length($verb);
	    $proc = 1;
	}
	if ($argp=~/[pP]/) {
	    $m_parms{NO_SHELL_RUN} = 1;
	    $proc = 1;
	}
	if ($argp=~/[lL]/) {
	    $m_parms{DO_LOOPS} = 1;
	    $proc = 1;
	}
	if ($argp=~/[gG]/) {
	    $m_parms{MD_LEVEL} = $m_parms{LOOP_MD_LEVEL} = 'refine_3';
	    $proc = 1;
	}
	if ($argp=~/[xX]/) {
	    $m_parms{MD_LEVEL} = 'refine_3'; # why bother here.
            $m_parms{LOOP_MD_LEVEL} = 'refine_5';
	    $proc = 1;
	}
	if ($argp=~/[cC]/) {
	    $m_parms{REFINE_HOT_ONLY} = 1;
	    $proc = 1;
	}
	if ($argp=~/[mM]/) {
	    my @rg = $argp=~/[mM]([0-9]+),([0-9]+)/;
	    if (scalar @rg==2) {
		if (($rg[0]>0)
		    and ($rg[0]<=$rg[1])) {
		    ($m_parms{STARTING_MODEL},
                     $m_parms{ENDING_MODEL}) = @rg;
                    $proc=1;
		} else {
		    die "-mX,Y where X<Y.\n$Usage";
		}
	    } else {
		die "Set model range like -mX,Y where integers X<=Y\n$Usage";
	    }
	}
	if ($argp=~/[nN]/) {
	    my @rg = $argp=~/[nN]([0-9]+),([0-9]+)/;
	    if (scalar @rg==2) {
		if (($rg[0]>0)
		    and ($rg[0]<=$rg[1])) {
		    ($m_parms{LOOP_MODEL_START},
                     $m_parms{LOOP_MODEL_END}) = @rg;
                    $proc=1;
		} else {
		    die "Set loop model range like -nX,Y where X<=Y.\n$Usage";
		}
	    } else {
		die "Set loop model range like -nX,Y where integers X<=Y\n$Usage";
	    }
	}
        if ($argp=~/[tT]/) {
            $m_parms{localdisk}="/tmp";
            $proc=1;
        }
	if ($argp=~/[hH]/) {
	    print $Usage."\n$Opts\n";
	    exit(0);
	}
	($proc == 1) and ($argp = shift @ARGV); # throw away args even if there was junk
    } 
}
if ($verb>1) {
    print "-- Modelling Parameters :\n";
    foreach my $q (keys %m_parms) { print "--  $q = ".$m_parms{$q}."\n"; }
}
my (@tplts, @srng, $s_len);
$modelname = $argp;
((defined($modelname))
 and (defined($s_profn = shift @ARGV))
 and defined($pdbid = shift @ARGV)) or die "$Usage";
$s_prof = blst_chk_read($s_profn) or die "Failed to read $s_profn as a blast profile!\n";
$s_len = seq_size(seqprof_get_seq($s_prof));

my ($pr,$rng) = $pdbid=~/([1-9][A-Za-z0-9]{3}[A-Z_])(:[-0-9,]+)?/;
if ($rng) {
    push @srng, [ $rng=~/(\d+)[,-](\d+)/];
    push @tplts, $pr;
} else {
    push @tplts, $pdbid;
    push @srng, [1, $s_len];
}    

# anything that remains are new searchpaths
# if they look like a directory ???

while (scalar @ARGV) {
    my $newd = shift @ARGV;
    if (-d $newd) {
        unshift @bindir, $newd;
        unshift @profdir, $newd;
    } else {
        if ($newd=~/[1-9][A-Za-z0-9]{3}[A-Z_]/) {
            ($pr,$rng) = $newd=~/([1-9][A-Za-z0-9]{3}[A-Z_])(:[-0-9,]+)?/;
            if (defined($rng)) {
                push @srng, [ $rng=~/(\d+)[-,](\d+)/];
                push @tplts, $pr;
            } else {
                push @tplts, $newd;
                push @srng, [1, $s_len];
            }    
        } else {
            print "# Ignoring $newd\n";
        }
    }
}
my (@psets,@coords, @cnam);
my ($thr_prof, $ps_offset);

$pr = 0;
do {
    $pdbid= $tplts[$pr];
    if (($srng[$pr][0] == 1)
        and ($srng[$pr][1] == $s_len)) {
        $thr_prof = $s_prof;
        $ps_offset = 0;
    } else {
        if ($thr_prof = seqprof_trim($s_prof, 
                                     $srng[$pr][0],
                                     $srng[$pr][1])) {
            $ps_offset = $srng[$pr][0]-1;
        } else {
            $thr_prof = $s_prof;
            warn("ignoring invalid range for $pdbid : [".(join ",",@{$srng[$pr]})."]\n");
        }
    }
        
    my ($coord, $coord_prof) = prof_bin_valid($pdbid);
    (defined($coord)) or die "$Usage";

    x_sq_set_alignparams ( \%prof_params );
    
    my $fx_param = param_fx_read ($prof_params{PARAM_FILE}) || die "Failed to read ".$prof_params{PARAM_FILE}."\n";
    my $submat = sub_mat_read ($prof_params{MATR_NAME}) || die "Fail on ".$prof_params{MATR_NAME}."\n";
    zero_shift_mat ($submat, $prof_params{zero_shift});
    
    my $pfile;
    if (defined($prof_params{rescore_paramsfile})) {
        $pfile = $XSQ_rescore_paramsfile;
    } else {
        $pfile = "$FindBin::Bin/../../params";
        $pfile = $pfile.'/cc_allat_p891+0.param';
    };
    
    my $rescore_params = param_rs_read ($pfile) || die "Failed to read Rescore params $pfile";
    
# ----------------------------------------
    
    my ($total_sw_score,
        $full_pairset,
        $tot_nw_score);
    
    my $aldata = do_prof_align( $thr_prof, $coord, $coord_prof, 'none', $fx_param, $submat, $rescore_params);
    $total_sw_score = shift @$aldata;
    $full_pairset = shift @$aldata;
    if ($ps_offset) {
        if ($verb>10) {
            print "unshifted alignment\n".pair_set_pretty_string($full_pairset,seqprof_get_seq($thr_prof), coord_get_seq($coord))."\n";
        }

        pair_set_shift($full_pairset, $ps_offset); # adjust for sequence range
    }
    $tot_nw_score = shift @$aldata;
    
    my $safename=$modelname;
    
    my $seq = seqprof_get_seq($s_prof);
    my $model = make_model($full_pairset, $seq, $coord);
    $safename=~s/[^A-Za-z0-9_]+//g;
    if ($safename != $modelname) { print "# modeller models are called $safename\n"; };
    my ($cov,$nothing) = pair_set_coverage($$aldata[6], seq_size($seq), coord_size($coord));
    if ($verb > 2) { 
        print "S & W coverage\n";
        $cov =~ s/1/X/g;
        $cov =~ s/0/\-/g;
        print $cov, "\n";
    };
    if ($verb) { 
        printf "Score is %.4g sw cover: %.2f nw cover %.2f\n",
        $total_sw_score, $$aldata[0], $$aldata[1];
        print "\n".pair_set_pretty_string($full_pairset,$seq, coord_get_seq($coord))."\n";
    }
    coord_2_spdb("$safename.$pdbid.pdb", $model, $$aldata[(scalar @{$aldata})-3],$seq);
    push @cnam, $pdbid;
    push @coords, $coord;
    push @psets, $full_pairset;
    
} while (++$pr<scalar @tplts);
if (($pr>1) && ($verb)) {
    # print the msa equivalent nicely
    my $lw = 72;
    my $ll = 1;
    my $f = 0;
    my @msa = Wurst::Modeller::pir_mulpairs( seqprof_get_seq($s_prof), \@coords, \@psets);
    my @ls;
    while (!$f) {
        for ($pr=0;$pr<scalar @cnam; $pr++) { 
            if (length $msa[$pr]<=($ll*$lw)) {
                $f=1;
            }
            push @ls, $cnam[$pr]."\t".(substr $msa[$pr], 72*($ll-1), 72)."\n";
        }
        push @ls, "$modelname\t".(substr $msa[$pr], 72*($ll-1), 72)."\n\n";
        $ll++;
    }
    print @ls;
}

# now make the models
my @gmod = Modeller_m_model( $modelname, seqprof_get_seq($s_prof), 
                             "Made from a bunch of template alignments",
                             \@cnam, \@coords, \@psets, \%m_parms);
print "# made these models\n".(join "\n",@gmod,"");
