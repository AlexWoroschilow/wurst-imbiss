#!/usr/bin/perl

# calcpend - single process calculation of pending structural alignment pairs

use strict;
use FindBin;
use lib "$FindBin::Bin/.";
use lib "$FindBin::Bin/../../src/Wurst/blib/arch";
use lib "$FindBin::Bin/../../src/Wurst/blib/lib";

use POSIX qw(EXIT_SUCCESS EXIT_FAILURE);

use Wurst;

use Wurst::SeqStrCmp;

use Wurst::AlDb;

DBI::trace(1); # Debugging

use Wurst::AlDb::StrAlgn;

use Wurst::AlDb::Alignment;
use Wurst::AlDbUtils;
use Parallel::MPI::SimD;
use Storable qw(nfreeze thaw);
use Compress::Zlib qw(uncompress compress);

my $verbosity = 1; # debug!
my $start_time = time(); # baseline.

# helper function again.


sub prof_bin_valid {
    my ($chainid, $binpath, $profpath) = @_;
    my ($coord, $prof);
    (-f $binpath."$chainid.bin") and ($coord = coord_read($binpath."$chainid.bin")) 
        and (-f $profpath."$chainid.prof") and ($prof = blst_chk_read($profpath."$chainid.prof"))
        and return ($coord,$prof);
    return undef;
};

sub Make_coordprof_lib( $ $; $ ) {
    my ($pdbpath, $profpath, $Protiterator) = @_; 
    my $libhash={};
    # Protiterator could specify a subset of proteins to make in the library.
    (not defined($Protiterator)) and $Protiterator = Wurst::AlDb::Protein->retrieve_all();
    my $prot;
    while ($prot = $Protiterator->next()) {
        my ($pdbid, $coord, $prof) = ($prot->pdbid, prof_bin_valid($prot->pdbid, $pdbpath, $profpath));
        if (defined($prof)) {
            $$libhash{$pdbid}=compress(nfreeze( [coord_pack($coord), seqprof_pack($prof)]));
        } else {
            warn "Serious: skipping $pdbid because there's no valid profiles or bin's for it\n";
            $$libhash{$pdbid}=compress(nfreeze([-1])); # just to say its part of the library
        }
    }
    return $libhash;
}

# note - this function writes back to destroy the packed library scalar
sub Inflated_coordprof_lib( $ ) {
#    my ($lhash) = @_; 
    my ($libhash) = @_;
#    my $libhash=thaw(uncompress(${$libmesg}));
    (defined($libhash)) 
        or die("Couldn't unpack library!");
#    ${$libmesg}=\0; # memory efficiency.
    my $prot;
    foreach $prot (keys %{$libhash}) {
        my $msg = $$libhash{$prot};
        $$libhash{$prot} = thaw(uncompress($msg));
        if (scalar @{$$libhash{$prot}} == 2) {
            $$libhash{$prot}->[1] = seqprof_unpack($$libhash{$prot}->[1]);
            $$libhash{$prot}->[0] = coord_unpack($$libhash{$prot}->[0]);
            if ($verbosity>21) { print "".(time()-$start_time)."Unpacked $prot\n" };
        } else {
            if ($verbosity>11) { print "".(time()-$start_time)."Skipping $prot in library\n"; };
        }
    }
    ($verbosity>11) and print "".(time()-$start_time)." ".`uname -n`."unpacked library.\n";
    return $libhash;
}




# globals
my ($paramHash, $AlType); # on a slave - this is a hash of values
my $Protein_lib; # hash of all protein coords and profiles
my ($sA,$sB,$fname1,$fname2); 

#Set_Verbosity(100); 
Set_Old_Format_Op(4.5);

my $Pdb_path = "/bm/wurst_server/FoldLibs/EvaNR12feb2004/"; # these will be overwritten
my $Profilespath="/bm/wurst_server/FoldLibs/EvaNR12feb2004_prof/";
my $max_alignments = 10000; # rough guess

my $Usage = "$0 [alignment-name] [max_alignments to make]\n";


# SimD functions

sub init_parse_args ( $ ) {
    my $a = shift @_;
    print "".(time()-$start_time)."$0\n";
    $SIG{TRAP} = 'IGNORE'; kill 'TRAP', $$;

    if (scalar @{$a}>2) {
        warn "Up to two arguments!\nUsage :\n$Usage\n";
        return(EXIT_FAILURE,undef);
    }
    
    if ((scalar @{$a}) and ($$a[0] =~ m/[^-0-9]+/)) {
        # looks like an alignment name
        my $al_name = shift @{$a};
        my $altype = Wurst::AlDb::AlType->search(name=>$al_name);
        if (defined($altype)) {
            $AlType = $altype;
        } else {
            warn "Unknown alignment type $al_name\nSorry. $Usage\n";
            return(EXIT_FAILURE,undef);
        }
    }
    if (scalar @{$a}) {
        $$a[0]=~s/[^-0-9]//g;
        $max_alignments = $$a[0];
    }
    if (not defined($AlType)) {
        # set altype from first pending entry.
        my $pend = Wurst::AlDb::PendingStrAl->retrieve_from_sql("1 limit 1");
        if (not $pend->count()) {
            warn("No pending pairs - nothing to be done.");
            return(EXIT_FAILURE, undef);
        }
        $AlType = $pend->next()->altype(); # take the most common one in the set (costs several minutes!)
    }
    ($verbosity>0) and print "".(time()-$start_time)." Initial AlType is ".$AlType->name()."\n";
    # construct init string - will be the needed protein profiles and coordinates
    $paramHash = Storable::thaw($AlType->params);
    $Profilespath = $$paramHash{profiles_path};
    $Pdb_path = $AlType->pdbpath;
    ($verbosity>10) and print "".(time()-$start_time)."making library\n";
    # $Protein_lib = Make_coordprof_lib($Pdb_path, $Profilespath);
    $Protein_lib = undef;
    ($verbosity>11) and print "".(time()-$start_time)."made.\n";
    $SIG{TRAP} = 'IGNORE'; kill 'TRAP', $$;
    return (EXIT_SUCCESS, [ $AlType->name, $Pdb_path, $paramHash, $Protein_lib]);
};


sub set_defaults ( $ ) {
    my $init_string = shift @_;
    $AlType = $$init_string[0];
    $Pdb_path = $$init_string[1];
    $paramHash = $$init_string[2]; # keep our own copy.
    if (defined($$init_string[3])) {
        $Protein_lib = Inflated_coordprof_lib($$init_string[3]);
    } else {
        $Protein_lib = {};
    };
    $Profilespath = $$paramHash{profiles_path};
    if (not(defined($AlType) and defined($paramHash) and ($Pdb_path or keys %{$Protein_lib}))) {
        warn(`uname -a`."\nDidn't receive correct init.\n");
        return (EXIT_FAILURE);
    }
    $init_string=\0; # free memory!
    x_sq_set_alignparams($paramHash);
    x_sq_prot_alignparams(1); # prevent x_prof updating with hardwired parameters.
    
    return(EXIT_SUCCESS);
}

sub make_pairs {
    # we already know about the parameters in $AlType;
    my $workset=[];
    ($verbosity>1) and print "".(time()-$start_time)." Counting pending alignments.\n";
    my $numpend = Wurst::AlDb::PendingStrAl->count_all();
    ($max_alignments < 0) and $max_alignments = $numpend;
    ($verbosity>1) and print "".(time()-$start_time)." Retrieving pending alignments.\n";
    my $pending_pair = Wurst::AlDb::PendingStrAl->retrieve_from_sql("altype=".$AlType." limit $max_alignments");
    
    if ($numpend) {
        print "".(time()-$start_time)." Will do ".(($numpend>$max_alignments) ? $max_alignments : $numpend)." alignments out of $numpend\n";
    } else {
        print "".(time()-$start_time)." Nothing to be done.\n";
    }
    my $al_pend;
    while (($al_pend = $pending_pair->next())and ($numpend--)) {
        my ($pdbA) = $al_pend->protA->pdbid;
        my ($pdbB) = $al_pend->protB->pdbid;
        if (($al_pend->protA->region!=0) || ($al_pend->protB->region!=0)) {
            die "this version doesn't deal with complex region selections!\n";
        }
        push @$workset, [$al_pend->Pendingstral_id, $pdbA, $pdbB];
    }
    
    return ($workset);
}

sub work_function {
    my $packet = shift @_;
    (defined($packet) and (scalar @{$packet}==3)) or return(BROKE_WORK);
    ($verbosity>11) and print "".(time()-$start_time)." ".`uname -n`." is working on ".$$packet[1]." vs ".$$packet[2]."\n";

    # see if we know about this protein
    my ($protA,$protB) = ($$packet[1], $$packet[2]);
    if (not ((exists $$Protein_lib{$protA}) and (defined($$Protein_lib{$protA})) and (scalar @{$$Protein_lib{$protA}}>1))) {
        my ($coord, $prof) = prof_bin_valid($protA, $Pdb_path, $Profilespath);
        ($coord and $prof) or(($verbosity>3) and print "".(time()-$start_time)." ".`uname -n`." can't get pdb $protA\n");
        ($coord and $prof) or return(BROKE_WORK);
        $$Protein_lib{$protA}=[$coord, $prof];
    }
    if (not ((exists $$Protein_lib{$protB}) and (defined($$Protein_lib{$protB})) and (scalar @{$$Protein_lib{$protB}}>1))) {
        my ($coord, $prof) = prof_bin_valid($protB, $Pdb_path, $Profilespath);
        ($coord and $prof) or (($verbosity>3) and print "".(time()-$start_time)." ".`uname -n`." can't get pdb $protB\n");
        ($coord and $prof) or return(BROKE_WORK);
        $$Protein_lib{$protB}=[$coord, $prof];
    }
    
    my @res_vec = x_prof_str_align($$Protein_lib{$$packet[1]}->[1],$$Protein_lib{$$packet[1]}->[0],$$Protein_lib{$$packet[2]}->[1],$$Protein_lib{$$packet[2]}->[0]);
    (scalar @res_vec < 10) and return (BROKE_WORK);
    return (HAPPY_SLAVE, [@res_vec]);
};

# this function doesn't need to cope with duplicates.
sub accept_results {
    my ($wpack, $result) = @_;
    my $al_pend = Wurst::AlDb::PendingStrAl->retrieve($$wpack[0]);
    if ($al_pend) {
        my $stralgn = Wurst::AlDb::add_StrAlgn( $al_pend->altype,  $al_pend->protA, $al_pend->protB, [ x_sq_struct_legend() ], $result);
        if (defined $stralgn) {
            $al_pend->delete();
        } else {
            warn((time()-$start_time)." Warning !  Failed to add alignment for ".$al_pend->protA->pdbid." against ".$al_pend->protB->pdbid."\n");
        }
    } else {
        ($verbosity>0) and print "".(time()-$start_time)." Alignment entry collision for ".$al_pend->protA->pdbid." against ".$al_pend->protB->pdbid."\n";
    };
}

sub tidy_up {
    print "".(time()-$start_time).":".(Wurst::AlDb::PendingStrAl->count_all())." Remaining pending alignments.";
}
#SimD_Debug(1);

$SIG{TRAP} = 'IGNORE'; kill 'TRAP', $$;

exit(SimD_main( \&init_parse_args, \&set_defaults, \&make_pairs, \&work_function, \&accept_results, \&tidy_up, @ARGV));
