#!/usr/bin/perl
# aldb_addset.pl
# 
# derived from the test_addpend.pl script
# this sets up all the pend records for a new set of proteins
# these are :
# all records for each new protein as chainB in a comparison with all existing proteins as chainA
#  and all records for the current new protein as chainA compared to all remaining new proteins as chain B

# It only makes sense for this to be a single process script

use strict;
use FindBin;
use lib "$FindBin::Bin/.";
use lib "$FindBin::Bin/../../src/Wurst/blib/arch";
use lib "$FindBin::Bin/../../src/Wurst/blib/lib";

use POSIX qw(EXIT_SUCCESS EXIT_FAILURE);

use Wurst;
use Wurst::SeqStrCmp; # needed ?

use Wurst::AlDb;
use Wurst::AlDbUtils; # for what ?
use Wurst::AlDb::StrAlgn;
use Wurst::AlDb::Alignment; # needed ?

## helper functions (for speed here)
sub prof_bin_valid {
    my ($chainid, $binpath, $profpath) = @_;
    my ($coord, $prof);
    (-f $binpath."$chainid.bin") and ($coord = coord_read($binpath."$chainid.bin")) 
        and (-f $profpath."$chainid.prof") and ($prof = blst_chk_read($profpath."$chainid.prof"))
        and return ($coord,$prof);
    return undef;
};


$SIG{TRAP} = 'IGNORE'; kill 'TRAP', $$;

# set up the database connection data here
# Wurst::AlDb->database ("dbi:mysql:database=wurstreal;host=localhost;port=3380;mysql_socket=/tmp/procter-mysql","root","jim", {});
#Wurst::AlDb->database("dbi:mysql:database=wursttest;host=localhost;port=3380;mysql_socket=/tmp/procter-mysql",
#                      "root",
#                      "jim", {AutoCommit=>0});
Wurst::AlDb->database("dbi:mysql:database=wurstreal;host=localhost;port=3380;mysql_socket=/local/procter-mysql",
                      "root",
                      "jim", {AutoCommit=>0});


# Default alignment info. If none exists

my $Pdb_path = "/bm/wurst_server/FoldLibs/pdb90/";
my $Profilespath="/bm/wurst_server/FoldLibs/pdb90_prof/";
my $altype = Wurst::AlDb::AlType->find_or_create(name=>"swprofile_fx");

if (!defined($altype->pdbpath)) { 
    # take current setting
    use Wurst::SeqStrCmp;
    Set_Verbosity(100);
    Set_Old_Format_Op(4.5);
    set_swprof_params();
    my $params = x_sq_get_alignparams();
    $altype->name($$params{align_type}) if (not $$params{align_type} eq $altype->name);
    $$params{profiles_path} = $Profilespath;
    $altype->pdbpath($Pdb_path);
    $altype->params(Storable::nfreeze($params));
    $altype->update;
} else {
    print "Setting alignment parameters from database entry\n";
    my $params = Storable::thaw($altype->params);
    x_sq_set_alignparams($params);
    $Profilespath = $$params{profiles_path};
    $Pdb_path = $altype->pdbpath;
}

my (%p_id, @existing_prots);
my (%np_id, @new_proteins,@np_com); # the to be accepted proteins
my @P_bits; # bits of a valid chain
my ($p_i, $i);

$p_i = Wurst::AlDb::Protein->retrieve_all();
while ($i=$p_i->next()) {
    $p_id{$i->pdbid} = $i->protein_id;
}

my $Usage = "Usage : $0 .. [-f listofpdbsinafile] pdbid ..... \n
 Scans list, finding new chains to add to the database.
 For each new chain, make PendingStrAlgn entries to all existing\n";

while (scalar @ARGV) {
    my $pdbid = shift @ARGV;

    if ((lc $pdbid) eq "-f") {
        my $pdblist = shift @ARGV;
        
        if (-f $pdblist) {
            my $p;
            (open PDBLIST, "$pdblist") or die "Can't open $pdblist\n$Usage";
            do {
                $p=<PDBLIST>;
                if (defined($p)) {
                    substr($p,0,4) = lc substr($p,0,4);
                    my ($chain) = $p=~/([1-9][a-z0-9]{3}[A-Z_0-9])/;
                    if (not (exists($p_id{$chain}) and defined($p_id{$chain}))) {
                        @P_bits = prof_bin_valid($chain, $Pdb_path, $Profilespath);
                        if (scalar @P_bits>1) {
                            my $np = Wurst::AlDb::Protein->create({pdbid=>$chain, region=>0});
                            push @np_com, $np;
                            push @new_proteins, $chain;
                            $np_id{$chain} = $np->protein_id;
                        }
                    }
                }
            } while ($p and not eof);
            close PDBLIST;
        }
    } else {
        substr($pdbid,0,4) = lc substr($pdbid,0,4);
        my ($chain) = $pdbid=~/([1-9][a-z0-9]{3}[A-Z_0-9])/;
        if (not (exists($p_id{$chain}) and defined($p_id{$chain}))) {
            @P_bits = prof_bin_valid($chain, $Pdb_path, $Profilespath);
            if (scalar @P_bits>1) {
                my $np = Wurst::AlDb::Protein->create({pdbid=>$chain, region=>0});
                push @np_com, $np;
                push @new_proteins, $chain;
                $np_id{$chain} = $np->protein_id;
            }
        }
    }
}

map { $_->dbi_commit(); } @np_com;
(scalar @new_proteins == 0) and die "No new proteins !\n$Usage";


# Make new Pending relations
# pattern this is done has been modified.

my ($newP, $otherP, @plist);
my $pendingr = 0;

if (not scalar keys %p_id) {
    # start the ball rolling...
    my $p = shift @new_proteins;
    $p_id{$p} = $np_id{$p};
}

my @pending_clist;

while (scalar @new_proteins) {
    $newP = shift @new_proteins;
    @plist = keys %p_id;


    while ($otherP = shift @plist) {
        if ((defined($p_id{$otherP})) and ($p_id{$otherP}<$np_id{$newP})) {
            my $pend = Wurst::AlDb::PendingStrAl->create({protA=>$p_id{$otherP}, protB=>$np_id{$newP}, altype => $altype});
            if (defined($pend)) {
                push @pending_clist, $pend;
                ++$pendingr;
                
            } else {
                warn ("Failed to create pending AlignStr (".$otherP.",".$newP.")\n");
            }
        }
    }
    if (scalar @pending_clist>20000) {
        Wurst::AlDb::->dbi_commit();
        # don't need to commit each one.
        @pending_clist = ();
    }
    $p_id{$newP} = $np_id{$newP};
};
(scalar @pending_clist) and (Wurst::AlDb::->dbi_commit());
@pending_clist=();

print "Added $pendingr new pending relations. There are now ".(Wurst::AlDb::PendingStrAl->count_all())."\n";



