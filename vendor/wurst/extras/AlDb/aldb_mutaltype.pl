#!/usr/bin/perl
# aldb_mutaltype
# reset/set the alignment table with a new alignment style
# this modifies an existing entry to new parameters (the name changes too!)
use strict;
use FindBin;
use lib "$FindBin::Bin/.";
use lib "$FindBin::Bin/../../src/Wurst/blib/arch";
use lib "$FindBin::Bin/../../src/Wurst/blib/lib";

use POSIX qw(EXIT_SUCCESS EXIT_FAILURE);

use Wurst;
#use lib "$FindBin::Bin//Wurst";
use Wurst::SeqStrCmp;

use Wurst::AlDb;
use Wurst::AlDb::StrAlgn;
#use Wurst::Protein;
use Wurst::AlDb::Alignment;


$SIG{TRAP} = 'IGNORE'; kill 'TRAP', $$;
my ($sA,$sB,$fname1,$fname2); 

Set_Verbosity(100);
Set_Old_Format_Op(4.5);
#my $Pdb_path = "/bm/wurst_server/FoldLibs/EvaNR12feb2004/";
#my $Profilespath="/bm/wurst_server/FoldLibs/EvaNR12feb2004_prof/";
my $Pdb_path = "/bm/wurst_server/FoldLibs/pdb90/";
my $Profilespath="/bm/wurst_server/FoldLibs/pdb90_prof/";

use Wurst::AlDb;

my $def_name = "swprofile_fx"; # this is the altype entry that is to be modified. - 

# here we should get or create the AlType entry.
# for now we just want to update the param hash with all parameters
# typically this might be for adding additional pdb paths, or something.

my $altype = Wurst::AlDb::AlType->retrieve(name=>$def_name);
die ("Can't find $def_name\n") unless defined($altype);
if (defined($altype->params())) {
    my $yesorno=shift @ARGV;
    (defined($yesorno)) or $yesorno="no";
    if ($yesorno eq "overwrite.") {
        
        if (Wurst::AlDb::StrAlgn->count_ofaltype($altype)>0) {
            
            $yesorno = shift @ARGV;
            if ($yesorno eq "remove.") {
                Wurst::AlDb::StrAlgn->search(altype=>$altype)->delete_all();
              } else {
                  die "StrAlignments exist with the old parameter set!\nAdd 'remove.' to arguments if you want to delete these.\n";
              }
        }
    } else {
        die ("Need 'overwrite.' as an argument to overwrite the existing $def_name aligntype.\n");
    }
}

set_swprof_params();
#set_prof_alignparams();

my $params = x_sq_get_alignparams();

# any more changes to be made ?

$altype->name($$params{align_type}) if (not $$params{align_type} eq $altype->name);

$$params{profiles_path} = $Profilespath;
$altype->pdbpath($Pdb_path);
$altype->params(Storable::nfreeze($params));
$altype->update;
