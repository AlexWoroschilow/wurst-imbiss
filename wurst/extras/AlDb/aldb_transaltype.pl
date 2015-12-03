#!/usr/bin/perl
# aldb_transaltype.pl
# Not much use for this except when fixing a broken database.
# update the entries in a set of StrAlignments with a new alignment style

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

use Wurst::AlDb;

my $def_name = "swprofile_fx"; # this is the altype entry that is to be modified.

# here we should get or create the AlType entry.
# for now we just want to update the param hash with all parameters
# typically this might be for adding additional pdb paths, or something.

my $altype = Wurst::AlDb::AlType->search(name=>$def_name);

my $strals = Wurst::AlDb::StrAlgn->retrieve_all();
my $i=0;
my $str;
while ($str=$strals->next()) {
    if (not defined($str->altype)) {
        $str->altype($altype);
        $str->update;
        $i++;
    }
}
