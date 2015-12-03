#!/usr/bin/perl
# aldb_rollback_altype.pl
# if we've used the wrong parameters for an alignment,
# its better to stop as soon as possible, and recreate
# the PendingStrAl records.
# Also useful when we've had to change the info. stored
# for any particular alignment (doh)

# aldb_setaltype
# reset/set the alignment table with a new alignment style

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
my $Pdb_path = "/bm/wurst_server/FoldLibs/EvaNR12feb2004/";
my $Profilespath="/bm/wurst_server/FoldLibs/EvaNR12feb2004_prof/";

use Wurst::AlDb;

Wurst::AlDb->database("dbi:mysql:database=wursttest;host=localhost;port=3380;mysql_socket=/local/procter-mysql",
                      "root",
                      "jim", {AutoCommit=>0});


my $def_name = shift @ARGV;
(not defined($def_name)) and die "Must give an existing altype name to retrieve the invalid alignments.\n";
#"swprofile_fx"; # this is the altype entry that is to be modified.

# here we should get or create the AlType entry.
# for now we just want to update the param hash with all parameters
# typically this might be for adding additional pdb paths, or something.

my ($altype) = Wurst::AlDb::AlType->search(name=>$def_name);
(not defined($altype)) and die "$def_name isn't in the list of AlTypes.\n";
my $strals = Wurst::AlDb::StrAlgn->search(altype=>$altype);
my $i=0;
my $str;
while ($str=$strals->next()) {
    
    my $newpend = Wurst::AlDb::PendingStrAl->create({protA=>$str->protA(),
                                                     protB=>$str->protB(),
                                                     altype=>$altype}
                                                    );
    if (defined($newpend)) {
        $str->delete(); # should delete alignments too.
        $i++;
    };
}
print "Committing.\n";
Wurst::AlDb->dbi_commit();
