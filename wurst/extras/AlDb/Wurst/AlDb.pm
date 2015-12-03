#!/usr/bin/perl

package Wurst::AlDb;

use strict;
use warnings;
use base 'Class::DBI';

use vars qw( $VERSION );
$VERSION = 0.1;
#
# really simple structure

# Protein square = { proteinchain } - complete relational squares
# protein chain Set = { proteinchain ... }  - every protein chain referred to in DB
# proteinchain = { pdb_chid, length, sequence, pdb_bin_name, blast_profile_name }
# Alpair = { Aproteinchain Bproteinchain 

# following the cookbook.

# assertion. We must have some DB path already defined in the script using the class
#my $ALDB_DATABASE = "dbi:SQLite:/home/procter/src/prj/wurst/extras/AlDb/test.dbl";
my $verbosity = 0; # silent running

sub verbosity {
    my $self = shift;
    (scalar @_) and $verbosity = shift @_;
    return $verbosity;
}

my $dbh;

my $ALDB_DATABASE = "dbi:mysql:database=wursttest;host=localhost;mysql_socket=/local/procter-mysql";
my $ALDB_CREDENTIAL = "calc";
my $ALDB_CREDENTIAL2 = "";
my %ALDB_ATTRIB = Class::DBI->_default_attributes;

sub database ($ $ $; %) {
    my $self = shift @_;
    my ($dbn, $dbc1,$dbc2, $dbattr) = @_;
    my %aldb_attrib = %ALDB_ATTRIB;
    (defined($dbh)) and croak("Existing connection to $ALDB_DATABASE.\n");
    if (defined $dbattr) {
        foreach my $i (keys %$dbattr) {
            $aldb_attrib{$i} = $$dbattr{$i};
        }
    }
    eval { $dbh = DBI->connect($dbn, $dbc1, $dbc2, \%aldb_attrib); };
    if ($dbh) {
        $ALDB_DATABASE = $dbn;
        $ALDB_CREDENTIAL = $dbc1;
        $ALDB_CREDENTIAL2 = $dbc2;
        %ALDB_ATTRIB = %aldb_attrib;
    } else {
        $self->SUPER::_croak("* Script error
* Couldn't connect to :
* '$dbn'
* with $dbc1 user credential.
".(join "\n",$@, $DBI::errstr,"\n"));
    }
    return $dbh;
};

my $Verbosity = 0;

sub db_Main {
    if (not defined($dbh)) { 
        if ($Verbosity>0) {warn( "Opening default test database\n"); } ;
        $dbh = DBI->connect($ALDB_DATABASE, $ALDB_CREDENTIAL, $ALDB_CREDENTIAL2, \%ALDB_ATTRIB); };
    if ($Verbosity>0) {warn( "Opened database\n"); } ;
    return $dbh;
}

((not defined("$ALDB_DATABASE")) or ($@)) 
    and croak("Wurst::AlDB : Database environment globs not defined!\nSet $ALDB_DATABASE='dbi:sqlite:dbname'\n and $ALDB_CREDENTIAL='wurst' (user and password sed for DB)\n");

# packaged methods for adding to child classes and other essential operations




package Wurst::AlDb::Protein;
use base 'Wurst::AlDb';

Wurst::AlDb::Protein->table("Proteins");
Wurst::AlDb::Protein->columns(Primary=> qw/protein_id/);
Wurst::AlDb::Protein->columns(Essential=> qw/pdbid/);
Wurst::AlDb::Protein->columns(Occasional=> qw/region/);
Wurst::AlDb::Protein->columns(TEMP => qw/altype StrA_outdeg Pend_outdeg/);

Wurst::AlDb::Protein->set_sql(Pend_outdeg => qq{select PendingStrAls.altype AS altype, Proteins.protein_id, count(PendingStrAls.protB) as Pend_outdeg from Proteins,PendingStrAls where Proteins.protein_id=PendingStrAls.protA and PendingStrAls.altype=? group by  PendingStrAls.protA});

Wurst::AlDb::Protein->set_sql(StrA_outdeg => qq{select StrAlignments.altype AS altype,Proteins.protein_id, count(StrAlignments.protB) as StrA_outdeg from Proteins,StrAlignments where Proteins.protein_id=StrAlignments.protA and StrAlignments.altype=? group by  StrAlignments.protA});

#Wurst::AlDb::Protein->set_sql(Outdegs => qq{select Proteins.protein_id, count(PendingStrAls.protB) as Pend_outdeg, count(StrAlignments.protB) as StrA_outdeg from AlTypes,Proteins,PendingStrAls,StrAlignments where StrAlignments.altype=AlTypes.AlType_id and StrAlignments.protA=Proteins.protein_id and PendingStrAls.altype=AlTypes.AlType_id and PendingStrAls.protA=Proteins.protein_id group by  PendingStrAls.protA,StrAlignments.protA});



package Wurst::AlDb::AlType;
use base 'Wurst::AlDb';
Wurst::AlDb::AlType->table("AlTypes");
Wurst::AlDb::AlType->columns(Primary=> qw/AlType_id/);
Wurst::AlDb::AlType->columns(Essential=> qw/name/);
Wurst::AlDb::AlType->columns(Details=> qw/params pdbpath/);

package Wurst::AlDb::PendingStrAl;
use base 'Wurst::AlDb';
Wurst::AlDb::PendingStrAl->table("PendingStrAls");
Wurst::AlDb::PendingStrAl->columns(Primary=> qw/Pendingstral_id/);
Wurst::AlDb::PendingStrAl->columns(Essential=> qw/protA protB altype/);
Wurst::AlDb::PendingStrAl->has_a(altype=>'Wurst::AlDb::AlType');
Wurst::AlDb::PendingStrAl->has_a(protA=>'Wurst::AlDb::Protein');
Wurst::AlDb::PendingStrAl->has_a(protB=>'Wurst::AlDb::Protein');
# get the top 10 popular AlType entries in PendingStrAl
Wurst::AlDb::PendingStrAl->set_sql(common_altype => qq{SELECT altype, COUNT(altype) FROM PendingStrAls GROUP BY altype LIMIT 10});


package Wurst::AlDb;

# complex normalization search - operates on Protein, PendingStrAl and StrAlgn
# for each protein (i) and each aligntype (t)
#  there should be i-1 entries in StrAlgn where protB == i and altype == t
#   and be N-i-1 entries in StrAlgn where protA == i and altype == t
#   all of which should be distinct in the other protein entry (count duplicates)
#  deficits should be present as PendingStrAl entries. If they aren't there then they should be created.


# /* DBI entries */
my $SQL_table_creation = '
DROP TABLE Alignments/
DROP TABLE StrAlignments/
DROP TABLE Proteins/
DROP TABLE AlTypes/
DROP TABLE PendingStrAls/

CREATE TABLE Alignments (
        Alignment_id INTEGER PRIMARY KEY,
        seq     INTEGER, 
        str     INTEGER,
        altype  INTEGER,
        Total   FLOAT,  
        SwTot   FLOAT,  
        NwCov   FLOAT,  
        SwCov   FLOAT,  
        SwSc1   FLOAT,
        SwSc2   FLOAT,
        NwSc1   FLOAT,  
        NwSc2   FLOAT,
        T3      FLOAT, 
        T4      FLOAT, 
        T5      FLOAT, 
        SeqId   FLOAT, 
        Seqcov TEXT,  
        Strcov TEXT,
        Seqcovsw TEXT,
        Strcovsw TEXT,
        Localscore BLOB, 
        StrAlgn_id INTEGER 
)
/
CREATE TABLE StrAlignments (
        StrAlgn_id INTEGER PRIMARY KEY,
        protA   INTEGER, 
        protB   INTEGER, 
        altype  INTEGER, 
        A2B     INTEGER, 
        B2A     INTEGER, 
        ATotal   FLOAT,  
        BTotal   FLOAT,  
        ANwTot   FLOAT, 
        ASwTot   FLOAT,  
        ASwCov   FLOAT,  
        ASwSc1   FLOAT, 
        ASwSc2   FLOAT, 
        ASwSqG   FLOAT,  
        ASwStG   FLOAT,  
        BNwTot   FLOAT, 
        BSwTot   FLOAT,  
        BSwCov   FLOAT,  
        BSwSc1   FLOAT, 
        BSwSc2   FLOAT, 
        BSwSqG   FLOAT,  
        BSwStG   FLOAT,  
        T3      FLOAT, 
        T4      FLOAT, 
        T5      FLOAT, 
        SeqId   FLOAT, 
        Lal     FLOAT, 
        Lalpfrg FLOAT, 
        fragset TEXT, 
        covA    TEXT, 
        covB    TEXT,
        LscorA  BLOB, 
        LscorB  BLOB  
)
/
CREATE TABLE Proteins (
        protein_id INTEGER PRIMARY KEY,
        pdbid   CHAR(5), 
        region  INTEGER
)
/
CREATE TABLE AlTypes (
        AlType_id INTEGER PRIMARY KEY,
        name CHAR(24),
        params BLOB,
        pdbpath TEXT
)
/
CREATE TABLE PendingStrAls (
        Pendingstral_id INTEGER PRIMARY KEY,
        protA INTEGER,
        protB INTEGER,
        altype INTEGER
)
/
';
my $mysql_tables = '
DROP TABLE Alignments/
DROP TABLE StrAlignments/
DROP TABLE Proteins/
DROP TABLE AlTypes/
DROP TABLE PendingStrAls/

CREATE TABLE Alignments (
        Alignment_id INTEGER NOT NULL PRIMARY KEY AUTO_INCREMENT,
        seq     INTEGER, 
        str     INTEGER,
        altype  INTEGER,
        Total   FLOAT,  
        SwTot   FLOAT,  
        NwCov   FLOAT,  
        SwCov   FLOAT,  
        SwSc1   FLOAT,
        SwSc2   FLOAT,
        NwSc1   FLOAT,  
        NwSc2   FLOAT,
        T3      FLOAT, 
        T4      FLOAT, 
        T5      FLOAT, 
        SeqId   FLOAT, 
        Seqcov TEXT,  
        Strcov TEXT,
        Seqcovsw TEXT,  
        Strcovsw TEXT,
        Localscore BLOB, 
        StrAlgn_id INTEGER 
)
/
CREATE TABLE StrAlignments (
        StrAlgn_id INTEGER NOT NULL PRIMARY KEY AUTO_INCREMENT,
        protA   INTEGER, 
        protB   INTEGER, 
        altype  INTEGER, 
        A2B     INTEGER, 
        B2A     INTEGER, 
        ATotal   FLOAT,  
        BTotal   FLOAT,  
        ANwTot   FLOAT, 
        ASwTot   FLOAT,  
        ASwCov   FLOAT,  
        ASwSc1   FLOAT, 
        ASwSc2   FLOAT, 
        ASwSqG   FLOAT,  
        ASwStG   FLOAT,  
        BNwTot   FLOAT, 
        BSwTot   FLOAT,  
        BSwCov   FLOAT,  
        BSwSc1   FLOAT, 
        BSwSc2   FLOAT, 
        BSwSqG   FLOAT,  
        BSwStG   FLOAT,  
        T3      FLOAT, 
        T4      FLOAT, 
        T5      FLOAT, 
        SeqId   FLOAT, 
        Lal     FLOAT, 
        Lalpfrg FLOAT, 
        fragset TEXT, 
        covA    TEXT, 
        covB    TEXT,
        LscorA  BLOB, 
        LscorB  BLOB  
)
/
CREATE TABLE Proteins (
        protein_id INTEGER NOT NULL PRIMARY KEY AUTO_INCREMENT,
        pdbid   CHAR(5), 
        region  INTEGER
)
/
CREATE TABLE AlTypes (
        AlType_id INTEGER NOT NULL PRIMARY KEY AUTO_INCREMENT,
        name CHAR(24),
        params BLOB,
        pdbpath TEXT
)
/
CREATE TABLE PendingStrAls (
        Pendingstral_id INTEGER NOT NULL PRIMARY KEY AUTO_INCREMENT,
        protA INTEGER,
        protB INTEGER,
        altype INTEGER
)
/
';

=head1 Wurst::AlDb - persistence for wurst alignments

Class::DBI database classes for storing huge numbers of alignments.

=cut
