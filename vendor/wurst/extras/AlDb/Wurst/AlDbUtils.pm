# AlDbUtils.pm
# 1-June 2004
# This is an overlay module, without its own package space.
# it adds added sugar to the existing Wurst::AlDb classes.
# hence the flat structure and lack of an EXPORT.

use strict;

package Wurst::AlDb;



# do stuff with results
sub add_StrAlgn( $ $ $ @ @) {
#    my $self = shift @_;
    my ($altype, $chainA, $chainB, $legset, $resset) = @_;
    my (%dathash, %A2B, %B2A);
    while (scalar @$legset) {
        my $l = shift @$legset;
        my ($htb,$key) = ($l=~m/(.2.)\.(.+)/);
        if (defined ($htb)) {
            (Wurst::AlDb->verbosity()>1) and (print "Placing $l in table $htb { $key } \n");
            ($htb=~/A2B/) and $A2B{$key} = ($l=~/scor/) ? (Storable::nfreeze(shift @$resset)) : (shift @$resset);
            ($htb=~/B2A/) and $B2A{$key} = ($l=~/scor/) ? (Storable::nfreeze(shift @$resset)) : (shift @$resset);
        } else {
            (Wurst::AlDb->verbosity()>1) and print "Placing $l in StrAlgn\n";
            $dathash{$l} = ($l=~/scor/) ? (Storable::nfreeze(shift @$resset)) : (shift @$resset);
        }
    };
    $A2B{seq}=($chainA);
    $A2B{str}=($chainB);
    $A2B{AlType}=$altype;
    $B2A{seq}=($chainB);
    $B2A{str}=($chainA);
    $B2A{AlType}=$altype;
    my ($a2b,$b2a, $al_str_obj);
    $a2b = Wurst::AlDb::Alignment->create(\%A2B);
    (defined($a2b)) and ($b2a = Wurst::AlDb::Alignment->create(\%B2A));
    (defined($b2a)) and ($al_str_obj = Wurst::AlDb::StrAlgn->create({protA=>$chainA, protB=>$chainB, altype=>$altype, (%dathash), A2B=>$a2b, B2A=>$b2a} )); # interpolation works ?
    if (defined($al_str_obj)) {
        $a2b->StrAlgn_id($al_str_obj);
        $b2a->StrAlgn_id($al_str_obj);
        $a2b->update;
        $b2a->update;
    } else {
        (defined($b2a)) and ($b2a->delete());
        (defined($a2b)) and ($a2b->delete());
    }
    return $al_str_obj;
}

# validation conditions

package Wurst::AlDb::AlType;
# Protein entry can only exist if at least one altype provides a path to find its data (.bin file,
# and .prof).

sub valid_pdbid {
    my $altype = shift;
    my ($binname) = shift;
    defined($binname) or return undef;
    my ($coord, $prof);
    unless ($altype->isa("Wurst::AlDb::AlType")) { die "need an AlType object!"; };
    my $params=Storable::thaw($altype->params);
    
    if (not ((-f ($altype->pdbpath)."$binname.bin") and ($coord=Wurst::coord_read(($altype->pdbpath)."$binname.bin")))) {
        warn("Can't find $binname.bin in ".($altype->pdbpath)."\n");
    }
    if (exists $$params{profiles_path}) {
        if ((-f $$params{profiles_path}."$binname.prof") 
            and ($prof = Wurst::blst_chk_read($$params{profiles_path}."$binname.prof"))) {
            (defined $coord) and return ((Wurst::AlDb::Protein->find_or_create(pdbid=>$binname, region=>0)), $coord,$prof);
        } else {
            warn("Can't find $binname.prof in ".($$params{profiles_path})."\n");
        }
    } elsif ($coord) {
        return ((Wurst::AlDb::Protein->find_or_create(pdbid=>$binname, region=>0)), $coord);
    };
    return undef;
};

1;

