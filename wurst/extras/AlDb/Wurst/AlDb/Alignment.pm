
package Wurst::AlDb::Alignment;
#use Scorset;
use base Wurst::AlDb;

use Storable qw(nfreeze thaw);

# House keeping

# what kind of alignment (profileStrStr, str/str, seqstr).
# what parameters were used (dep. on align type)

# which proteins form a square of relations
# - cached because it requires a number of counts otherwise


# The dataset normally stored in a stream file

# Protein map - index of proteins - arbitrary order

# Key formed from proteinA_proteinB where A<B
# score data

# localScore-Aseq-Bstr
# localScore-Bseq-Astr

# Alignments
# 

# id protA-protB
# protA-protB-Aldata
# protB-protA-Aldata
# Conserved-Aldata
# Distancemeasurevalue (0-many) [ dmtype, value ]

# protein data derived from relations
# cluster-membership [ clustering id, cluster ]
# conservation pattern

Wurst::AlDb::Alignment->table("Alignments");
Wurst::AlDb::Alignment->columns(Primary=> qw/Alignment_id/);
Wurst::AlDb::Alignment->columns(Essential=> qw/Total SwTot NwCov SwCov SwSc1 SwSc2 NwSc1 NwSc2 T3 T4 T5 SeqId/);
Wurst::AlDb::Alignment->columns(Pairs => qw/seq str Seqcov Strcov Seqcovsw Strcovsw Localscore/);
Wurst::AlDb::Alignment->columns(Details => qw/altype StrAlgn_id/);

# guess this is really the only way to do the storable hook now
# This is an ugly hack - it requires the use of $Alignment->Localscore rather than $Alignment->get(Localscore)
# not so terrible, but euggh.
use Compress::Zlib;
sub Localscore {
    my $this = shift;
    if (@_) {
        $this->set("Localscore"=>Storable::nfreeze( shift @_ ));
    } else {

        return Storable::thaw($this->get("Localscore"));
    }
}

use Carp qw(croak cluck);

sub _croak {
    my ($self, $message, %info) = @_;
    croak($message)
        unless $message=~/^Can't insert .*/; #';
    cluck ($message);
    return;
}


# Stuff that nearly works but leaves an array reference blessed as a universal class reference
# Wurst::AlDb::Alignment->has_a(Localscore=> 'UNIVERSAL',
#                               inflate => sub  { 
#                                   print "inflating\n".(join ",",(@_),"\n");
#                                   return bless Storable::thaw( shift @_ ); } ,
#                               deflate => sub  { 
#                                   print "deflating\n".(join ",",(@_),"\n");
#                                   return Storable::nfreeze( shift @_ ); }  );
# has many's

Wurst::AlDb::Alignment->has_a(seq=>Wurst::AlDb::Protein);
Wurst::AlDb::Alignment->has_a(str=>Wurst::AlDb::Protein);

Wurst::AlDb::Alignment->has_a(altype=>Wurst::AlDb::AlType);
Wurst::AlDb::Alignment->has_a(StrAlgn_id=>Wurst::AlDb::StrAlgn);


                                
#                                seqstr_id 
#seq str
#        altype  INTEGER,
        
# Stralgn_id INTEGER 
