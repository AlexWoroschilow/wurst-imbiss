
package Wurst::AlDb::StrAlgn;

use base Wurst::AlDb;

Wurst::AlDb::StrAlgn->table("StrAlignments");
Wurst::AlDb::StrAlgn->columns(Primary=> qw/StrAlgn_id/);
Wurst::AlDb::StrAlgn->columns(Essential=> qw/protA protB ATotal BTotal T3 T4 T5 SeqId Lal Lalpfrg/);
Wurst::AlDb::StrAlgn->columns(Pairs => qw/fragset covA covB LscorA LscorB/);
Wurst::AlDb::StrAlgn->columns(PairsDetails => qw/A2B B2A/);
Wurst::AlDb::StrAlgn->columns(Details => qw/ANwTot ASwTot ASwCov ASwSc1 ASwSc2 ASwSqG ASwStG BNwTot BSwTot BSwCov BSwSc1 BSwSc2 BSwSqG BSwStG A2B B2A altype/);

use Compress::Zlib;
sub LscorA {
    my $this = shift;
    if (@_) {
        $this->set("LscorA"=>( Storable::nfreeze( shift @_ ) ));
    } else {
        return Storable::thaw($this->get("LscorA"));
    }
}

sub LscorB {
    my $this = shift;
    if (@_) {
        $this->set("LscorB"=>(Storable::nfreeze( shift @_ )));
    } else {
        return Storable::thaw($this->get("LscorB"));
    }
}

Wurst::AlDb::StrAlgn->set_sql(count_ofaltype => 'SELECT COUNT(*) FROM StrAlignments WHERE altype = ? ' );
sub count_ofaltype {
    my $this = shift;
    my $altype = shift;
    if (defined ($altype) ) {
        return $this->sql_count_ofaltype($altype->AlType_id())->select_val;
    }
    $this->_carp("count_ofaltype needs an altype id!\n");
    return undef;
}

use Carp qw(cluck croak);

sub _croak {
    my ($self, $message, %info) = @_;
    croak($message)
        unless $message=~/^Can't insert .*/; #';
    cluck ($message);
    return;
}



# find a better way for multiple aligntypes
#sub Count_AlTypeSet {
#    my $this = shift;
#    my @ret_vals;
#    if (@_) {
#        my $altype = shift @_;
#        while ($altype->isa("Wurst::AlDb::AlType")) {
#            Wurst::AlDb::StrAlgn->set_sql("SELECT COUNT(*) FROM ".$this->table()." WHERE altype=".$altype->Altype_id());
              
#Wurst::AlDb::StrAlgn->has_a(LscorA => 'UNIVERSAL', inflate => 
#                            sub { return bless Storable::thaw( shift @_ ); },
#                            deflate => 
#                            sub { return Storable::nfreeze { shift @_ }; }
#                            );

#Wurst::AlDb::StrAlgn->has_a(LscorB => 'UNIVERSAL', inflate => 
#                            sub { return bless Storable::thaw( shift @_ ); },
#                            deflate => 
#                            sub { return Storable::nfreeze { shift @_ }; }
#                            );

#Wurst::AlDb::StrAlgn->has_many(Alignments=>"Wurst::AlDb::Alignment");
Wurst::AlDb::StrAlgn->has_a(protA=>Wurst::AlDb::Protein);
Wurst::AlDb::StrAlgn->has_a(protB=>Wurst::AlDb::Protein);

Wurst::AlDb::StrAlgn->has_a(altype=>Wurst::AlDb::AlType);
Wurst::AlDb::StrAlgn->has_a(A2B => Wurst::AlDb::Alignment);
Wurst::AlDb::StrAlgn->has_a(B2A => Wurst::AlDb::Alignment);

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
# coveragestr-Seq-to-Str
# coveragestr-Str-to-Seq
# an Alignment
# Al-seqid
# Geom (T3,T4,T5)
# Score-Rescore (rs-paramtype, rs-val)
# Score-Total
# Sw-score-tot
# Sw-score+rs
# Nw-score-tot
# Nw-score+rs
# 


# id protA-protB
# protA-protB-Aldata
# protB-protA-Aldata
# Conserved-Aldata
# Conserved-Fragset

# Distancemeasurevalue (0-many) [ dmtype, value ]

# protein data derived from relations
# cluster-membership [ clustering id, cluster ]
# conservation pattern
