use v5.12.1;

my %nodes;
my $nodeA;
my $nodeB;
my $i = 0;
my $TMscore;
my $rmsd;
open INPUT, "<$ARGV[0]";
open LABELS, ">labels.txt";
while (<INPUT>){
    given($_){
        when(/^\/smallfiles\/public\/no_backup\/bm\/pdb_all_vec_6mer_struct\/([0-9a-z]{4}[A-Z0-9])/) { $nodes{$1} = $i++;
                                          print LABELS "$i\t$1\n";}
    } 
}
close INPUT;
open INPUT, "<$ARGV[0]";
while (<INPUT>){
    chomp $_;
    given($_){
      #when(/>([0-9a-z]{4}[A-Z0-9])$/) { $node = $1; 
        when(/^\/smallfiles\/public\/no_backup\/bm\/pdb_all_vec_6mer_struct\/([0-9a-z]{4}[A-Z0-9])/) { 
            $nodeA = $1; 
            $i = 1;
            $TMscore = 0;
        }
        when(/[0-9]+\s+\/smallfiles\/public\/no_backup\/bm\/pdb_all_vec_6mer_struct\/([0-9a-z]{4}[A-Z0-9a-z])\s+/) {
            $nodeB = $1;
        }
        when(/rmsd: ([0-9]+.[0-9]+), /){
            $rmsd = $1;
        }
        when(/tm_scr ([0-9]+.[0-9]+)/){
            $TMscore = $1;
        }
        when(/^1-exp\(-HSS/){
            unless(($nodeA eq $nodeB) or ($TMscore eq 0) ){
                #say "$nodes{$node}\t$nodes{$1}\t$i";
                #say "$nodeA\t$nodeB\t$rmsd\t$TMscore";
                say "$nodeA\t$nodeB\t$TMscore";
                $i++;
            }
        }
    }
}

