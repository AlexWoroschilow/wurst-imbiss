# 27 March 2002
# For testing pdb reading routines.
# $Id: pdb_read.pl,v 1.2 2010/12/01 13:26:20 mosisch Exp $

use FindBin;
use lib "$FindBin::Bin/../blib/arch";
use lib "$FindBin::Bin/../blib/lib";

use strict;
use POSIX qw(EXIT_SUCCESS EXIT_FAILURE);

use Wurst;

# ----------------------- constants ---------------------------------
use vars qw ($pdbdir @pdbfiles @chains);

sub set_constants ()
{
    $pdbdir = '/projects/biodata/pdb/data/structures/divided/pdb';
    @pdbfiles = ('pdb5pgm.ent.Z', '../struct/pdb1lea.ent.Z', 'pdb2pgk.ent.Z');
    @chains   = ('C',                 ' ',          '_');
}

# ----------------------- path_from_pdb -----------------------------
sub path_from_pdb ($)
{
    my ($pdbfile) = @_;
    my ($vol, $dir, $fname) = File::Spec->splitpath( $pdbfile );
    my $twocode = $fname;
    $twocode =~ s/^pdb//;
    $twocode = substr ($twocode, 1, 2);
    my $path = $pdbdir;
    if ( ! ($pdbfile =~ /\//)) {   # pdb file does not have leading path
        $path = "${path}/${twocode}/";
        $path = "${path}/${pdbfile}";
        $path = File::Spec->canonpath ($path);
    } else {
        $path = $pdbfile;
    }
    if (! -f $path) {
        print STDERR "Cannot find $path\n";
        return undef;
    } else {
        return ($path);
    }
}

# ----------------------- get_uncompressed --------------------------
# If the file seems to be compressed, uncompress it to a temp file.
# Return a name that can be deleted. Either a temporary file, or just
# make a link to the original file.
sub get_uncompressed ($)
{
    my $pdbfile = shift;
    my $gunzip = 'gunzip';
    my $temppdb = "delete_me.$$";
    if ( ! -f $pdbfile ) {
        die "Cannot find $pdbfile as pdbfile"; }
    if (( $pdbfile =~ /\.gz$/) || ($pdbfile =~ /\.Z$/ )) {
        use File::Copy;
        if ( -f $temppdb) {
            unlink ($temppdb); }
        copy ( $pdbfile, "$temppdb.gz");    
        my @cmdline = ('gunzip', "$temppdb.gz");
        system (@cmdline) && die "$gunzip returned ", $? >> 8, ". Stopping";
    } else {
        symlink ($pdbfile, $temppdb) || die "link failed. $!";
    }
    return $temppdb;
}

# ----------------------- one_file ----------------------------------
sub one_file ($ $)
{
    use File::Basename;
    $DB::single = 1;
    my ($pdbfile, $chain) = @_;
    my $temppdb = get_uncompressed ($pdbfile);
    print "temppdb is $temppdb\n";
    my $mdl = pdb_read ($temppdb, '', $chain) || die "pdb_read fail: $!";
    unlink ($temppdb) || die "unlink fail: $!";
    my ($file, $path, $ext) = fileparse ($pdbfile, '\..+$');
    my $out = "del_me_${file}.pdb";
    print "Read ", coord_size($mdl), " residues OK.\n";
    coord_2_pdb ($out, $mdl) || die "Failed writing $out";
}

# ----------------------- mymain  -----------------------------------
sub mymain ()
{
    set_constants();
    func_int; # force loading of wurst
#   $SIG{TRAP} = 'IGNORE';# kill 'TRAP', $$;
    for (my $i = 0; $i < @pdbfiles; $i++) {
        my $p = $pdbfiles[$i];
        my $chain = $chains[$i];
        my $path = path_from_pdb ($p);
        if ( ! defined ($path)) {
            print STDERR "Giving up on $p\n";
            next;
        }
        print "Now doing $path $chain\n";
        one_file ($path, $chain);
    }
    return EXIT_SUCCESS;
}
exit mymain();
