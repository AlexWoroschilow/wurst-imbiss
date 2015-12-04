# 6 Jan 2002
# Run a huge pile of scripts in this directory
# rcsid = $Id: all.pl,v 1.1 2007/09/28 16:57:02 mmundry Exp $

use FindBin;
use strict;


# ----------------------- run_one     -------------------------------
sub
run_one ($)
{
    my $f = shift;
    my $forker = fork();
    if ($forker == 0) {         # child
        open (STDOUT, "/dev/null") || die "redirecting STDOUT; $! ";
        my $result = do $f;
        if ($result != 0) {
            exit (0);
        } else {
            print STDERR "child: BROKE $f got $result\n";
            exit (1);
        }
    } else {
        my $x = wait;
        if ($? != 0) {
            print STDERR "parent: Fail on $f, returned $?\n";
        } else {
            print "parent: ok $f\n" }
    }

}


# ----------------------- file_fun    -----------------o--------------
# Return true if this file looks like it will be fun to
# execute. That is, it is a .pl file and not this routine.
sub
file_fun ($)
{
    my $f = shift;
    my $ret;
    if ( $f =~ m/\.pl$/) {
        $ret = 1; }
    else {
        return 0; }
    $0 =~ s/.*\///;
    $f =~ s/.*\///;
    if ( $f eq $0) {
        return 0;}
    return 1;
}

# ----------------------- main_bit    -------------------------------
sub
main_bit ()
{

    use FindBin;
    my $orig_dir = 'foo';

    my $orig_dir = `pwd`;
    my $dir = $FindBin::Bin;
    chdir $dir || die "Cant cd to $dir: $!\n";
    opendir (DIR,$dir) || die "Dir open fail on $dir: $!";


    while (my $f = readdir (DIR)) {
        if (file_fun ($f)) {
            run_one ($f); } }
 
    closedir (DIR);

    chdir $orig_dir;
    return 0;
}

# ----------------------- main        -------------------------------
exit main_bit();

