# 13 March 2002
# Check functions for reading substitution matrices
# rcsid = $Id: sub_mat.pl,v 1.1 2007/09/28 16:57:01 mmundry Exp $

eval 'use warnings';
if ($@) {
    print "Warnings not present\n";}

use FindBin;

my $paths_place = 'paths.inc';
do $paths_place || die "Failed to read $paths_place: $!";
if ($@) { die "Compilation fail on $paths_place:\n$@" };

use vars qw ($MATRIX_DIR $MATRIX_FILE);


use strict;
my $tdir = "$FindBin::Bin/../src/Wurst";
*WURST_DIRS = ["$tdir/arch",
               "$tdir/lib" ];
# But after installation
# WURST_DIRS = [ "$ENV{HOME}/pl/lib" ];

use lib "$FindBin::Bin../blib/arch";
use lib "$FindBin::Bin../blib/lib";
#use lib "$ENV{HOME}";  # Where the wurst stuff lives after installation
use Wurst;

use POSIX qw(EXIT_SUCCESS EXIT_FAILURE);

# ----------------------- mymain       ------------------------------
sub mymain ()
{
    my $matname = $MATRIX_DIR . '/' . $MATRIX_FILE;
    my $subst_mat = sub_mat_read ($matname) || 
                   die "Fail getting sub mat $matname";
    print "sub mat starts at\n", sub_mat_string ($subst_mat), "\n";
    my $s2 = sub_mat_shift ($subst_mat, -2);
    print "New sub mat with min element -2, is \n",
    sub_mat_string ($s2);
    print "Another one at -4.5\n",
             sub_mat_string(sub_mat_shift($subst_mat, -4.5));
    return (EXIT_SUCCESS);
}

exit mymain();
