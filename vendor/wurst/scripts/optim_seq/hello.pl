# This script must run without errors before there is any
# serious chance of anything else running.
use strict;

package Alignfunc;
use vars qw($VERSION @ISA @EXPORT);
require Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(align_cost init_align_cost);

use lib "$ENV{HOME}/pl/lib";
use Wurst;
use POSIX qw(EXIT_SUCCESS EXIT_FAILURE);

print "Hello. \
It is now worth trying to run another script from this directory.\n";
