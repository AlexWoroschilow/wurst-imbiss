# Date
# what I do

use FindBin;
use lib "$FindBin::Bin../blib/arch";
use lib "$FindBin::Bin../blib/lib";

use strict;
use POSIX qw(EXIT_SUCCESS EXIT_FAILURE);

use Wurst;
# ----------------------- mymain  -----------------------------------
sub mymain ()
{
    return EXIT_SUCCESS;
}
exit mymain();
