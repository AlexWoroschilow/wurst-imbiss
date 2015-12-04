# 7 June 2003
# This function is now used by a few different applications,
# so it wants its own file.  It gets its defaults from another
# package (Al_defaults).
# $Id: Pair_err.pm,v 1.1 2007/09/28 16:57:02 mmundry Exp $
{
use strict;
package Pair_err;
use vars qw($VERSION @ISA @EXPORT);
require Exporter;
@ISA = qw(Exporter);
@EXPORT = qw (pair_err);

# ----------------------- activate     ------------------------------
# This is modelled after the activation (logistic) function
# used in neural nets.
sub activate ($)
{
    my $x = shift;
    use Al_defaults qw ($SWITCH_POINT $X_SCALE);
    my $offset = $SWITCH_POINT;
    my $xscl = $X_SCALE;
    my $a = $offset - $x;
    $a *= $xscl;
    $a = exp ($a);
    return ( 1.0 / ( 1.0 + $a));
}

# ----------------------- pair_err     ------------------------------
# For one pair of structures, get the error.
# The fourth argument, if set, means that we weight the alignment
# according to how much of the smaller structure it accounts for. This
# will gently favour larger models.
# Unfortunately, the argument has to be optional because older scripts
# did not know about this.
sub pair_err ( $ $ $ ; $)
{
    my ($pair_set, $c1, $c2, $size_wt) = @_;
    use Wurst qw(coord_get_seq make_model coord_size dme_thresh coord_name);
    my $tmp_c = make_model ($pair_set, coord_get_seq ($c1), $c2);
    use Al_defaults qw ($MIN_MODEL);
    if ( !$tmp_c) {
        return 0; }
    my $smallsize = coord_size ($c1);
    {
        my $size2 = coord_size ($c2);
        if ($size2 < $smallsize) {
            $smallsize = $size2; }
    }
    my $model_size = coord_size ($tmp_c);
    my $model_frac = $model_size / $smallsize;
    if (( $model_frac ) < $MIN_MODEL) {
        return 0;
    }


    my $fraction = 0.0;
    if ( ! dme_thresh ( $fraction, $tmp_c, $c1, 4.0)) {
        print STDERR "dme error comparing coordinates from..\nModel ",
        coord_name ($tmp_c), "\n";
        print " of size ", coord_size ($tmp_c), " and \n";
        print coord_name ($c1), " of size ", coord_size ($c1), "\n";
        return 0;
    }

    my $x = - activate ($fraction);
    if ( defined ($size_wt) ) {
        $x *= $model_frac; }
    return ($x);
}
}
1;
