package Wurst;
# rcsid = "$Id: Wurst.pm,v 1.14 2015/09/24 12:31:33 wurst Exp $"
use strict;
use warnings;
use Carp;

package Wurst;
use vars qw($VERSION @ISA @EXPORT @EXPORT_OK $AUTOLOAD);

require Exporter;
require DynaLoader;
require AutoLoader;


@ISA = qw(Exporter DynaLoader );
# Items to export into callers namespace by default. Note: do not export
# names by default without a very good reason. Use EXPORT_OK instead.
# Do not simply export all your public functions/methods/constants.

# Please please delete "str_to_class" XXXX do not forget !!!

@EXPORT = qw(
        addvec
        addseq
        aa_2_prob_vec
        aa_strct_2_prob_vec
        aa_strct_clssfcn_read
        aa_strct_dump
        aa_strct_nclass
        aa_strct_size
        att_size_calpha
        ac_read_calpha 
        ac_size_calpha
        ac_nclass_calpha
        ac_class_wt_calpha
        ac_dump
        ac_nclass
        ac_read
        ac_size
        bit_set_bitset
        bit_set_calloc
        bit_set_clean
        bit_set_free
        bit_set_getsize
        bit_set_ifbitset
        blst_chk_read
        calpha_clssfcn_return
        calpha_clssfcn_destroy
        calpha_strct_2_prob_vec
        calpha_strct_ComputeMembership
        check_data_read
        check_cmpnd
        cleanveclist
        cleanseqlist
        conserv_str
        computeMembership
        computeMembershipAA
        computeMembershipAAProf
        coord_get_sec_s
        coord_get_numbering
        coord_has_sec_s
        coord_ca_ca_dist
        coord_c_n_dist
        coord_calc_psi
        coord_calc_phi
        coord_calc_omega
        coord_calc_tau
        coord_calc_theta
        coord_psi
        coord_phi
        coord_omega
        coord_tau
        coord_theta
        coord_copy
        coord_deletion
        coord_geo_gap
        coord_get_seq
        coord_name
        coord_name_thr
        coord_name_thr_destroy
        coord_read
        coord_rmsd
        coord_size
        coord_2_bin
        coord_2_pdb
        coord_2_pnlty
        coord_2_spdb
        deleteveclist
        deleteseqlist
        dme_thresh
        dpt_get_n
        dpt_get_val
        dpt_read
        dpt_set_val
        dpt_string
        initveclist
        initseqlist 
        fragmt_2_prob
        find_alt_path_score_simple
        file_pre_read
        file_un_read
        func_int
        func_float
        func_char
        funcs1_char
        funcs2_char
        get_apdb
        get_clssfcn
        get_conserv
        get_float_value
        get_focf
        get_nbor
        get_rmsd
        getconservvec
        get_solv_data
        get_torsion_data
        get_seq_id
        get_seq_id_simple
        get_seq_conserv
        get_threshs 
        get_pair_set_m
        get_pair_set_n
        pair_set_get_strNum
        make_model
        model_pdb_num        
        merge_alignments
        merge_align
        merge_localigns
        model_res_num
        multal_string
        split_multal
        paa_strct_clss6
        pair_set_chimera
        pair_set_coverage
        pair_set_gap
        pair_set_get_startpos
        pair_set_string
        pair_set_stringI
        pair_set_stringN        
        pair_set_pretty_string
        pair_set_print_indices
        pair_set_score
        pair_set_extend
        pair_set_sel_geti
        pair_set_print_prepare
        selected_pair_set_get
        pair_set_sel_delmaxdistance
        pair_set_sel_print
        pair_set_sel_printXXXXX
        pair_set_sel_len
        printstr
        pair_set_dest
        pair_set_get_alignment_indices
        pair_set_get_alignmentlength
        pair_set_get_netto_alignmentlength
        pair_set_get_residue
        param_fx_read
        param_rs_read
        pdb_read
        printconserv        
        prob_vec_info
        prob_vec_read        
        prob_vec_size
        prob_vec_length
        prob_vec_copy
        prob_vec_write
        prob_vec_destroy
        prof_aa_2_prob_vec
        prof_aa_strct_2_prob_vec
        ReadRescoreParam
        remove_seq
        remove_gaps
        score_mat_add
        score_mat_info
        score_mat_new
        score_mat_read
        score_mat_remove_alignment
        score_mat_scale
        score_mat_shift
        score_mat_string
        score_mat_sum_smpl
        score_mat_sum_sec
        score_mat_sum_full
        score_mat_write        
        score_metric
        sec_s_data_read
        sec_s_data_string
        seq_from_string
        seq_get_1
        seq_num
        seq_print_many
        seq_print
        seq_read
        seq_read_many
        seq_size
        seqprof_get_seq
        seqprof_str
        scor_set_fromvec
        scor_set_scale
        scor_set_simpl
        scor_set_to_vec
        score_dpt
        score_fx
        score_fx_prof
        score_prof_prof
        score_pvec
        score_pvec_gl
        score_rs
        score_sec
        score_smat
        score_sprof
        str_to_class
        strct_2_prob_vec
        struct_2_prob_vec
        sub_mat_read
        sub_mat_string
        sub_mat_shift
        sub_mat_scale
        sub_mat_get_by_c
        sub_mat_get_by_i
        sub_mat_set_by_c
        sub_mat_set_by_i
        tm_score
        tm_score_s
        write_fpc_bin
        wurstli_hello
        $EXT_LONG
        $EXT_SHORT
        $N_AND_W
        $S_AND_W
);
$VERSION = '0.01';
use vars qw($EXT_LONG $EXT_SHORT $N_AND_W $S_AND_W $N_O_R);
*EXT_LONG  = \0;
*EXT_SHORT = \1;
*N_AND_W   = \0;
*S_AND_W   = \1;

bootstrap Wurst $VERSION;

# Preloaded methods go here.
sub END {
    free_scratch();
}

use strict;

# ----------------------- score_mat_sum_smpl ------------------------
# OK Mr Smartypants. There are two very nasty tricks here.
#   1. The first argument will be modified. $rmat is a matrix to be
#      passed back to the interpreter.
#   2. On calling the score_m.. function, we pass in $$rmat, not $rmat.
# After exiting, $rmat will have been mutated into the correct class.
#   More tricks
# * The C code expects some float arrays, but usually we do not
#   want to use this feature. The magic below with the $null variable
#   is simply to pass NULL pointers to the C code.
# * The last arguement is a "pair_set" - the result of a previous
#   alignment calculation. It is optional, hence the two versions of
#   of the call.
sub score_mat_sum_smpl (\$ $ $ $ $ $ $; $)
{
    my ($rmat, $scores, $pgap_open, $pgap_widen, $qgap_open, $qgap_widen,
        $align_type, $bias_set) = @_;
    my $null = 0;
    bless (\$null, 'floatPtr');
    my $ret;
    if ( ! defined ($bias_set)) {
        $ret =
            score_mat_sum_full ($$rmat, $scores, $pgap_open, $pgap_widen,
                                $qgap_open, $qgap_widen, \$null, \$null,
                                $align_type);

    } else {
        $ret =
            score_mat_sum_full ($$rmat, $scores, $pgap_open, $pgap_widen,
                                $qgap_open, $qgap_widen, \$null, \$null,
                                $align_type, $bias_set);
    }
    return $ret;
}

sub score_mat_sum_smpl_circlpermut (\$ $ $ $ $ $ $ $ $; $)
{
    my ($rmat, $scores, $pgap_open, $pgap_widen, $qgap_open, $qgap_widen,
        $align_type, $pvec_len1, $pvec_len2, $bias_set) = @_;
    my $null = 0;
    bless (\$null, 'floatPtr');
    my $ret;
    if ( ! defined ($bias_set)) {
        $ret =
            score_mat_sum_full_circlpermut ($$rmat, $scores, $pgap_open, $pgap_widen,
                                $qgap_open, $qgap_widen, \$null, \$null,
                                $align_type, $pvec_len1, $pvec_len2);

    }
    
    return $ret;
}


# ----------------------- wurstli_hello      ------------------------
sub wurstli_hello
{
    print "wurstli wurstli\n";
}

# ----------------------- score_mat_sum_sec  ------------------------
sub score_mat_sum_sec ($ $ $ $ $ $ $ $ $; $)
{
    my ($rmat, $scores, $coord, $sec_pnlty,
        $pgap_open, $pgap_widen, $qgap_open, $qgap_widen,
        $align_type, $bias_set) = @_;
    my $null = 0;
    bless (\$null, 'floatPtr');

    if (! defined ($sec_pnlty)) {
        print STDERR "score_mat_sum_sec: called with \"sec_pnlty\" undef\n";
        print STDERR "score_mat_sum_sec: setting to 1.0\n";
        $sec_pnlty = 1;
    }
    my $mult = coord_2_pnlty ($coord, $sec_pnlty);
    my $ret;
    if ( ! defined ($bias_set)) {
        $ret =
            score_mat_sum_full ($rmat, $scores, $pgap_open, $pgap_widen,
                                $qgap_open, $qgap_widen, $mult, \$null,
                                $align_type);
    } else {
       $ret =
            score_mat_sum_full ($rmat, $scores, $pgap_open, $pgap_widen,
                                $qgap_open, $qgap_widen, $mult, \$null,
                                $align_type, $bias_set);
    }
    return $ret;
}


# Autoload methods go after =cut, and are processed by the autosplit program.

1;
__END__
# Below is the stub of documentation for your module. You better edit it!


=head1 NAME

Wurst - Perl extension for playing with alignment methods

=head1 SYNOPSIS

  use Wurst;
  blah blah blah

=head1 DESCRIPTION

The full description is in the file, F<wurst.pod>.

=head1 SEE ALSO

perl(1).

=cut
