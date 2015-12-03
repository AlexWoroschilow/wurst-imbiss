# SeqStrCmp.pm
# comparison routines that end up filling a database with protein similarity info
# 3rd June 2004
# This file was
# Seq_Str_Cmp.pm
# 20th May 2002
# Bundled Wurst function for computing measure of similarity
# between two sequence/structure pairs based on a given 
# forcefield and parameterization.

#!/usr/bin/perl 

use FindBin;
# use lib "/home/500/jxp500/Wurst/wurst/src/Wurst";
use lib "$FindBin::Bin/../src/Wurst/blib/arch";
use lib "$FindBin::Bin/../src/Wurst/blib/lib";

use strict;
use POSIX qw(EXIT_SUCCESS EXIT_FAILURE);

use Wurst;

{
  # * Local near-Constants
  # * - these can be overridden with new
  # * - values but would ideally contain
  # * - the latest parameter set
  
  use vars qw ($XSQ_PARAM_FILE $XSQ_MATR_NAME  $VERBOSITY);
  *XSQ_PARAM_FILE  = \"$FindBin::Bin/../../params/g9.param";
  *XSQ_MATR_NAME = \"$FindBin::Bin/../../matrix/blosum62.50";
  $VERBOSITY = 1;
  sub Set_Verbosity ( $ ) {
    $VERBOSITY = shift @_;
  }
  use vars qw ($OLD_FORMAT_OP);
  $OLD_FORMAT_OP = 2; # Always give alignment string back in results. - still current (ie version 1-2 - scor_sets not included by default
  sub Set_Old_Format_Op ( $ ) {
    $OLD_FORMAT_OP = shift @_;
  }
  

  # \"
  # From ~torda/bridget/one_align.pl 24th April 2002
  # ----------------------- Defaults ----------------------------------
  # These may look demented, (too many significant digits) but they
  # are pasted from an optimisation run.
  use vars qw (
	       $XSQ_align_type
	       $XSQ_sw1_pgap_open
	       $XSQ_sw1_qgap_open
	       $XSQ_sw1_pgap_widen
	       $XSQ_sw1_qgap_widen
               $XSQ_wt1 $XSQ_wt2 $XSQ_wt3
               $XSQ_off1

	       $XSQ_sw2_pgap_open
	       $XSQ_sw2_qgap_open
	       $XSQ_sw2_pgap_widen
	       $XSQ_sw2_qgap_widen
	       $XSQ_sw2_sec_pnlty

	       );
  use vars qw (
	       $XSQ_nw_pgap_open
	       $XSQ_nw_qgap_open
	       $XSQ_nw_pgap_widen
	       $XSQ_nw_qgap_widen
	       $XSQ_nw_sec_pnlty
               $XSQ_seq_wt
               $XSQ_final_nw
	       );
  
# Now, some magic constants for rescoring
  use vars qw (
	       $XSQ_NW_RS_SCR2 $XSQ_NW_GAP_GEO $XSQ_NW_SEQ_GAP $XSQ_NW_STR_GAP $XSQ_NW_STR_WDN
	       $XSQ_SW_RS_SCR2 $XSQ_SW_GAP_GEO $XSQ_SW_SEQ_GAP $XSQ_SW_STR_GAP $XSQ_SW_STR_WDN
	       );
  ( $XSQ_NW_RS_SCR2, $XSQ_NW_GAP_GEO, $XSQ_NW_SEQ_GAP, $XSQ_NW_STR_GAP, $XSQ_NW_STR_WDN) =
      ( 0.41587,     -3.27610,    -1.02167 ,   -11.52223,   -3.08488);
  ( $XSQ_SW_RS_SCR2, $XSQ_SW_GAP_GEO, $XSQ_SW_SEQ_GAP, $XSQ_SW_STR_GAP, $XSQ_SW_STR_WDN) =
      ( 0.47417,     -5.16621,    -1.35113,    -16.20915,  -6.57821);
# Params : 27th June
# for no secondary structure, two stage gotoh, final NW alignment.
        $XSQ_sw1_pgap_open  = 9.67385398757248;
        $XSQ_sw1_qgap_open  = 18.9017440741075;
        $XSQ_sw1_pgap_widen = 0.743565089047655;
        $XSQ_sw1_qgap_widen = 1.6537099328616;
        $XSQ_sw2_pgap_open  = 8.34897;
        $XSQ_sw2_qgap_open  = 18.9417;
        $XSQ_sw2_pgap_widen = 0.846406;
        $XSQ_sw2_qgap_widen = 1.75334;
        $XSQ_sw2_sec_pnlty  = 2.58267;

        $XSQ_nw_pgap_open  =  9.85439;
        $XSQ_nw_qgap_open  = 19.7665;
        $XSQ_nw_pgap_widen =  0.80877;
        $XSQ_nw_qgap_widen =  1.48625;
        $XSQ_nw_sec_pnlty  =  2.65248;
      
  # init values for profile mixing parameters
        $XSQ_seq_wt = 0.0;
  # and switches for the kinds of alignment
        $XSQ_final_nw = 1;
  # Local persistence for static data
  
  use vars qw ($XSQ_al_first $XSQ_FX_Parameters $XSQ_Sub_mat $XSQ_rescore_params $XSQ_rescore_paramsfile);
  use vars qw ($XSQ_params_protect);
  $XSQ_al_first = 1;
  $XSQ_FX_Parameters = undef;
  $XSQ_Sub_mat = undef;
  $XSQ_rescore_params = undef;
  $XSQ_rescore_paramsfile = undef;
  
  $XSQ_params_protect = 0; # this is set if the alignment routines should not change the current numbers

  sub x_sq_prot_alignparams {
      ($XSQ_params_protect = (shift @_)) if (@_);
      return $XSQ_params_protect;
  }
  sub x_sq_set_alignparams ( % ) {
      my $newparms = shift @_;
      (keys %{$newparms} == values %{$newparms}) or return EXIT_FAILURE;
      (exists $$newparms{PARAM_FILE}) and *XSQ_PARAM_FILE  = \ $$newparms{PARAM_FILE};
      (exists $$newparms{MATR_NAME}) and *XSQ_MATR_NAME = \ $$newparms{MATR_NAME};
      (exists $$newparms{align_type}) and $XSQ_align_type = $$newparms{align_type};
      (exists $$newparms{sw1_pgap_open}) and   $XSQ_sw1_pgap_open = $$newparms{sw1_pgap_open};
      (exists $$newparms{sw1_qgap_open}) and   $XSQ_sw1_qgap_open = $$newparms{sw1_qgap_open};
      (exists $$newparms{sw1_pgap_widen}) and   $XSQ_sw1_pgap_widen = $$newparms{sw1_pgap_widen};
      (exists $$newparms{sw1_qgap_widen}) and   $XSQ_sw1_qgap_widen =  $$newparms{sw1_qgap_widen};
      (exists $$newparms{sw2_pgap_open}) and   $XSQ_sw2_pgap_open = $$newparms{sw2_pgap_open};
      (exists $$newparms{sw2_qgap_open}) and   $XSQ_sw2_qgap_open = $$newparms{sw2_qgap_open};
      (exists $$newparms{sw2_pgap_widen}) and   $XSQ_sw2_pgap_widen = $$newparms{sw2_pgap_widen};
      (exists $$newparms{sw2_qgap_widen}) and   $XSQ_sw2_qgap_widen =  $$newparms{sw2_qgap_widen};
      (exists $$newparms{sw2_sec_pnlty})and $XSQ_sw2_sec_pnlty = $$newparms{sw2_sec_pnlty};
      (exists $$newparms{nw_pgap_open}) and   $XSQ_nw_pgap_open = $$newparms{nw_pgap_open};
      (exists $$newparms{nw_qgap_open}) and   $XSQ_nw_qgap_open = $$newparms{nw_qgap_open};
      (exists $$newparms{nw_pgap_widen}) and   $XSQ_nw_pgap_widen = $$newparms{nw_pgap_widen};
      (exists $$newparms{nw_qgap_widen}) and   $XSQ_nw_qgap_widen =  $$newparms{nw_qgap_widen};
      (exists $$newparms{nw_sec_pnlty})and $XSQ_nw_sec_pnlty = $$newparms{nw_sec_pnlty};
      (exists $$newparms{rescore_paramsfile}) and $XSQ_rescore_paramsfile = $$newparms{rescore_paramsfile};
      (exists $$newparms{final_nw}) and $XSQ_final_nw = $$newparms{final_nw};
      (exists ($$newparms{nw_rs_scr2})) and ($XSQ_NW_RS_SCR2 = $$newparms{nw_rs_scr2} );
      (exists ($$newparms{sw_rs_scr2})) and ($XSQ_SW_RS_SCR2 = $$newparms{sw_rs_scr2} );
      (exists ($$newparms{sw_gap_geo})) and ($XSQ_SW_GAP_GEO = $$newparms{sw_gap_geo} );
      (exists ($$newparms{sw_seq_gap})) and ($XSQ_SW_SEQ_GAP = $$newparms{sw_seq_gap} );
      (exists ($$newparms{sw_str_gap})) and ($XSQ_SW_STR_GAP = $$newparms{sw_str_gap} );
      (exists ($$newparms{sw_str_wdn})) and ($XSQ_SW_STR_WDN = $$newparms{sw_str_wdn} );
      (exists ($$newparms{nw_gap_geo})) and ($XSQ_NW_GAP_GEO = $$newparms{nw_gap_geo} );
      (exists ($$newparms{nw_seq_gap})) and ($XSQ_NW_SEQ_GAP = $$newparms{nw_seq_gap} );
      (exists ($$newparms{nw_str_gap})) and ($XSQ_NW_STR_GAP = $$newparms{nw_str_gap} );
      (exists ($$newparms{nw_str_wdn})) and ($XSQ_NW_STR_WDN = $$newparms{nw_str_wdn} );

      # Reset persistent data
      $XSQ_FX_Parameters = undef;
      $XSQ_Sub_mat = undef;
      $XSQ_rescore_params = undef;
      $XSQ_al_first = 1;
      return EXIT_SUCCESS;
  }
  sub x_sq_get_alignparams  {
      my $newparms = {};
      $$newparms{PARAM_FILE} = $XSQ_PARAM_FILE ;
      $$newparms{MATR_NAME} = $XSQ_MATR_NAME;
      $$newparms{align_type} = $XSQ_align_type;
      $$newparms{sw1_pgap_open} =   $XSQ_sw1_pgap_open;
      $$newparms{sw1_qgap_open} =   $XSQ_sw1_qgap_open;
      $$newparms{sw1_pgap_widen} =   $XSQ_sw1_pgap_widen;
      $$newparms{sw1_qgap_widen} =   $XSQ_sw1_qgap_widen;
      $$newparms{sw2_pgap_open} =   $XSQ_sw2_pgap_open;
      $$newparms{sw2_qgap_open} =   $XSQ_sw2_qgap_open;
      $$newparms{sw2_pgap_widen} =   $XSQ_sw2_pgap_widen;
      $$newparms{sw2_qgap_widen} =   $XSQ_sw2_qgap_widen;
      $$newparms{nw_pgap_open} =   $XSQ_nw_pgap_open;
      $$newparms{nw_qgap_open} =   $XSQ_nw_qgap_open;
      $$newparms{nw_pgap_widen} =   $XSQ_nw_pgap_widen;
      $$newparms{nw_qgap_widen} =   $XSQ_nw_qgap_widen;
      $$newparms{seq_wt} = $XSQ_seq_wt;
      $$newparms{nw_rs_scr2} = $XSQ_NW_RS_SCR2;
      $$newparms{sw_rs_scr2} = $XSQ_SW_RS_SCR2;
      $$newparms{sw_gap_geo} = $XSQ_SW_GAP_GEO; 
      $$newparms{sw_seq_gap} = $XSQ_SW_SEQ_GAP; 
      $$newparms{sw_str_gap} = $XSQ_SW_STR_GAP; 
      $$newparms{sw_str_wdn} = $XSQ_SW_STR_WDN;
      $$newparms{nw_gap_geo} = $XSQ_NW_GAP_GEO; 
      $$newparms{nw_seq_gap} = $XSQ_NW_SEQ_GAP; 
      $$newparms{nw_str_gap} = $XSQ_NW_STR_GAP; 
      $$newparms{nw_str_wdn} = $XSQ_NW_STR_WDN;

      $$newparms{rescore_paramsfile} = $XSQ_rescore_paramsfile;
      $$newparms{zero_shift} = $XSQ_zero_shift;
      $$newparms{final_nw} = $XSQ_final_nw;
      return $newparms;
  }
  
# From libsrch.pl
  sub get_scores ($ $ $ $ $) {
      my ($pair_set, $coord, $seq, $rescore_params, $to_use) = @_;
      
      my ($scr_tot, $coverage, $score1, $score2, $geo_gap,$too_small);
      my ($str1, $str2, $crap);
      my ($open_cost, $widen_cost, $nseq_gap);
      ($crap, $score1) = pair_set_score( $pair_set );
      # Also get the coverage onto the structure because we can recreate the
      # alignment that way.
      ($str1, $str2) =
          pair_set_coverage ($pair_set, seq_size ($seq), coord_size ($coord));
      $coverage = ($str1 =~ tr/1/1/);       # This is coverage as an integer
      $too_small = ($coverage<=2);          # trap a bug
      $coverage = $coverage / seq_size ($seq); #and as fraction of sequence
      
      my ($k_scr2, $k_gap_geo, $k_seq_gap, $k_str_gap, $k_str_wdn);
      if ($to_use eq 's_and_w') {
          ( $k_scr2,     $k_gap_geo,  $k_seq_gap,  $k_str_gap,  $k_str_wdn) =
              ( $XSQ_SW_RS_SCR2, $XSQ_SW_GAP_GEO, $XSQ_SW_SEQ_GAP, $XSQ_SW_STR_GAP, $XSQ_SW_STR_WDN);
      } else {
          ( $k_scr2,     $k_gap_geo,  $k_seq_gap,  $k_str_gap,  $k_str_wdn) =
              ( $XSQ_NW_RS_SCR2, $XSQ_NW_GAP_GEO, $XSQ_NW_SEQ_GAP, $XSQ_NW_STR_GAP, $XSQ_NW_STR_WDN);
      }

      if ($too_small || ($coverage  < .05 )) {
          $score2 = 0;
          $geo_gap = 0;
          $nseq_gap = 0;
          $open_cost = 0;
          $widen_cost = 0;
      } else {
          my $model = make_model ($pair_set, $seq, $coord); # Build model
          $score2 = score_rs ($model, $rescore_params);
          ($crap, $geo_gap, $crap, $nseq_gap) = coord_geo_gap ($model, 1, 10);
          ($open_cost, $widen_cost) = pair_set_gap($pair_set, 1, 1);
      }
      
      $scr_tot = $score1 +
          $k_scr2    * $score2 +
          $k_gap_geo * $geo_gap +
          $k_seq_gap * $nseq_gap +
          $k_str_gap * $open_cost +
          $k_str_wdn * $widen_cost;
      return ( $scr_tot, $coverage, $score1, $score2, $nseq_gap, $open_cost, $str1,$str2); # Return the coverage string, too.
  }
                  
                  
# ----------------------- do_align ----------------------------------
sub do_align ($ $ $ $ $ $)
{
    my ($seq, $coord, $sec_s_data, $fx_params, $submat, $rescore_params)
        = @_;
    
#   Now we have all files we need. Start scoring.
#   Begin by giving us four empty score matrices.
    my $sim_scr_mat =                                 # For similarity matrix
        score_mat_new (seq_size ($seq), seq_size (coord_get_seq ($coord)));
    my $sec_scr_mat =                                 # For sec struct matrix
	score_mat_new (seq_size ($seq), seq_size (coord_get_seq ($coord)));
    my $fx_scr_mat =                                  # For fx9 func scores
        score_mat_new (seq_size ($seq), seq_size (coord_get_seq ($coord)));
    my $total_scr_mat =                               # Where we put totals
        score_mat_new (seq_size ($seq), seq_size (coord_get_seq ($coord)));
    score_smat ($sim_scr_mat, $seq, coord_get_seq ($coord), $submat);
    if ($XSQ_sec_s_scale != 0.0) {
	score_sec  ($sec_scr_mat, $sec_s_data, $coord); }
    score_fx   ($fx_scr_mat, $seq, $coord, $fx_params);
    score_mat_shift ($fx_scr_mat, $XSQ_fx_mat_shift);
#   We have three score, matrices. Now add them together.
    $total_scr_mat = score_mat_add ($fx_scr_mat, $sim_scr_mat, $XSQ_sim_scale);
    if ($XSQ_sec_s_scale != 0.0) {
        $total_scr_mat =
            score_mat_add ($total_scr_mat, $sec_scr_mat,$XSQ_sec_s_scale); }
#   This actually does the alignment.
    my ($sw_scr_tot, $sw_coverage, $sw_score1, $sw_score2,
        $sw_seq_gap, $sw_strct_gap);
    my ($nw_scr_tot, $nw_coverage, $nw_score1, $nw_score2,
        $nw_seq_gap, $nw_strct_gap, $cover_string, $s_cover_string);

    my $sw_pair_set =
        score_mat_sum_sec (my $result_mat, $total_scr_mat,
                           $coord, $XSQ_sw1_sec_pnlty,
                           $XSQ_sw1_pgap_open, $XSQ_sw1_pgap_widen,
                           $XSQ_sw1_qgap_open, $XSQ_sw1_qgap_widen,
                           $S_AND_W);
    $sw_pair_set =
        score_mat_sum_sec (my $result_mat2, $total_scr_mat,
                           $coord, $XSQ_sw2_sec_pnlty,
                           $XSQ_sw2_pgap_open, $XSQ_sw2_pgap_widen,
                           $XSQ_sw2_qgap_open, $XSQ_sw2_qgap_widen,
                           $S_AND_W, $sw_pair_set);
    
    ($sw_scr_tot, $sw_coverage, $sw_score1, $sw_score2,
     $sw_seq_gap, $sw_strct_gap, $cover_string, $s_cover_string) =
         get_scores ($sw_pair_set, $coord, $seq, $rescore_params, 's_and_w');
    my $nw_pair_set =
            score_mat_sum_sec (my $result_mat3, $total_scr_mat,
                               $coord, $XSQ_nw_sec_pnlty,
                               $XSQ_nw_pgap_open, $XSQ_nw_pgap_widen,
                               $XSQ_nw_qgap_open, $XSQ_nw_qgap_widen,
                               $N_AND_W, $sw_pair_set);
    ($nw_scr_tot, $nw_coverage, $nw_score1, $nw_score2,
     $nw_seq_gap, $nw_strct_gap, $cover_string, $s_cover_string) =
         get_scores ($nw_pair_set, $coord, $seq, $rescore_params, 'n_and_w');
    
    my $scor_set = scor_set_simpl($nw_pair_set, $total_scr_mat); # annotation for the alignment
    
    my @r =
        ($nw_scr_tot, $nw_pair_set, $sw_scr_tot,
         $sw_coverage, $nw_coverage, $sw_score1, $sw_score2,
         $nw_score1, $nw_score2, $sw_pair_set, $total_scr_mat, $cover_string, $s_cover_string, $scor_set);
    
    return (\@r);
}


sub x_sq_struct_align ( $ $ )
{
    my ($coord1,$coord2) = @_;
    my ($param, $submat, $rescore_params);
    if ($XSQ_al_first) {
        $param = param_fx_read ($XSQ_PARAM_FILE) || die "Fail on $XSQ_PARAM_FILE";
        $submat = sub_mat_read ($XSQ_MATR_NAME) || die "Fail on $XSQ_MATR_NAME";
        $submat = sub_mat_shift ($submat, $XSQ_sim_mat_bottom);
        my $pfile;
        if (not defined($XSQ_rescore_paramsfile))  {
            $pfile = "$FindBin::Bin/../../params";
            $pfile = $pfile.'/cc_allat_p891+0.param';
            $XSQ_rescore_paramsfile = $pfile;
        }
        $rescore_params = param_rs_read ($XSQ_rescore_paramsfile) || die "Rescore params";
        $XSQ_FX_Parameters = $param;
        $XSQ_Sub_mat = $submat;
        $XSQ_rescore_params = $rescore_params;
        $XSQ_al_first = 0;
    } else {
        $param = $XSQ_FX_Parameters;
        $submat = $XSQ_Sub_mat;
        $rescore_params = $XSQ_rescore_params;
    }
    
    (not $coord1 or not $coord2) and return;
    
    my ($seq1,$seq2) = (coord_get_seq ($coord1), coord_get_seq($coord2));
    
    my $A1S2 = do_align( $seq1, $coord2, "nothing here", $param, $submat, $rescore_params);
    my $A2S1 = do_align( $seq2, $coord1, "nothing here", $param, $submat, $rescore_params);
    # returns all coverage strings and local score annotation
    ((scalar @$A1S2 != 14) || (scalar @$A2S1 != 14)) and return (2,"Alignments failed (AtoB, BtoA)", (scalar @$A1S2), (scalar @$A2S1));
    # Find Common and Coverage
    my (@arrset);
    @arrset = ps_cmp_and_model ($A1S2->[1], $A1S2->[10], $coord1, $A2S1->[1], $A2S1->[10], $coord2);
    (scalar @arrset != 6) and return (3, "Comparison of alignments failed");
    
    
    if ($VERBOSITY > 10) {
        # debug VERbosity
        print "Alignment A onto B\n".(pair_set_pretty_string $A1S2->[1],$seq1,$seq2)."\n";
        print "Alignment B onto A\n".(pair_set_pretty_string $A2S1->[1],$seq2,$seq1)."\n";
        print "Conserved alignment: A,B :\n".(pair_set_pretty_string $arrset[5],$seq1,$seq2)."\n";      
    };
    
    my $BsmodelofA = make_model($A1S2->[1], $seq1, $coord2);
    my $AsmodelofB = make_model($A2S1->[1], $seq2, $coord1);
    
    #Construct Result Set
    # From :
    # Symmetric Relations, A-Symmetric Relations
    #  \_> |C| |CId|       \_> SeqAStructB, SeqBStructA
    #      |CperFrag|           |->   NW_Score_Total, nw Coverage, swScrTotal, sw Coverage
    #      Cons_Score_Smpl            Geom(Model, OrigStructure)
    #      Geom(ConsA,ConsB)                
    #      Alignment (Coverage, fraglist)
    # Decide on the results set
    #  - practical point - conserved coverage string could be used to filter alignments to 'unique' areas of each representative
    #  - original pairsets are useful?
    my ($thresh, @res, $temp);
    
    # Fraction aligned and average frag length
    my $frags = $arrset[4]=~tr/)/)/; # Number of Fragments
my $av_cons_al 
    = 2.0*coord_size($arrset[1])/(coord_size($coord1)+coord_size($coord2)); # Consistent Rali
my $av_cons_perfrag = $av_cons_al;
($frags>0) and ($av_cons_perfrag = $av_cons_al/$frags);
my $num_cons = coord_size($arrset[1]);
my ($cons_cov_stringA,$cons_cov_stringB);
($cons_cov_stringA, $cons_cov_stringB) = 
    pair_set_coverage ($arrset[5], coord_size ($coord1), coord_size ($coord2));
my $cons_SqId = 0.0; #### !!! JBPNote
my $cons_ps_B2A = pair_set_xchange $arrset[5];
# Compute Rescores for each conserved model structure
my ($sc_cons_b, $sc_cons_a) = pair_set_score($arrset[5]);
my ($ABcswt,$ABcnwt, $ABcnwc, $ABnws1, $ABnws2, $ABnwsgap, $ABnwstrgp);
my ($BAcswt,$BAcnwt, $BAcnwc, $BAnws1, $BAnws2, $BAnwsgap, $BAnwstrgp,$crap);
($ABcnwt, $ABcnwc, $ABnws1, $ABnws2, $ABnwsgap, $ABnwstrgp,$crap,$crap) = 
    get_scores ($arrset[5], $coord2, $seq1, $rescore_params, 'n_and_w'); # We expect to have 'fragment like' alignments
($ABcswt, $ABcnwc, $ABnws1, $ABnws2, $ABnwsgap, $ABnwstrgp,$crap,$crap) = 
    get_scores ($arrset[5], $coord2, $seq1, $rescore_params, 's_and_w'); # We expect to have 'fragment like' alignments
($BAcnwt, $BAcnwc, $BAnws1, $BAnws2, $BAnwsgap, $BAnwstrgp,$crap,$crap) = 
    get_scores ($cons_ps_B2A, $coord1, $seq2, $rescore_params, 'n_and_w');
($BAcswt, $BAcnwc, $BAnws1, $BAnws2, $BAnwsgap, $BAnwstrgp,$crap,$crap) = 
    get_scores ($cons_ps_B2A, $coord1, $seq2, $rescore_params, 's_and_w');
#      if ($VERBOSITY>10) {
#	print "n_and_w rescores :\n\tA Side : ".(join ",",(get_scores ($arrset[5], $coord2, $seq1, $rescore_params, 'n_and_w')))."\n";
#	print "\tB Side : ".(join ",",(get_scores ($cons_ps_B2A, $coord1, $seq2, $rescore_params, 'n_and_w')))."\n";
#      }

# Scores from each alignment
# nw_scr_tot, .. , sw_scr_tot, sw_coverage, nw_coverage, sw_score1, sw_score2, nw_score1, nw_scoer2, .., ..;
my @Sqa_scoreset = @$A1S2;
my @Sqb_scoreset = @$A2S1;
@Sqa_scoreset = @Sqa_scoreset[0,2,3,4,5,6,7,8];
@Sqb_scoreset = @Sqb_scoreset[0,2,3,4,5,6,7,8];

# Fill out result set :
# Overall :  Numbers (useful for clustering) followed by Alignment details (Fragments, coverage strings, etc)
# Symmetric relations first
push @res, $av_cons_al, $av_cons_perfrag;#, $cons_SqId;
# Geometric Similarity Scores
# Seq-Str Consistent Residues
foreach $temp (3.0, 4.0, 5.0) {
    (coord_size($arrset[1])>=2) and (dme_thresh $thresh, $arrset[1],$arrset[2],$temp);
    (coord_size($arrset[1])<2) and (($thresh = 0.0),1.0);
    push @res, $thresh;
}
# Relations for Seq(A)->Struct(B)
push @res, @Sqa_scoreset, $ABcnwt, $ABcswt, $ABcnwc, $ABnws1, $ABnws2, $ABnwsgap, $ABnwstrgp;
# Deformation of As structure when sequence is aligned to structure B
if (coord_size($arrset[0])>=2) {
    foreach $temp (3.0, 4.0, 5.0) {
        dme_thresh $thresh, $arrset[0],$BsmodelofA,$temp;
        push @res, $thresh;
    }
} else {
    push @res, 0.0,0.0,0.0;
};

# for Seq B to Str A
push @res, @Sqb_scoreset, $BAcnwt, $BAcswt, $BAcnwc, $BAnws1, $BAnws2, $BAnwsgap, $BAnwstrgp;
# Deformation of Bs structure when sequence is aligned to structure A
foreach $temp (3.0, 4.0, 5.0) {
    (coord_size($arrset[3])>=2) and (dme_thresh $thresh, $arrset[3],$AsmodelofB,$temp);
    (coord_size($arrset[3])<2) and (($thresh = 0.0),1.0);
    push @res, $thresh;
}

# Stings and things
push @res, $arrset[4],$cons_cov_stringA,$cons_cov_stringB,$A1S2->[11],$A2S1->[11];
(push @res, $A1S2->[12],$A2S1->[12]) if ($OLD_FORMAT_OP>1);
(push @res, scor_set_to_vec($A1S2->[13]),scor_set_to_vec($A2S1->[13])) if ($OLD_FORMAT_OP>2); # we don't normally put these in yet.
return (@res);

}
  
sub x_sq_struct_failed ( $ $ ) {
    my $silent = shift @_;
    my $aset = shift @_;
    
    if (scalar (@$aset) <= 3) {
        # print an error code if we can work it out
        (scalar (@$aset) == 3) and (not $silent) and (print STDERR "Seq_Str_Cmp.pm: ".$aset->[1]." ".$aset->[2]."\n");
        (not scalar @$aset) and (not $silent) and (print STDERR "Seq_Str_Cmp.pm: unknown error!!!!");
        return 't';
    };
    
    return undef;
}

# From the wurst_server script may31st2004
sub set_may31prof_alignparams {
    my %prof_params;
    %prof_params = (
                    "align_type" => 'profile_fx',
                    "sw1_pgap_open" => 1.213,
                    "sw1_pgap_widen" => 0.194,
                    "sw1_qgap_open"  =>  2.791,    
                    "sw1_qgap_widen" =>  0.944,
                    "sec_s_scale"    =>   0.0,
                    "fx_mat_shift"       =>   1.995,
                    "seq_wt"         =>   0.814,
                    "zero_shift"     =>   0.185,
                    
                    "nw_pgap_open"  => 1.213,
                    "nw_qgap_open"  =>    2.791,
                    "nw_pgap_widen" =>    0.194,
                    "nw_qgap_widen" =>    0.944,
                    "MATR_NAME" => "$FindBin::Bin/../../matrix/708tmp.mat",
                    "rescore_paramsfile" => "$FindBin::Bin/../../params/cc_allat_p891+0.param"
                    );
    x_sq_set_alignparams ( \%prof_params );
    # rescoring seems to be the same as ever!
}
# From the wurst_server script 10th June 2004
sub set_prof_params {
# From the wurst_server script 10th June 2004
    my %prof_params = (
                       "align_type" => 'prof_fx_jun10_2004',
                       "sw1_pgap_open" => 1.140,
                       "sw1_pgap_widen" => 0.190,
                       "sw1_qgap_open"  =>  2.733,    
                       "sw1_qgap_widen" =>  0.952,
                       "sec_s_scale"    =>   0.0,
                       "fx_mat_shift"       =>   1.987,
                       "seq_wt"         =>   0.830,
                       "zero_shift"     =>   0.163,
                       
                       "nw_pgap_open"  => 1.140,
                       "nw_qgap_open"  =>    2.733,
                       "nw_pgap_widen" =>    0.190,
                       "nw_qgap_widen" =>    0.952,
                       "nw_rs_scr2" => .13,
                       "sw_rs_scr2" => .18,
                       "sw_gap_geo" => -1.12,
                       "sw_seq_gap" => -0.57,
                       "sw_str_gap" => -2.18,
                       "sw_str_wdn" => -0.21,
                       "nw_gap_geo" => -1.27,
                       "nw_seq_gap" => -0.09,
                       "nw_str_gap" => -2.35,
                       "nw_str_wdn" => -0.93,
#                   "sim_mat_bottom" => 
                       "MATR_NAME" => "/home/procter/src/prj/wurst/matrix/719.out",
                       "rescore_paramsfile" => "/home/procter/src/prj/Cleaning/wurst/params/cc_allat_p891+0.param",
                       "PARAM_FILE" => "/home/procter/src/prj/Cleaning/wurst/params/g9.param"
                       );
    x_sq_set_alignparams(\%prof_params);
    

    
}
sub set_swprof_params {
    set_prof_alignparams();
    my $params = x_sq_get_alignparams();
    x_sq_set_alignparams ( { "align_type" => "sw_".$$params{align_type},
                             "final_nw" => 0 });
}

# ----------------------- zero_shift_mat ----------------------------
# This shifts a substitution matrix, not an alignment score matrix.
sub zero_shift_mat ($ $)
{
    my ($sub_mat, $shift) = @_;
    for (my $i = 0; $i < 20; $i++) {
        for (my $j = $i; $j < 20 ; $j++) {
            my $t = sub_mat_get_by_i ($sub_mat, $i, $j) + $shift;
            sub_mat_set_by_i ($sub_mat, $i, $j, $t);
        }
    }
}

sub do_prof_align ($ $ $ $ $ $ $) {
    my ($profile, $coord, $tmplt_prof, $sec_s_data,  $fx_params, $submat, $rescore_params)
        = @_;
    
    my $seq = seqprof_get_seq ($profile);
    my $seq_size   = seq_size ($seq);
    my $coord_size = coord_size ($coord);
    my $fx_mat     =  score_mat_new ($seq_size, $coord_size);
    my $sim_scr_mat   = score_mat_new ($seq_size, $coord_size);
    my $seq_strct_mat = score_mat_new ($seq_size, $coord_size);
    $DB::single = 1;
    score_fx_prof   ($fx_mat, $profile, $coord, $fx_params);
    $fx_mat = score_mat_scale ($fx_mat, $XSQ_fx_mat_shift);
    $fx_mat = score_mat_shift ($fx_mat, $XSQ_fx_mat_shift);

    my $sim_scr_mat =                                 # For similarity matrix
        score_mat_new (seq_size ($seq), seq_size (coord_get_seq ($coord)));
    score_prof_prof ($sim_scr_mat, $profile, $tmplt_prof, $submat);

    my $fx_wt = 1.0 - $XSQ_seq_wt;      # Weights of sequence and structure terms
    my $total_mat =                                  # Where we put totals
        score_mat_new (seq_size ($seq), seq_size (coord_get_seq ($coord)));

#   We have two score, matrices. Now add them together.
    $total_mat = score_mat_add ($total_mat, $fx_mat, $fx_wt);
    $total_mat = score_mat_add ($total_mat, $sim_scr_mat, $XSQ_seq_wt);


#   This actually does the alignment.
    my ($sw_scr_tot, $sw_coverage, $sw_score1, $sw_score2,
        $sw_seq_gap, $sw_strct_gap, $sw_covstring, $sw_s_covstring);
    my ($nw_scr_tot, $nw_coverage, $nw_score1, $nw_score2,
        $nw_seq_gap, $nw_strct_gap, $cover_string, $s_cover_string);

    my $sw_pair_set =
        score_mat_sum_smpl (my $result_mat, $total_mat,
                            $XSQ_sw1_pgap_open, $XSQ_sw1_pgap_widen,
                            $XSQ_sw1_qgap_open, $XSQ_sw1_qgap_widen,
                            $S_AND_W);

    ($sw_scr_tot, $sw_coverage, $sw_score1, $sw_score2,
     $sw_seq_gap, $sw_strct_gap, $sw_covstring, $sw_s_covstring) =
         get_scores ($sw_pair_set, $coord, $seq, $rescore_params, 's_and_w');

    my $nw_pair_set;
    if ($XSQ_final_nw) {
        $nw_pair_set =
            score_mat_sum_sec (my $result_mat2, $total_mat,
                               $coord, $XSQ_nw_sec_pnlty,
                               $XSQ_nw_pgap_open, $XSQ_nw_pgap_widen,
                               $XSQ_nw_qgap_open, $XSQ_nw_qgap_widen,
                               $N_AND_W, $sw_pair_set);
        ($nw_scr_tot, $nw_coverage, $nw_score1, $nw_score2,
         $nw_seq_gap, $nw_strct_gap, $cover_string, $s_cover_string) =
             get_scores ($nw_pair_set, $coord, $seq, $rescore_params, 'n_and_w');
    } else {
        $nw_pair_set = $sw_pair_set;
        ($nw_scr_tot, $nw_coverage, $nw_score1, $nw_score2,
         $nw_seq_gap, $nw_strct_gap, $cover_string, $s_cover_string)
            = ($sw_scr_tot, $sw_coverage, $sw_score1, $sw_score2,
               $sw_seq_gap, $sw_strct_gap, $sw_covstring, $sw_s_covstring)
        }
    my $scor_set = scor_set_simpl($nw_pair_set, $total_mat); # annotation of local scores
    my @r =
        ($sw_scr_tot, $nw_pair_set, $nw_scr_tot,
         $sw_coverage, $nw_coverage, $sw_score1, $sw_score2,
         $nw_score1, $nw_score2, $sw_pair_set, $total_mat, $cover_string, $s_cover_string, $scor_set, $sw_covstring, $sw_s_covstring);
    return ([@r]);
}

sub calc_seqcons ( $ $ $ $ ) {
    my ($cov_s_A, $cov_s_B, $seqA, $seqB) = @_;
    return (0.0) unless ($cov_s_A=~/1/);
    my $SqId = 0.0;
    my $numAl = 0;
    my @c_covA = split '',$cov_s_A;
    my @c_resA = split '', $seqA;
    my @c_covB = split '',$cov_s_B;
    my @c_resB = split '', $seqB;
    my ($ca,$ra,$cb,$rb) = ('0','','0','');
    while ((scalar @c_covA) && (scalar @c_covB)) {
        while ((scalar @c_covA) and (not $ca)) {
            ($ca,$ra)  = ((shift @c_covA), (shift @c_resA));
        }
        while ((scalar @c_covB) and (not $cb)) {
            ($cb,$rb)  = ((shift @c_covB), (shift @c_resB));
        }
        if ($ca and $cb) {
            $numAl++;
            if ($ra eq $rb) { $SqId++; };
            $ca=0;
            $cb=0;
        }
    }
    $SqId/=$numAl;
};


sub x_prof_str_align ( $ $ $ $ ) {
    my ($profA, $coordA, $profB, $coordB) = @_;
    (not $coordA or not $coordB or not $profA or not $profB) 
        and (warn "null parameters for x_prof_str_align",return);
    # use same semaphores for rereading coords or not
    my ($param, $submat, $rescore_params);
    if ($XSQ_al_first) {
        # set parameters correctly if not already done
        if (not $XSQ_params_protect) { set_prof_alignparams (); };
        $param = param_fx_read ($XSQ_PARAM_FILE) || die "Fail on $XSQ_PARAM_FILE";
        $submat = sub_mat_read ($XSQ_MATR_NAME) || die "Fail on $XSQ_MATR_NAME";
        zero_shift_mat ($submat, $XSQ_zero_shift);
        my $pfile;
        if (defined($XSQ_rescore_paramsfile)) {
            $pfile = $XSQ_rescore_paramsfile;
        } else {
            $pfile = "$FindBin::Bin/../../params";
            $pfile = $pfile.'/cc_allat_p891+0.param';
        }
        $rescore_params = param_rs_read ($pfile) || die "Failed to read Rescore params $pfile";
        $XSQ_FX_Parameters = $param;
        $XSQ_Sub_mat = $submat;
        $XSQ_rescore_params = $rescore_params;
        $XSQ_al_first = 0;
    } else {
        $param = $XSQ_FX_Parameters;
        $submat = $XSQ_Sub_mat;
        $rescore_params = $XSQ_rescore_params;
    }
    
    my $seqA = seqprof_get_seq($profA);
    my $seqB = seqprof_get_seq($profB);
    
    
    my $AlAtoB = do_prof_align($profA,$coordB,$profB,"nothing here", $param, $submat, $rescore_params);
#    ((push $@,"Failed on an alignment : ".coord_name($coordA)." onto ".coord_name($coordB)."\n"), return undef) if (not (scalar @{$AlAtoB}));
    my $AlBtoA = do_prof_align($profB,$coordA,$profA,"nothing here", $param, $submat, $rescore_params);
#    ((push $@,"Failed on nx alignment : ".coord_name($coordB)." onto ".coord_name($coordA)."\n"), return undef) if (not (scalar @{$AlAtoB}));
    ((scalar @$AlAtoB != 16) || (scalar @$AlBtoA != 16)) and return (2,"Alignments failed (AtoB, BtoA)", (scalar @$AlAtoB), (scalar @$AlBtoA));
    # sequence id for alignments
    my $A2B_seqid = calc_seqcons( $$AlAtoB[11], $$AlAtoB[12], (seq_print $seqA), (seq_print $seqB));
    my $B2A_seqid = calc_seqcons( $$AlBtoA[11], $$AlBtoA[12], (seq_print $seqB), (seq_print $seqA));
    # Analyse consistency and stuff
    
    # Construct the results
    my @arrset = ps_cmp_and_model ($AlAtoB->[1], $AlAtoB->[10], $coordA, $AlBtoA->[1], $AlBtoA->[10], $coordB);
    (scalar @arrset != 6) and return (3, "Comparison of alignments failed");
    
    
    if ($VERBOSITY > 10) {
        # debug VERbosity
        print "Alignment A onto B\n".(pair_set_pretty_string $AlAtoB->[1],$seqA, $seqB)."\n";
        print "Alignment B onto A\n".(pair_set_pretty_string $AlBtoA->[1],$seqB, $seqA)."\n";
        print "Conserved alignment: A,B :\n".(pair_set_pretty_string $arrset[5],$seqA, $seqB)."\n";      
    };
    
    my $BsmodelofA = make_model($AlAtoB->[1], $seqA, $coordB);
    my $AsmodelofB = make_model($AlBtoA->[1], $seqB, $coordA);
    
    my ($thresh, @res, $temp);
    
    # Fraction aligned and average frag length
    (0) and (tr/(/(/); # stupid formatting!
                 my $frags = $arrset[4]=~tr/)/)/; # Number of Fragments (())
    my $av_cons_al 
        = 2.0*coord_size($arrset[1])/(coord_size($coordA)+coord_size($coordB)); # Consistent Rali
    my $av_cons_perfrag = $av_cons_al;
    ($frags>0) and ($av_cons_perfrag = $av_cons_al/$frags);
    
    my $num_cons = coord_size($arrset[1]);
    
    my ($cons_cov_stringA,$cons_cov_stringB);
    ($cons_cov_stringA, $cons_cov_stringB) = 
        pair_set_coverage ($arrset[5], coord_size ($coordA), coord_size ($coordB));
    
    my $cons_SqId=0.;
    if ($num_cons) {
        $cons_SqId = calc_seqcons( $cons_cov_stringA, $cons_cov_stringB, (seq_print $seqA), (seq_print $seqB));
    }
    
    my $cons_ps_B2A = pair_set_xchange $arrset[5];
# Compute Rescores for each conserved model structure
    my ($sc_cons_b, $sc_cons_a) = pair_set_score($arrset[5]);
    my ($ABcswt,$ABcnwt, $ABcnwc, $ABnws1, $ABnws2, $ABnwsgap, $ABnwstrgp);
    my ($BAcswt,$BAcnwt, $BAcnwc, $BAnws1, $BAnws2, $BAnwsgap, $BAnwstrgp,$crap);
    ($ABcnwt, $ABcnwc, $ABnws1, $ABnws2, $ABnwsgap, $ABnwstrgp,$crap,$crap) = 
        get_scores ($arrset[5], $coordB, $seqA, $rescore_params, 'n_and_w'); # We expect to have 'fragment like' alignments
    ($ABcswt, $ABcnwc, $ABnws1, $ABnws2, $ABnwsgap, $ABnwstrgp,$crap,$crap) = 
        get_scores ($arrset[5], $coordB, $seqA, $rescore_params, 's_and_w'); # We expect to have 'fragment like' alignments
    ($BAcnwt, $BAcnwc, $BAnws1, $BAnws2, $BAnwsgap, $BAnwstrgp,$crap,$crap) = 
        get_scores ($cons_ps_B2A, $coordA, $seqB, $rescore_params, 'n_and_w');
    ($BAcswt, $BAcnwc, $BAnws1, $BAnws2, $BAnwsgap, $BAnwstrgp,$crap,$crap) = 
        get_scores ($cons_ps_B2A, $coordA, $seqB, $rescore_params, 's_and_w');
#      if ($VERBOSITY>10) {
#	print "n_and_w rescores :\n\tA Side : ".(join ",",(get_scores ($arrset[5], $coord2, $seq1, $rescore_params, 'n_and_w')))."\n";
#	print "\tB Side : ".(join ",",(get_scores ($cons_ps_B2A, $coord1, $seq2, $rescore_params, 'n_and_w')))."\n";
#      }
    
# Scores from each alignment
# nw_scr_tot, .. , sw_scr_tot, sw_coverage, nw_coverage, sw_score1, sw_score2, nw_score1, nw_scoer2, .., ..;
    my @Sqa_scoreset = @$AlAtoB;
    my @Sqb_scoreset = @$AlBtoA;
    @Sqa_scoreset = @Sqa_scoreset[0,2,3,4,5,6,7,8];
    @Sqb_scoreset = @Sqb_scoreset[0,2,3,4,5,6,7,8];
    
# Fill out result set :
# Overall :  Numbers (useful for clustering) followed by Alignment details (Fragments, coverage strings, etc)
# Symmetric relations first
    push @res, $av_cons_al, $av_cons_perfrag; 
    (push @res, $cons_SqId) if ($OLD_FORMAT_OP>3); # Only put in for very latest (2004)
# Geometric Similarity Scores
# Seq-Str Consistent Residues
    foreach $temp (3.0, 4.0, 5.0) {
        (coord_size($arrset[1])>=2) and (dme_thresh $thresh, $arrset[1],$arrset[2],$temp);
        (coord_size($arrset[1])<2) and (($thresh = 0.0),1.0);
        push @res, $thresh;
    }
# Relations for Seq(A)->Struct(B)
    push @res, @Sqa_scoreset, $ABcnwt, $ABcswt, $ABcnwc, $ABnws1, $ABnws2, $ABnwsgap, $ABnwstrgp;
# Deformation of As structure when sequence is aligned to structure B
    if (coord_size($arrset[0])>=2) {
        foreach $temp (3.0, 4.0, 5.0) {
            dme_thresh $thresh, $arrset[0],$BsmodelofA,$temp;
            push @res, $thresh;
        }
    } else {
        push @res, 0.0,0.0,0.0;
    };
    
# for Seq B to Str A
    push @res, @Sqb_scoreset, $BAcnwt, $BAcswt, $BAcnwc, $BAnws1, $BAnws2, $BAnwsgap, $BAnwstrgp;
# Deformation of Bs structure when sequence is aligned to structure A
    foreach $temp (3.0, 4.0, 5.0) {
        (coord_size($arrset[3])>=2) and (dme_thresh $thresh, $arrset[3],$AsmodelofB,$temp);
        (coord_size($arrset[3])<2) and (($thresh = 0.0),1.0);
        push @res, $thresh;
    }
    
# Stings and things
    push @res, $arrset[4],$cons_cov_stringA,$cons_cov_stringB;
    (push @res, $AlAtoB->[14],$AlBtoA->[14],$AlAtoB->[15],$AlBtoA->[15]) if ($OLD_FORMAT_OP>4); # the sw-coverages.
    push @res, $AlAtoB->[11],$AlBtoA->[11];
    (push @res, $AlAtoB->[12],$AlBtoA->[12]) if ($OLD_FORMAT_OP>1);
    (push @res, scor_set_to_vec($AlAtoB->[13]),scor_set_to_vec($AlBtoA->[13])) if ($OLD_FORMAT_OP>2); # we don't normally put these in yet.
    (push @res, $A2B_seqid, $B2A_seqid) if ($OLD_FORMAT_OP>4); # for the database - completeness
    return (@res);
}

sub x_sq_struct_legend () {
    (return (
             "Lal(Cons)",
             "LalPerFrg",
             "T3(Cons)",
             "T4(Cons)",
             "T5(Cons)",
             "A->BTotal",
             "A->BSwTot",
             "A->BSwCov",
             "A->BNwCov",
             "A->BSwSc1",
             "A->BSwSc2",
             "A->BNwSc1",
             "A->BNwSc2",
             "C(A)NwTot",
             "C(A)SwTot",
             "C(A)CovSW",
             "C(A)SwSc1",
             "C(A)SwSc2",
             "C(A)SwSqG",
             "C(A)SwStG",
             "T3(A->B)",
             "T4(A->B)",
             "T5(A->B)",
             
             "B->ATotal",
             "B->ASwTot",
             "B->ASwCov",
             "B->ANwCov",
             "B->ASwSc1",
             "B->ASwSc2",
             "B->ANwSc1",
             "B->ANwSc2",
             "C(B)NwTot",
             "C(B)SwTot",
             "C(B)CovSW",
             "C(B)SwSc1",
             "C(B)SwSc2",
             "C(B)SwSqG",
             "C(B)SwStG",
             "T3(B->A)",
             "T4(B->A)",
             "T5(B->A)",
             "Cfragset",
             "covConsA",
             "covConsB",
             
             "covAligA",
             "covAligB"
             )) if ($OLD_FORMAT_OP==1);
    
    (return (
            "Lal(Cons)",
            "LalPerFrg",
            "T3(Cons)",
            "T4(Cons)",
            "T5(Cons)",
            "A->BTotal",
            "A->BSwTot",
            "A->BSwCov",
            "A->BNwCov",
            "A->BSwSc1",
            "A->BSwSc2",
            "A->BNwSc1",
            "A->BNwSc2",
            "C(A)NwTot",
            "C(A)SwTot",
            "C(A)CovSW",
            "C(A)SwSc1",
            "C(A)SwSc2",
            "C(A)SwSqG",
            "C(A)SwStG",
            "T3(A->B)",
            "T4(A->B)",
            "T5(A->B)",
            
            "B->ATotal",
            "B->ASwTot",
            "B->ASwCov",
            "B->ANwCov",
            "B->ASwSc1",
            "B->ASwSc2",
            "B->ANwSc1",
            "B->ANwSc2",
            "C(B)NwTot",
            "C(B)SwTot",
            "C(B)CovSW",
            "C(B)SwSc1",
            "C(B)SwSc2",
            "C(B)SwSqG",
            "C(B)SwStG",
            "T3(B->A)",
            "T4(B->A)",
            "T5(B->A)",
            "Cfragset",
            "covConsA",
            "covConsB",
            
            "covAligA",
            "covAligB",
            "covBstrA",
            "covAstrB")) if ($OLD_FORMAT_OP==2);
    
    return (
            "Lal(Cons)",
            "LalPerFrg",
            "T3(Cons)",
            "T4(Cons)",
            "T5(Cons)",
            "A->BTotal",
            "A->BSwTot",
            "A->BSwCov",
            "A->BNwCov",
            "A->BSwSc1",
            "A->BSwSc2",
            "A->BNwSc1",
            "A->BNwSc2",
            "C(A)NwTot",
            "C(A)SwTot",
            "C(A)CovSW",
            "C(A)SwSc1",
            "C(A)SwSc2",
            "C(A)SwSqG",
            "C(A)SwStG",
            "T3(A->B)",
            "T4(A->B)",
            "T5(A->B)",
            
            "B->ATotal",
            "B->ASwTot",
            "B->ASwCov",
            "B->ANwCov",
            "B->ASwSc1",
            "B->ASwSc2",
            "B->ANwSc1",
            "B->ANwSc2",
            "C(B)NwTot",
            "C(B)SwTot",
            "C(B)CovSW",
            "C(B)SwSc1",
            "C(B)SwSc2",
            "C(B)SwSqG",
            "C(B)SwStG",
            "T3(B->A)",
            "T4(B->A)",
            "T5(B->A)",
            "Cfragset",
            "covConsA",
            "covConsB",
            
            "covAligA",
            "covAligB",
            "covBstrA",
            "covAstrB",
            "locBstrA",
            "locAstrB") if ($OLD_FORMAT_OP==3);

    return (
            "Lal",
            "Lalpfrg",
            "SeqId",
            "T3",
            "T4",
            "T5",
            "A2B.Total",
            "A2B.SwTot",
            "A2B.SwCov",
            "A2B.NwCov",
            "A2B.SwSc1",
            "A2B.SwSc2",
            "A2B.NwSc1",
            "A2B.NwSc2",
            "ANwTot",
            "ASwTot",
            "ASwCov",
            "ASwSc1",
            "ASwSc2",
            "ASwSqG",
            "ASwStG",
            "A2B.T3",
            "A2B.T4",
            "A2B.T5",
            
            "B2A.Total",
            "B2A.SwTot",
            "B2A.SwCov",
            "B2A.NwCov",
            "B2A.SwSc1",
            "B2A.SwSc2",
            "B2A.NwSc1",
            "B2A.NwSc2",
            "BNwTot",
            "BSwTot",
            "BSwCov",
            "BSwSc1",
            "BSwSc2",
            "BSwSqG",
            "BSwStG",
            "B2A.T3",
            "B2A.T4",
            "B2A.T5",
            "fragset",
            "covA",
            "covB",
            
            "A2B.Seqcovsw",# added
            "B2A.Seqcovsw",# .
            "A2B.Strcovsw",# .
            "B2A.Strcovsw",# .
            "A2B.Seqcov",
            "B2A.Seqcov",
            "A2B.Strcov",
            "B2A.Strcov",
            "A2B.Localscore",
            "B2A.Localscore", # ($OLD_FORMAT_OP==4.5);
            "A2B.SeqId",
            "B2A.SeqId")
        #  return (@leg);
    }
# Big Z = normal cstring
# Little z = compressed cstring (ie long redundant)
# f - float.
# need another. v - vector of floats.
sub x_sq_struct_packstring() {
    (return "fffff"."ffffffff"."fffffff"."fff"."ffffffff"."fffffff"."fff"."Zzzzz") if ($OLD_FORMAT_OP==1);
    (return "fffff"."ffffffff"."fffffff"."fff"."ffffffff"."fffffff"."fff"."Zzzzzzz") if ($OLD_FORMAT_OP==2);
    (return "fffff"."ffffffff"."fffffff"."fff"."ffffffff"."fffffff"."fff"."Zzzzzzzvv") if ($OLD_FORMAT_OP==3);
    (return "ffffff"."ffffffff"."fffffff"."fff"."ffffffff"."fffffff"."fff"."Zzzzzzzzzzzvvff") if ($OLD_FORMAT_OP==4);
    (return "ffffff"."ffffffff"."fffffff"."fff"."ffffffff"."fffffff"."fff"."Zzzzzzzzzzzvvff") if ($OLD_FORMAT_OP<5);

}
1; # to make the module be true
};
