package Salamisrvini;

use strict;
use warnings;

use vars qw($VERSION @ISA @EXPORT @EXPORT_OK $AUTOLOAD);

require Exporter;
require DynaLoader;
require AutoLoader;


@ISA = qw(Exporter DynaLoader);
$VERSION = '0.01';

@EXPORT = qw($LIB_BASE
             $LIB_LIB
             $LIB_ARCH
             $PVEC_STRCT_DIR
             $PVEC_CA_DIR
             $CA_CLASSFILE
             $CLASSFILE
             $LOCAL_STORAGE
             $LOG_BASE
             $TOP_TEMP
             $JOBDIR
             $RESDIR
             $JSCRIPTSPATH
             $RESURL
             $OUTURL
             @DFLT_STRUCT_DIRS
             $BIN_SUFFIX
             $MPI_LIB_BASE
             $INPUT_CLST_LIST
             $PDB_TOP_DIR
             $OUTPUT_BIN_DIR
             $OUTPUT_LIB_LIST
             @GUNZIP
             $TMPDIR
             $STRUCT_DIR
             );

# bootstrap Salamisrvini $VERSION;

use vars qw($LIB_BASE
            $LIB_LIB
            $LIB_ARCH
            $PVEC_STRCT_DIR
            $PVEC_CA_DIR
            $CA_CLASSFILE
            $CLASSFILE
            $LOCAL_STORAGE
            $LOG_BASE
            $TOP_TEMP
            $JOBDIR
            $RESDIR
            $JSCRIPTSPATH
            $RESURL
            $OUTURL
            @DFLT_STRUCT_DIRS
            $BIN_SUFFIX
            $MPI_LIB_BASE
            $INPUT_CLST_LIST
            $PDB_TOP_DIR
            $OUTPUT_BIN_DIR
            $OUTPUT_LIB_LIST
            @GUNZIP
            $TMPDIR
            );

#*LIB_BASE  = \'/smallfiles/public/bm/salamiServer/v02/wurst2';
*LIB_BASE = \'/home/other/wurst/salamiServer/v02/wurst';
*LIB_LIB   = \"$LIB_BASE/blib/lib";
*LIB_ARCH  = \"$LIB_BASE/blib/arch";

*MPI_LIB_BASE = \'/work/other/wurst/salamiServer/update2/Parallel-MPI-Simple-0.10';   #@TODO


*PVEC_STRCT_DIR  = \'/work/public/no_backup/bm/pdb_all_vec_6mer_struct';#'/smallfiles/public/no_backup/bm/pdb_all_vec_6mer_struct';    #initialize in local Salamisrvini.pm;
*PVEC_CA_DIR    = \'/work/public/no_backup/bm/pdb90_vec_7mer_ca_mod_new';#/smallfiles/public/no_backup/bm/pdb90_vec_7mer_ca_mod_new';  #initialize in local Salamisrvini.pm;
*CA_CLASSFILE  = \'/work/public/no_backup/bm/F7_ca_mod';#/smallfiles/public/no_backup/bm/F7_ca_mod';        #initialize in local Salamisrvini.pm;

*CLASSFILE     = \'/work/public/no_backup/bm/classfile';#/smallfiles/public/no_backup/bm/classfile';        #initialize in local Salamisrvini.pm;


# This is a top level temporary directory. We can have a cron job
# run around and clean it up at night. Each job creates its own
# temporary directory under this one.
*TOP_TEMP   = \'/home/other/wurst/wurst_delete_able_temp';

*RESDIR   = \'/home/other/wurst/salamiServer/results/';  #test testetstes
*OUTURL = \'/home/other/wurst/public_html/salami/results/jobs';
*RESURL = \'http://flensburg.zbh.uni-hamburg.de/~wurst/salami/results';
# Where we will write logs to
*LOG_BASE   = \'log';
*JSCRIPTSPATH = \'/home/other/wurst/public_html/salami/results/';

*DFLT_STRUCT_DIRS = ['/work/public/no_backup/bm/pdb_all_bin', '.' ];
*STRUCT_DIR       = \'/work/public/no_backup/bm/pdb_all_bin';
*BIN_SUFFIX       = \'.bin';
$LOCAL_STORAGE    = '/dev/shm'; # This one may be modified

#==========================================================
##   pdb and vectors update settings
##==========================================================
*INPUT_CLST_LIST= \'/work/public/no_backup/pdb/derived_data/NR/clusters90.txt';
*PDB_TOP_DIR    = \'/work/public/no_backup/pdb/data/structures/divided/pdb';
*OUTPUT_BIN_DIR = \'/work/public/no_backup/bm/pdb_all_bin';
*OUTPUT_LIB_LIST= \'/work/public/no_backup/bm/pdb_lib.list';
*GUNZIP = ['/usr/bin/gunzip', '-f'];
*TMPDIR = \'/work/other/wurst/salamiServer/update2/TEMP';
