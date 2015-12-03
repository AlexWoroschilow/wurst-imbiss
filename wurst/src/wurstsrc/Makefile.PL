#Builds an archive (static library) containing Wurst symbols
use ExtUtils::MakeMaker;
use Config;

# Here you can add additional compiler flags
$additionalccflags="";
# additionalgccflags are only used if the cc is a gcc
$additionalgccflags="";

$ccflags="$Config{'ccflags'} $additionalccflags";
$depend="";
$dependpostamble="";
if ($Config{'ccname'} eq "gcc" ) {
 $ccflags="$ccflags $additionalgccflags";
 $dependtarget="depend::\n\trm -f depend\n\t$Config{'cc'} -MM *.c > depend\n";
 $dependpostamble="include depend";
}

$LIB_EXT = ".so";

WriteMakefile(
    'NAME'	=> 'LibWurst',
    'LINKTYPE'	=> 'dynamic',
    'LDDLFLAGS' => "--shared -L../../blib/arch/auto/LibMiniGSL/ -lMiniGSL",
#    'LDDLFLAGS' => "--shared -L/usr/lib64 -L../../blib/arch/auto/LibMiniGSL/ -lglut -lGL -lGLU -lMiniGSL",
    'LIBS' => '-lMiniGSL',
#    'LIBS' => '-lGL -lGLU -lMiniGSL',
    'INC'  => "-I/usr/include",
    'CCFLAGS'	=> "$ccflags",
    'OBJECT'	=> '$(O_FILES)',
    'clean'	=> {'FILES' => 'libwurst$(LIB_EXT) *~'},
    macro	=> {
		    ADD_AR_OPT	  => ($add_ar_flags || ''),
		    MY_CC_PRE_TARGET => ($Using_ms_vc
					 ? '-Fo'
					 : ($Using_Borland
					    ? '-o': '-o $(MY_EMPTY_STR)')),
		    MY_AR_PRE_TARGET => ($Using_ms_vc
					 ? '-out:'
					 : ($Using_Borland
					    ? '' : 'cr $(MY_EMPTY_STR)')),
		    MY_AR_OBJECT => ($Using_Borland
				     ? '$(OBJECT:^"+")'
				      : '$(OBJECT)'),
		    MY_EMPTY_STR  => '',
		   },
);

sub MY::top_targets {
'
all :: libwurst$(LIB_EXT)

tags:
	etags $(C_FILES)

splint:
	splint -preproc +posixlib -realcompare $(C_FILES)

static ::       libwurst$(LIB_EXT)

libwurst$(LIB_EXT): $(OBJECT) gsldir/Makefile
	-$(RM_F) libwurst$(LIB_EXT)
	$(AR) $(ADD_AR_OPT) $(MY_AR_PRE_TARGET)libwurst$(LIB_EXT) $(MY_AR_OBJECT)
	$(RANLIB) libwurst$(LIB_EXT)

'.$dependtarget.'
';
}



sub MY::postamble {$dependpostamble;}
