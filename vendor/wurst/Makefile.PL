# top level Makefile
use ExtUtils::MakeMaker;

$LIB_EXT = ".so";

WriteMakefile(
    'NAME'     => 'Wurst',
    'clean' => {FILES => "*% *.bak *.bck Base/*% out*.au"},
    'linkext' => { LINKTYPE => 'dynamic'},
#	'LDDLFLAGS' => "-lglut -lGL -lGLU -lminigsl",
	'LDDLFLAGS' => "-lminigsl",
#	'LIBS' => '-lGL -lGLU',
	'INC'  => "-I/usr/include",
    'dist'     => { COMPRESS => 'gzip -f9', SUFFIX => '.gz' },
    'VERSION_FROM'  => "lib/Wurst.pm",
#	'OPTIMIZE' => "-O3 -m64 -fPIC"
	'OPTIMIZE' => "-g"
);


