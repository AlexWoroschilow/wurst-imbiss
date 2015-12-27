CC=gcc
LD=${CC} 
CFLAGS= -O3 -g -DDONT_SEE_RCS -Wall -I./inc/ -I./vendor/cini/include/ -I./vendor/zlog/include/  -I./vendor/gnuplot/src/ -I./lib/ -I./vendor/wurst/src/wurstsrc/ -L./vendor/zlog/lib/ -L./vendor/cini/lib/ -L./vendor/wurst/src/wurstsrc/ -L./inc/ -L./lib/  
LDFLAGS= -lm -lwurst -lc -lpthread  -lzlog -lconfigini
CTAGS=ctags > tags
LIBS=-lob


WURSTIMBISSOBJ   = ./lib/memman.o\
	./lib/list.o\
	./lib/vtprogressbar.o\
	./lib/sort.o\
	./lib/cantor.o\
	./lib/stack.o\
	./lib/mathematics.o\
	./lib/stringutils.o\
	./lib/boyermoore.o\
	./lib/fileio.o\
	./lib/dpalign.c\
	./lib/hash.c\
	./inc/falphabet.o\
	./inc/createalphabet.o\
	./inc/encodeprobvec.o\
	./inc/intsequence.o\
	./inc/sufarray.o\
	./inc/depictseqs.o\
	./inc/mm.o\
	./inc/sufmatch.o\
	./inc/multiseq.o\
	./inc/salami.o\
	./inc/probseqscore.c\
	./inc/blaststat.c\
	./inc/imbissblast.c\
	./inc/imbiss_common.c\
	./inc/imbiss_getopt.c\
	./inc/logger.c\
	wurstimbiss.o\

all: wurstimbiss.x

wurstimbiss.x: ${WURSTIMBISSOBJ}
	gcc $(CFLAGS) ${WURSTIMBISSOBJ} -o $@ $(LDFLAGS)

clean: 
	rm ./inc/*.o
	rm ./lib/*.o
	rm *.o


