CC=gcc
LD=${CC} 
CFLAGS= -O3 -DDONT_SEE_RCS -Wall -pedantic -I./inc/ -I./vendor/cini/include/ -I./vendor/zlog/include/  -I./vendor/gnuplot/src/ -I./lib/ -I./vendor/wurst/src/wurstsrc/ -L./vendor/zlog/lib/ -L./vendor/cini/lib/ -L./vendor/wurst/src/wurstsrc/ -L./inc/ -L./lib/  
LDFLAGS= -lm -lwurst -lc -lpthread  -lzlog -lconfigini
CTAGS=ctags > tags
LIBS=-lob



PROBVECTESTOBJ = ./inc/prob_vec.o\
	./inc/mprintf.o\
	./lib/memman.o\
	./lib/vtprogressbar.o\
	./inc/e_malloc.o\
	./inc/fio.o\
	./inc/matrix.o\
	./inc/mgc_num.o\
	./inc/scratch.o\
	./lib/sort.o\
	./lib/cantor.o\
	./lib/stack.o\
	./lib/mathematics.o\
	./lib/stringutils.o\
	./lib/fileio.o\
	./inc/falphabet.o\
	./inc/createalphabet.o\
	./inc/encodeprobvec.o\
	./inc/intsequence.o\
	probvectest.o
				


CONVERTPROBVECOBJ = ./inc/prob_vec.o\
	./inc/mprintf.o\
	./lib/memman.o\
	./lib/vtprogressbar.o\
	./inc/e_malloc.o\
	./inc/fio.o\
	./inc/matrix.o\
	./inc/mgc_num.o\
	./inc/scratch.o\
	./lib/sort.o\
	./lib/cantor.o\
	./lib/stack.o\
	./lib/mathematics.o\
	./lib/stringutils.o\
	./lib/fileio.o\
	./inc/falphabet.o\
	./inc/createalphabet.o\
	./inc/encodeprobvec.o\
	./inc/intsequence.o\
	convertprobvec.o

READPROBSEQOBJ   = ./inc/prob_vec.o\
	./inc/mprintf.o\
	./lib/memman.o\
	./lib/vtprogressbar.o\
	./inc/e_malloc.o\
	./inc/fio.o\
	./inc/matrix.o\
	./inc/mgc_num.o\
	./inc/scratch.o\
	./lib/sort.o\
	./lib/cantor.o\
	./lib/stack.o\
	./lib/mathematics.o\
	./lib/stringutils.o\
	./lib/fileio.o\
	./inc/falphabet.o\
	./inc/createalphabet.o\
	./inc/encodeprobvec.o\
	./inc/intsequence.o\
	readprobseq.o

BMPROBSEQOBJ   = ./lib/memman.o\
	 ./inc/prob_vec.o\
	./inc/mprintf.o\
	./lib/vtprogressbar.o\
	./inc/e_malloc.o\
	./inc/fio.o\
	./inc/matrix.o\
	./inc/mgc_num.o\
	./inc/scratch.o\
	./lib/sort.o\
	./lib/cantor.o\
	./lib/stack.o\
	./lib/mathematics.o\
	./lib/stringutils.o\
	./lib/boyermoore.o\
	./lib/fileio.o\
	./inc/falphabet.o\
	./inc/createalphabet.o\
	./inc/encodeprobvec.o\
	./inc/intsequence.o\
	bmprobseq.o


SUFPROBSEQOBJ   = ./vendor/gnuplot/gnuplot_i.o\
	./lib/memman.o\
	./lib/list.o\
	./inc/prob_vec.o\
	./inc/mprintf.o\
	./lib/vtprogressbar.o\
	./inc/e_malloc.o\
	./inc/fio.o\
	./inc/matrix.o\
	./inc/mgc_num.o\
	./inc/scratch.o\
	./lib/sort.o\
	./lib/cantor.o\
	./lib/stack.o\
	./lib/mathematics.o\
	./lib/stringutils.o\
	./lib/boyermoore.o\
	./lib/fileio.o\
	./lib/dpalign.o\
	./inc/falphabet.o\
	./inc/createalphabet.o\
	./inc/encodeprobvec.o\
	./inc/intsequence.o\
	./inc/sufarray.o\
	./inc/depictseqs.o\
	./inc/mm.o\
	./inc/sufmatch.o\
	./inc/multiseq.o\
	sufprobseq.o


ALPHABETTESTOBJ = alphabettest.o\
	./lib/stringutils.o\
	./lib/fileio.o\
	./lib/sort.o\
	./lib/stack.o\
	./lib/cantor.o\
	./lib/mathematics.o\
	./inc/createalphabet.o\
	./inc/falphabet.o\

WURSTIMBISSOBJ   = ./vendor/gnuplot/gnuplot_i.o\
	./lib/memman.o\
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
	wurstimbiss.o\

IMSUBSTOBJ    = ./inc/prob_vec.o\
	./inc/mprintf.o\
	./lib/memman.o\
	./lib/vtprogressbar.o\
	./inc/e_malloc.o\
	./inc/fio.o\
	./inc/matrix.o\
	./inc/mgc_num.o\
	./inc/scratch.o\
	./lib/sort.o\
	./lib/cantor.o\
	./lib/stack.o\
	./lib/mathematics.o\
	./lib/stringutils.o\
	./lib/fileio.o\
	./inc/falphabet.o\
	./inc/createalphabet.o\
	./inc/encodeprobvec.o\
	./inc/intsequence.o\
	./inc/imsubst.o\
	createsubmatrix.o



all: readprobseq.x bmprobseq.x probvectest.x convertprobvec.x sufprobseq.x wurstimbiss.x createsubmatrix.x


mathematicstest.x: ${MATHEMATICSTESTOBJ}
	gcc $(CFLAGS) ${MATHEMATICSTESTOBJ} -o $@

alphabettest.x: ${ALPHABETTESTOBJ}
	gcc $(CFLAGS) ${ALPHABETTESTOBJ} -o $@

probvectest.x: ${PROBVECTESTOBJ}
	gcc $(CFLAGS) ${PROBVECTESTOBJ} -o $@ $(LDFLAGS)

readprobseq.x: ${READPROBSEQOBJ}
	gcc $(CFLAGS) ${READPROBSEQOBJ} -o $@ $(LDFLAGS)


sufprobseq.x: ${SUFPROBSEQOBJ}
	gcc $(CFLAGS) ${SUFPROBSEQOBJ} -o $@ $(LDFLAGS)


bmprobseq.x: ${BMPROBSEQOBJ}
	gcc $(CFLAGS) ${BMPROBSEQOBJ} -o $@ $(LDFLAGS)



createsubmatrix.x: ${IMSUBSTOBJ}
	gcc $(CFLAGS) ${IMSUBSTOBJ} -o $@ $(LDFLAGS)


convertprobvec.x: ${CONVERTPROBVECOBJ}
	gcc $(CFLAGS) ${CONVERTPROBVECOBJ} -o $@ $(LDFLAGS)


wurstimbiss.x: ${WURSTIMBISSOBJ}
	gcc $(CFLAGS) ${WURSTIMBISSOBJ} -o $@ $(LDFLAGS)

clean: 
	rm ./inc/*.o
	rm ./lib/*.o
	rm *.o


