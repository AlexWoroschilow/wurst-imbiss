#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include "memory.h"
#include "stringutils.h"
#include "fileio.h"
#include "intsequence.h"
#include "salami.h"
#include "massert.h"

/* ----------------------------- printSequence -------------------------------
 *    
 * prints a beautiful IntSequence to a buffer 
 * 
 */

char * printSequence(void *space, IntSequence *s, Uint cols) {
	Uint i, c, k, pos = 0, l = 0, width = 0;
	char *buf, *buf2;
	stringset_t *entries;

	entries = initStringset(space);
	for (i = 0; i < s->length; i++) {
		buf = ALLOCMEMORY(space, NULL, char, 32);
		buf = my_itoa(s->sequence[i], buf, 10);
		addString(space, entries, buf, strlen(buf));
		if (SETSTRLEN(entries, i)> width)
		width = SETSTRLEN(entries, i);
	}

	/*add spacer*/
	width++;
	c = cols / (width + 1);
	l = (s->descrlen + 2) + (s->namelen + 2) + (s->length * (width + 1)) + ((entries->noofstrings / c) * (5 + 1 + 1));

	buf = ALLOCMEMORY(space, NULL, char, l);
	memset(&buf[pos++], '>', 1);
	memmove(&buf[pos], s->url, s->urllen);
	pos += s->urllen;
	memset(&buf[pos++], '\n', 1);
	memset(&buf[pos++], '>', 1);
	memmove(&buf[pos], s->alphabetname, s->namelen);
	pos += s->namelen;
	memset(&buf[pos++], '\n', 1);
	for (i = 0; i < entries->noofstrings; i++) {
		if ((i % c) == 0) {
			memset(&buf[pos++], '\n', 1);
			buf2 = ALLOCMEMORY(space, NULL, char, 5);
			buf2 = my_itoa(i, buf2, 10);
			memset(&buf[pos], ' ', 5 - strlen(buf2));
			pos += 5 - strlen(buf2);
			memmove(&buf[pos], buf2, strlen(buf2));
			pos += strlen(buf2);
			memset(&buf[pos++], '\t', 1);
			FREEMEMORY(space, buf2);
		}
		k = (width - SETSTRLEN(entries, i));
		memset(&buf[pos], ' ', k);
		pos += k;

		memmove(&buf[pos], SETSTR(entries, i), SETSTRLEN(entries, i));
		pos += SETSTRLEN(entries, i);
	}

	buf[pos] = '\0';

	destructStringset(space, entries);
	return buf;
}

/* ----------------------------- printAlignment -------------------------------
 *    
 * prints a beautiful alignment to a buffer 
 * 
 */

char *
printAlignment(void *space, int *align, Uint size, IntSequence *a, IntSequence *b, Uint cols) {
	Uint i, c, k, pos1 = 0, pos2 = 0, pos3 = 0, l = 0, width = 0, m = 0, d = size;
	char *buf, *bufa, *bufb, *nobuf;
	stringset_t *first, *second;

	first = initStringset(space);
	for (i = 0; i < a->length; i++) {
		buf = ALLOCMEMORY(space, NULL, char, 32);
		buf = my_itoa(a->sequence[i], buf, 10);
		addString(space, first, buf, strlen(buf));
		if (SETSTRLEN(first, i)> width)
		width = SETSTRLEN(first, i);
	}

	second = initStringset(space);
	for (i = 0; i < b->length; i++) {
		buf = ALLOCMEMORY(space, NULL, char, 32);
		buf = my_itoa(b->sequence[i], buf, 10);
		addString(space, second, buf, strlen(buf));
		if (SETSTRLEN(second, i)> width)
		width = SETSTRLEN(second, i);
	}

	/*add spacer*/
	width++;
	c = cols / (width + 1);
	l = (a->descrlen + 2) + (a->namelen + 2) + (a->length * (width + 1)) + ((first->noofstrings / c) * (5 + 1 + 1));
	m = (b->length * (width + 1)) + ((second->noofstrings / c) * (5 + 1 + 1));

	bufa = ALLOCMEMORY(space, NULL, char, l);
	bufb = ALLOCMEMORY(space, NULL, char, m);
	buf = ALLOCMEMORY(space, NULL, char, (l + m) * 2);

	memset(bufa, 0, l);
	memset(bufb, 0, m);
	memset(buf, 0, (l + m) * 2);

	memset(&bufa[pos1++], '>', 1);
	memmove(&bufa[pos1], a->url, a->urllen);
	pos1 += a->urllen;
	memset(&bufa[pos1++], '\n', 1);
	memset(&bufa[pos1++], '>', 1);
	memmove(&bufa[pos1], a->alphabetname, a->namelen);
	pos1 += a->namelen;
	memset(&bufa[pos1++], '\n', 1);

	for (i = 0; i < first->noofstrings; i++) {
		if ((i % c) == 0) {
			memset(&bufa[pos1++], '\n', 1);
			memset(&bufb[pos2++], '\n', 1);

			memmove(&buf[pos3], bufa, pos1);
			pos3 += pos1;
			memmove(&buf[pos3], bufb, pos2);
			pos3 += pos2;
			memset(&buf[pos3++], '\n', 1);

			pos1 = 0;
			pos2 = 0;

			nobuf = ALLOCMEMORY(space, NULL, char, 5);
			nobuf = my_itoa(i, nobuf, 10);
			memset(&bufa[pos1], ' ', 5 - strlen(nobuf));
			pos1 += 5 - strlen(nobuf);
			memmove(&bufa[pos1], nobuf, strlen(nobuf));
			pos1 += strlen(nobuf);

			if (d > 1)
				nobuf = my_itoa(align[d - 1], nobuf, 10);
			else
				nobuf = my_itoa(align[1], nobuf, 10);

			memset(&bufb[pos2], ' ', 5 - strlen(nobuf));
			pos2 += 5 - strlen(nobuf);
			memmove(&bufb[pos2], nobuf, strlen(nobuf));
			pos2 += strlen(nobuf);

			memset(&bufa[pos1++], '\t', 1);
			memset(&bufb[pos2++], '\t', 1);

			FREEMEMORY(space, nobuf);
		}

		k = (width - SETSTRLEN(first, i));
		memset(&bufa[pos1], ' ', k);
		pos1 += k;
		memmove(&bufa[pos1], SETSTR(first, i), SETSTRLEN(first, i));
		pos1 += SETSTRLEN(first, i);

		if (d > 1 && align[d - 2] - 1 == i) {
			k = (width - SETSTRLEN(second, align[d - 1] - 1));
			memset(&bufb[pos2], ' ', k);
			pos2 += k;
			memmove(&bufb[pos2], SETSTR(second, align[d - 1] - 1), SETSTRLEN(second, align[d - 1] - 1));
			pos2 += SETSTRLEN(second, align[d - 1] - 1);
			d -= 2;
		} else {
			k = width - 1;
			memset(&bufb[pos2], ' ', k);
			pos2 += k;
			memset(&bufb[pos2], '-', 1);
			pos2++;
		}
	}

	memset(&bufa[pos1++], '\n', 1);
	memset(&bufb[pos2++], '\n', 1);

	memmove(&buf[pos3], bufa, pos1);
	pos3 += pos1;
	memmove(&buf[pos3], bufb, pos2);
	pos3 += pos2;
	memset(&buf[pos3++], '\n', 1);

	buf[pos3] = '\0';

	destructStringset(space, first);
	destructStringset(space, second);
	FREEMEMORY(space, bufa);
	FREEMEMORY(space, bufb);

	return buf;
}

/* ------------------------------ dumpSequence -------------------------------
 *  dumps the sequence to the screen.  
 * 
 */

void dumpSequence(IntSequence *s) {
	Uint i;
	printf("sequence:\n");
	for (i = 0; i < s->length; i++) {
		printf("%d", s->sequence[i]);
		if (i != (s->length - 1))
			printf("-");
	}
	printf("\n");
	printf("info:\n");
	for (i = 0; i < s->length; i++) {
		printf("%d", s->info[i]);
		if (i != (s->length - 1))
			printf("-");
	}
	printf("\n");
}

/*------------------------------- loadSequence -------------------------------
 *  loads a sequence from a file.   
 * 
 */

IntSequence* loadSequence(void *space, char *filename) {
	FILE *infile;
	Uint *sequence, *info;
	IntSequence *s;

	infile = fopen(filename, "r");
	if (infile == NULL) {
		fprintf( stderr, "couldn't open file '%s'; %s\n", filename, strerror(errno));
		exit(EXIT_FAILURE);
	}

	s = initSequence(space);
	massert(fread(s, sizeof(IntSequence), 1, infile), "Can not read sequence input file");

	s->description = ALLOCMEMORY(space, NULL, char, s->descrlen + 1);
	s->alphabetname = ALLOCMEMORY(space, NULL, char, s->namelen + 1);
	s->url = ALLOCMEMORY(space, NULL, char, s->urllen + 1);
	sequence = ALLOCMEMORY(space, NULL, Uint, s->length);
	info = ALLOCMEMORY(space, NULL, Uint, s->length);

	massert(fread(s->description, sizeof(char), s->descrlen + 1, infile), "Can not read sequence description");
	massert(fread(s->alphabetname, sizeof(char), s->namelen + 1, infile), "Can not read sequence alphabet name");
	massert(fread(s->url, sizeof(char), s->urllen + 1, infile), "Can not read sequence url");
	massert(fread(sequence, sizeof(Uint), s->length, infile), "Can not read sequence");
	massert(fread(info, sizeof(Uint), s->length, infile), "Can not read sequence info");

	s->sequence = sequence;
	s->info = info;

	if (fclose(infile) == EOF) /* close input file */
	{
		fprintf( stderr, "couldn't close file '%s'; %s\n", filename, strerror(errno));
		exit(EXIT_FAILURE);
	}

	return s;
}

/*------------------------------- saveSequence -------------------------------
 *
 *  saves the sequences to a file
 *
 */

void saveSequence(IntSequence *s, char *filename) {
	FILE *outfile;

	outfile = fopen(filename, "w");
	if (outfile == NULL) {
		fprintf( stderr, "couldn't open file '%s'; %s\n", filename, strerror(errno));
		exit(EXIT_FAILURE);
	}

	fwrite(s, sizeof(IntSequence), 1, outfile);
	fwrite(s->description, sizeof(char), s->descrlen + 1, outfile);
	fwrite(s->alphabetname, sizeof(char), s->namelen + 1, outfile);
	fwrite(s->url, sizeof(char), s->urllen + 1, outfile);
	fwrite(s->sequence, sizeof(Uint), s->length, outfile);
	fwrite(s->info, sizeof(Uint), s->length, outfile);

	if (fclose(outfile) == EOF) {
		fprintf( stderr, "couldn't close file '%s'; %s\n", filename, strerror(errno));
		exit(EXIT_FAILURE);
	}
	return;
}

/*---------------------------- createSequenceHash ----------------------------
 *    
 * creates a hash table (array) for IntSequences
 * 
 */

IntSequence**
createSequenceHash(void *space, Uint hashsize) {
	IntSequence** hashTable;
	hashTable = ALLOCMEMORY(space, NULL, IntSequence*, hashsize);
	memset(hashTable, 0, sizeof(IntSequence*) * hashsize);
	return hashTable;
}

/*----------------------------- findSequenceHash -----------------------------
 *    
 * find a slot in the sequence hash
 * 
 */

void findSequenceHash() {
	return;
}

/*---------------------------- lookupSequenceHash ----------------------------
 *    
 * lookup if a key is already present in a hash table
 * 
 */

void lookupSequenceHash() {
	return;
}

/*------------------------------- initSequence -------------------------------
 *    
 * allocates and initializes a sequence struct
 * 
 */

IntSequence*
initSequence(void *space) {
	IntSequence *s;

	s = ALLOCMEMORY(space, NULL, IntSequence, 1);
	s->description = NULL;
	s->descrlen = 0;
	s->alphabetname = NULL;
	s->namelen = 0;
	s->sequence = NULL;
	s->info = NULL;
	s->length = 0;
	s->urllen = 0;
	s->url = NULL;
	return s;
}

void destructSequence(void *space, IntSequence *sequence) {

	FREEMEMORY(space, sequence->sequence);
	FREEMEMORY(space, sequence->description);
	FREEMEMORY(space, sequence->alphabetname);
	FREEMEMORY(space, sequence->url);
	FREEMEMORY(space, sequence->info);
	FREEMEMORY(space, sequence);
	return;
}

IntSequence* sequence_init(void *space) {
	return initSequence(space);
}

IntSequence ** sequence_load_csv(void *imbiss, void *space, char* filename, char *delimeter, Uint *linecount,
		IntSequence* (*loader)(void *imbiss, void *space, char *filename)) {

	Uint i;
	stringset_t **fn = readcsv(space, filename, delimeter, linecount);
	IntSequence ** sequences = ALLOCMEMORY(space, NULL, IntSequence *, *linecount);
	for (i = 0; i < *linecount; i++) {
		const char * file = (const char *) SETSTR(fn[i], 0);
		sequences[i] = loader(imbiss, space, (char *) file);
	}

	for (i = 0; i < *linecount; i++) {
		destructStringset(space, fn[i]);
	}
	FREEMEMORY(space, fn);
	return sequences;
}

void sequence_dump(IntSequence *s) {
	sequence_dump_salami(s);
}

void sequence_dump_salami(IntSequence *s) {
	dumpSequence(s);
}

void sequence_dump_aacid(IntSequence *s) {
	Uint i;
	printf("sequence:\n");
	for (i = 0; i < s->length; i++) {
		printf("%c", s->sequence[i]);
		if (i != (s->length - 1))
			printf("-");
	}
	printf("\n");
	printf("info:\n");
	for (i = 0; i < s->length; i++) {
		printf("%d", s->info[i]);
		if (i != (s->length - 1))
			printf("-");
	}
	printf("\n");
}

char * sequence_print(void *space, IntSequence *s, Uint cols) {
	return printSequence(space, s, cols);
}

Uint * sequence_salami_to_uint(struct salami_sequence * sequence) {
	unsigned long i;
	Uint * sequence_uint = ALLOCMEMORY(space, NULL, Uint, sequence->length);
	for (i = 0; i < sequence->length; i++) {
		sequence_uint[i] = (Uint) sequence->sequence[i];
	}
	return sequence_uint;
}

char * sequence_code(char *url) {
	int i;
	Uint length = (strlen(url) - 1);
	for (i = length; i >= 0; i--) {
		if (url[i - 1] == '/') {
			break;
		}
	}
	return strdup((const char *) url + i);
}

IntSequence* sequence_aacid_load(void *imbiss, void *space, char *filename) {
	FILE *infile;

	massert(((infile = fopen(filename, "r"))!= NULL), "Couldn't open file");

	IntSequence *sequence = sequence_init(space);
	massert(fread(sequence, sizeof(*sequence), 1, infile), "Can not read sequence input file");

	sequence->description = ALLOCMEMORY(space, NULL, char, sequence->descrlen + 1);
	sequence->alphabetname = ALLOCMEMORY(space, NULL, char, sequence->namelen + 1);
	sequence->url = ALLOCMEMORY(space, NULL, char, sequence->urllen + 1);

	massert(fread(sequence->description, sizeof(char), sequence->descrlen + 1, infile), "Can not read sequence description");
	massert(fread(sequence->alphabetname, sizeof(char), sequence->namelen + 1, infile), "Can not read sequence alphabet name");
	massert(fread(sequence->url, sizeof(char), sequence->urllen + 1, infile), "Can not read sequence url");

	Uint *info = ALLOCMEMORY(space, NULL, Uint, sequence->length);
	massert(fread(info, sizeof(Uint), sequence->length, infile), "Can not read sequence file");

	struct salami_sequence * sequence_salami = salami_sequence_string(imbiss, sequence);
	massert((sequence_salami != NULL), "Salami sequence can not be null");

	sequence->length = sequence_salami->length;
	sequence->sequence = sequence_salami_to_uint(sequence_salami);
	sequence->info = info;

	massert((fclose(infile) != EOF), "Couldn't close file");

	return sequence;
}

IntSequence* sequence_salami_load(void *imbiss, void *space, char *filename) {
	return loadSequence(space, filename);
}

void sequence_destruct(void *space, IntSequence *sequence) {
	destructSequence(space, sequence);
}
