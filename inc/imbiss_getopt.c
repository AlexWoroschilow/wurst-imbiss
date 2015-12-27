#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <getopt.h>
#include <sys/types.h>
#include <sys/times.h>
#include <time.h>
#include <string.h>
#include <ncurses.h>
#include <memory.h>
#include <stddef.h>

const char * getopt_configfile(int argc, char** argv, const char * defaultfile) {
	int c;
	while ((c = getopt(argc, argv, "c:")) != -1) {
		switch (c) {
		case 'c':
			return (const char *) optarg;
		}
	}
	return defaultfile;
}
