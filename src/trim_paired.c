#include <assert.h>
#include <ctype.h>
#include <stdlib.h>
#include <limits.h>
#include <zlib.h>
#include <stdio.h>
#include "kseq.h"
#include <getopt.h>

KSEQ_INIT(gzFile, gzread)

#ifndef PROGRAM_NAME
#define PROGRAM_NAME "trim_paired"
#endif

#ifndef AUTHORS
#define AUTHORS "Nikhil Joshi, UC Davis Bioinformatics Core\n"
#endif

#ifndef VERSION
#define VERSION 0.0
#endif

#define ILLUMINA_TYPE 1
#define PHRED_TYPE 2
#define SANGER_TYPE 3


/* Options drawn from GNU's coreutils/src/system.h */
/* These options are defined so as to avoid conflicting with option
values used by commands */
enum {
  GETOPT_HELP_CHAR = (CHAR_MIN - 2),
  GETOPT_VERSION_CHAR = (CHAR_MIN - 3)
};
#define GETOPT_HELP_OPTION_DECL \
"help", no_argument, NULL, GETOPT_HELP_CHAR
#define GETOPT_VERSION_OPTION_DECL \
"version", no_argument, NULL, GETOPT_VERSION_CHAR
#define case_GETOPT_HELP_CHAR \
case GETOPT_HELP_CHAR: \
usage(EXIT_SUCCESS); \
break;
#define case_GETOPT_VERSION_CHAR(Program_name, Version, Authors) \
case GETOPT_VERSION_CHAR: \
fprintf(stdout, "%s version %0.3f\nCopyright (c) 2011 The Regents " \
"of University of California, Davis Campus.\n" \
"%s is free software and comes with ABSOLUTELY NO WARRANTY.\n"\
"Distributed under the MIT License.\n\nWritten by %s\n", \
Program_name, Version, Program_name, Authors); \
exit(EXIT_SUCCESS); \
break;
/* end code drawn from system.h */


int qual_threshold = 20;
int length_threshold = 20;

static struct option long_options[] = {
	{"pe-file1", required_argument, 0, 'f'},
	{"pe-file2", required_argument, 0, 'r'},
	{"qual-type", required_argument, 0, 't'},
	{"output-pe1", required_argument, 0, 'o'},
	{"output-pe2", required_argument, 0, 'p'},
	{"output-single", required_argument, 0, 's'},
	{"qual-threshold", optional_argument, 0, 'q'},
	{"length-threshold", optional_argument, 0, 'l'},
	{GETOPT_HELP_OPTION_DECL},
	{GETOPT_VERSION_OPTION_DECL},
	{NULL, 0, NULL, 0}
};

void usage (int status) {

	fprintf (stdout, "\nUsage: %s -f <paired-end fastq file 1> -r <paired-end fastq file 2> -t <quality type> -o <trimmed pe file 1> -p <trimmed pe file 2> -s <trimmed singles file>\n\
\n\
Options:\n\
-f, --pe-file1, Input paired-end fastq file 1 (required, must have same number of records as pe2)\n\
-r, --pe-file2, Input paired-end fastq file 2 (required, must have same number of records as pe1)\n\
-t, --qual-type, Type of quality values (illumina, phred, sanger) (required)\n\
-o, --output-pe1, Output trimmed fastq file 1 (required)\n\
-p, --output-pe2, Output trimmed fastq file 2 (required)\n\
-s, --output-single, Output trimmed singles fastq file (required)\n\
-q, --qual-threshold, Threshold for trimming based on average quality in a window. Default 20.\n\
-l, --length-threshold, Threshold to keep a read based on length after trimming. Default 20.\n\
--help, display this help and exit\n\
--version, output version information and exit\n\n", PROGRAM_NAME);

	exit (status);
}

int get_quality_num (char qualchar, int qualtype) {
	if (qualtype == ILLUMINA_TYPE) return ((int) qualchar - 64);
	else if (qualtype == PHRED_TYPE) return ((int) qualchar);
	else if (qualtype == SANGER_TYPE) return ((int) qualchar - 33);
}

int main (int argc, char *argv[]) {

	gzFile pe1=NULL;
	gzFile pe2=NULL;
	kseq_t *fqrec1;
	kseq_t *fqrec2;
	int l1,l2;
	FILE *outfile1=NULL;
	FILE *outfile2=NULL;
	FILE *single=NULL;
	int debug=0;
	int optc;
	extern char *optarg;
	int qualtype=-1;

	while (1) {
		int option_index = 0;
		optc = getopt_long (argc, argv, "df:r:t:o:p:s:q:l:", long_options, &option_index);

		if (optc == -1) break;

		switch (optc) {
			if (long_options[option_index].flag != 0) break;

			case 'f':
				pe1 = gzopen (optarg, "r");
				if (!pe1) {
					fprintf (stderr, "Could not open fastq file '%s'.\n", optarg);
					return EXIT_FAILURE;
				}
				break;

			case 'r':
				pe2 = gzopen (optarg, "r");
				if (!pe2) {
					fprintf (stderr, "Could not open fastq file '%s'.\n", optarg);
					return EXIT_FAILURE;
				}
				break;

			case 't':
				if (!strcmp (optarg, "illumina")) qualtype = ILLUMINA_TYPE;
				else if (!strcmp (optarg, "phred")) qualtype = PHRED_TYPE;
				else if (!strcmp (optarg, "sanger")) qualtype = SANGER_TYPE;
				else {
					fprintf (stderr, "Error: Quality type '%s' is not a valid type.\n", optarg);
					return EXIT_FAILURE;
				}
				break;

			case 'o':
				outfile1 = fopen (optarg, "w");
				if (!outfile1) {
					fprintf (stderr, "Could not open output file '%s'.\n", optarg);
					return EXIT_FAILURE;
				}
				break;

			case 'p':
				outfile2 = fopen (optarg, "w");
				if (!outfile2) {
					fprintf (stderr, "Could not open output file '%s'.\n", optarg);
					return EXIT_FAILURE;
				}
				break;

			case 's':
				single = fopen (optarg, "w");
				if (!single) {
					fprintf (stderr, "Could not open output file '%s'.\n", optarg);
					return EXIT_FAILURE;
				}
				break;

			case 'q':
				qual_threshold = atoi (optarg);
				if (qual_threshold < 0) {
					fprintf (stderr, "Quality threshold must be >= 0\n");
					return EXIT_FAILURE;
				}
				break;

			case 'l':
				length_threshold = atoi (optarg);
				if (length_threshold < 0) {
					fprintf (stderr, "Length threshold must be >= 0\n");
					return EXIT_FAILURE;
				}
				break;

			case 'd':
				debug = 1;
				break;

			case_GETOPT_HELP_CHAR;
			case_GETOPT_VERSION_CHAR(PROGRAM_NAME, VERSION, AUTHORS);

			case '?':
				usage (EXIT_FAILURE);
				break;

			default:
				usage (EXIT_FAILURE);
				break;
		}
	}


	if (!pe1 || !pe2 || !outfile1 || !outfile2 || !single || qualtype == -1) {
		usage (EXIT_FAILURE);
	}


	fqrec1 = kseq_init (pe1);
	fqrec2 = kseq_init (pe2);

	int window_size,i,j,window_start,window_total,p1flag,p2flag,p1cut,p2cut;

	while ((l1 = kseq_read (fqrec1)) >= 0) {

		l2 = kseq_read (fqrec2);
		if (l2 < 0) {
			fprintf (stderr, "Error: PE file 2 is shorter than PE file 1. Disregarding rest of PE file 1.\n");
			break;
		}

		window_size = (int) (0.1 * fqrec1->seq.l);
		window_start=0;
		window_total=0;
		p1flag=1;
		p1cut = fqrec1->seq.l;

		/* if the seq length is less then 10bp, */
		/* then make the window size the length of the seq */
		if (window_size == 0) window_size = fqrec1->seq.l;

		for (i=0; i<window_size; i++) {
			window_total += get_quality_num (fqrec1->qual.s[i], qualtype);
		}
if (debug) printf ("window total: %d, window_size: %d\n", window_total, window_size);

		for (i=0; i<fqrec1->qual.l; i++) {

if (debug) printf ("window total / window size: %f\n", (double)window_total / (double)window_size); 

			/* if the average quality in the window is less than the threshold */
			/* or if the window is the last window in the read */
			if (((double)window_total / (double)window_size < qual_threshold) || 
				(window_start+window_size > fqrec1->qual.l)) {

				/* at what point in the window does the quality dip below the threshold? */
				for (j=window_start; j<window_start+window_size; j++) {
					if (get_quality_num (fqrec1->qual.s[j], qualtype) < qual_threshold) {
						p1cut = j;
						if (p1cut < length_threshold) {p1flag = 0;}
						break;
					}
				}

				break;
			}

			/* instead of sliding the window, subtract the first qual and add the next qual */
			window_total -= get_quality_num (fqrec1->qual.s[window_start], qualtype);
			window_total += get_quality_num (fqrec1->qual.s[window_start+window_size], qualtype);
			window_start++;

if (debug) printf ("window total: %d\n", window_total);
		}


		window_size = (int) (0.1 * fqrec2->seq.l);
		window_start=0;
		window_total=0;
		p2flag=1;
		p2cut = fqrec2->seq.l;

		/* if the seq length is less then 10bp, */
		/* then make the window size the length of the seq */
		if (window_size == 0) window_size = fqrec2->seq.l;

		for (i=0; i<window_size; i++) {
			window_total += get_quality_num (fqrec2->qual.s[i], qualtype);
		}
if (debug) printf ("window total: %d, window_size: %d\n", window_total, window_size);

		for (i=0; i<fqrec2->qual.l; i++) {

if (debug) printf ("window total / window size: %f\n", (double)window_total / (double)window_size); 

			/* if the average quality in the window is less than the threshold */
			/* or if the window is the last window in the read */
			if (((double)window_total / (double)window_size < qual_threshold) || 
				(window_start+window_size > fqrec2->qual.l)) {

				/* at what point in the window does the quality dip below the threshold? */
				for (j=window_start; j<window_start+window_size; j++) {
					if (get_quality_num (fqrec2->qual.s[j], qualtype) < qual_threshold) {
						p2cut = j;
						if (p2cut < length_threshold) {p2flag = 0;}
						break;
					}
				}

				break;
			}

			/* instead of sliding the window, subtract the first qual and add the next qual */
			window_total -= get_quality_num (fqrec2->qual.s[window_start], qualtype);
			window_total += get_quality_num (fqrec2->qual.s[window_start+window_size], qualtype);
			window_start++;

if (debug) printf ("window total: %d\n", window_total);
		}

		if (p1flag == 1 && p2flag == 1) {
			fprintf (outfile1, "@%s\n", fqrec1->name.s);
			fprintf (outfile1, "%.*s\n", p1cut, fqrec1->seq.s);
			fprintf (outfile1, "+%s\n", fqrec1->name.s);
			fprintf (outfile1, "%.*s\n", p1cut, fqrec1->qual.s);

			fprintf (outfile2, "@%s\n", fqrec2->name.s);
			fprintf (outfile2, "%.*s\n", p2cut, fqrec2->seq.s);
			fprintf (outfile2, "+%s\n", fqrec2->name.s);
			fprintf (outfile2, "%.*s\n", p2cut, fqrec2->qual.s);
		}

		else if (p1flag == 1 && p2flag == 0) {
			fprintf (single, "@%s\n", fqrec1->name.s);
			fprintf (single, "%.*s\n", p1cut, fqrec1->seq.s);
			fprintf (single, "+%s\n", fqrec1->name.s);
			fprintf (single, "%.*s\n", p1cut, fqrec1->qual.s);
		}

		else if (p1flag == 0 && p2flag == 1) {
			fprintf (single, "@%s\n", fqrec2->name.s);
			fprintf (single, "%.*s\n", p2cut, fqrec2->seq.s);
			fprintf (single, "+%s\n", fqrec2->name.s);
			fprintf (single, "%.*s\n", p2cut, fqrec2->qual.s);
		}
	}

	if (l1 < 0) {
		l2 = kseq_read (fqrec2);
		if (l2 >= 0) {
			fprintf (stderr, "Error: PE file 1 is shorter than PE file 2. Disregarding rest of PE file 2.\n");
		}
	}

	kseq_destroy (fqrec1);
	kseq_destroy (fqrec2);
	gzclose (pe1);
	gzclose (pe2);
	fclose (outfile1);
	fclose (outfile2);
	fclose (single);
	return 0;
}
