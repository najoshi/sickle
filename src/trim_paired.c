#include <assert.h>
#include <ctype.h>
#include <stdlib.h>
#include <zlib.h>
#include <stdio.h>
#include "kseq.h"
#include <getopt.h>
#include "sickle.h"

KSEQ_INIT(gzFile, gzread)

int paired_qual_threshold = 20;
int paired_length_threshold = 20;

static struct option paired_long_options[] = {
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

void paired_usage (int status) {

	fprintf (stderr, "\nUsage: %s pe -f <paired-end fastq file 1> -r <paired-end fastq file 2> -t <quality type> -o <trimmed pe file 1> -p <trimmed pe file 2> -s <trimmed singles file>\n\
\n\
Options:\n\
-f, --pe-file1, Input paired-end fastq file 1 (required, must have same number of records as pe2)\n\
-r, --pe-file2, Input paired-end fastq file 2 (required, must have same number of records as pe1)\n\
-t, --qual-type, Type of quality values (illumina, phred, sanger) (required)\n\
-o, --output-pe1, Output trimmed fastq file 1 (required)\n", PROGRAM_NAME);

	fprintf (stderr, "-p, --output-pe2, Output trimmed fastq file 2 (required)\n\
-s, --output-single, Output trimmed singles fastq file (required)\n\
-q, --qual-threshold, Threshold for trimming based on average quality in a window. Default 20.\n\
-l, --length-threshold, Threshold to keep a read based on length after trimming. Default 20.\n\
--help, display this help and exit\n\
--version, output version information and exit\n\n");

	exit (status);
}

int paired_main (int argc, char *argv[]) {

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
	int window_size,i,j,window_start,window_total,p1flag,p2flag,p1cut,p2cut;

	while (1) {
		int option_index = 0;
		optc = getopt_long (argc, argv, "df:r:t:o:p:s:q:l:", paired_long_options, &option_index);

		if (optc == -1) break;

		switch (optc) {
			if (paired_long_options[option_index].flag != 0) break;

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
				if (!strcmp (optarg, "illumina")) qualtype = ILLUMINA;
				else if (!strcmp (optarg, "phred")) qualtype = PHRED;
				else if (!strcmp (optarg, "solexa")) qualtype = SOLEXA;
				else if (!strcmp (optarg, "sanger")) qualtype = SANGER;
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
				paired_qual_threshold = atoi (optarg);
				if (paired_qual_threshold < 0) {
					fprintf (stderr, "Quality threshold must be >= 0\n");
					return EXIT_FAILURE;
				}
				break;

			case 'l':
				paired_length_threshold = atoi (optarg);
				if (paired_length_threshold < 0) {
					fprintf (stderr, "Length threshold must be >= 0\n");
					return EXIT_FAILURE;
				}
				break;

			case 'd':
				debug = 1;
				break;

			case_GETOPT_HELP_CHAR(paired_usage);
			case_GETOPT_VERSION_CHAR(PROGRAM_NAME, VERSION, AUTHORS);

			case '?':
				paired_usage (EXIT_FAILURE);
				break;

			default:
				paired_usage (EXIT_FAILURE);
				break;
		}
	}


	if (!pe1 || !pe2 || !outfile1 || !outfile2 || !single || qualtype == -1) {
		paired_usage (EXIT_FAILURE);
	}


	fqrec1 = kseq_init (pe1);
	fqrec2 = kseq_init (pe2);

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
			if (((double)window_total / (double)window_size < paired_qual_threshold) || 
				(window_start+window_size > fqrec1->qual.l)) {

				/* at what point in the window does the quality dip below the threshold? */
				for (j=window_start; j<window_start+window_size; j++) {
					if (get_quality_num (fqrec1->qual.s[j], qualtype) < paired_qual_threshold) {
						p1cut = j;
						if (p1cut < paired_length_threshold) {p1flag = 0;}
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
			if (((double)window_total / (double)window_size < paired_qual_threshold) || 
				(window_start+window_size > fqrec2->qual.l)) {

				/* at what point in the window does the quality dip below the threshold? */
				for (j=window_start; j<window_start+window_size; j++) {
					if (get_quality_num (fqrec2->qual.s[j], qualtype) < paired_qual_threshold) {
						p2cut = j;
						if (p2cut < paired_length_threshold) {p2flag = 0;}
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
	return EXIT_SUCCESS;
}
