#include <assert.h>
#include <ctype.h>
#include <stdlib.h>
#include <zlib.h>
#include <stdio.h>
#include <getopt.h>
#include "sickle.h"

int single_qual_threshold = 20;
int single_length_threshold = 20;

static struct option single_long_options[] = {
  {"fastq-file", required_argument, 0, 'f'},
  {"output-file", required_argument, 0, 'o'},
  {"qual-type", required_argument, 0, 't'},
  {"qual-threshold", optional_argument, 0, 'q'},
  {"length-threshold", optional_argument, 0, 'l'},
  {GETOPT_HELP_OPTION_DECL},
  {GETOPT_VERSION_OPTION_DECL},
  {NULL, 0, NULL, 0}
};

void single_usage (int status) {
  
  fprintf (stderr, "\nUsage: %s se -f <fastq sequence file> -t <quality type> -o <trimmed fastq file>\n\
\n\
Options:\n\
-f, --fastq-file, Input fastq file (required)\n\
-t, --qual-type, Type of quality values (illumina, phred, sanger) (required)\n", PROGRAM_NAME);

  fprintf (stderr, "-o, --output-file, Output trimmed fastq file (required)\n\
-q, --qual-threshold, Threshold for trimming based on average quality in a window. Default 20.\n\
-l, --length-threshold, Threshold to keep a read based on length after trimming. Default 20.\n\
--help, display this help and exit\n\
--version, output version information and exit\n\n");

  exit (status);
}

int single_main (int argc, char *argv[]) {

	gzFile se=NULL;
	kseq_t *fqrec;
	int l;
	FILE *outfile=NULL;
	int debug=0;
	int optc;
	extern char *optarg;
	int qualtype=-1;
	int p1cut;

	while (1) {
		int option_index = 0;
		optc = getopt_long (argc, argv, "df:t:o:q:l:",single_long_options, &option_index);

		if (optc == -1) break;

		switch (optc) {
			if (single_long_options[option_index].flag != 0) break;

			case 'f':
				se = gzopen (optarg, "r");
				if (!se) {
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
				outfile = fopen (optarg, "w");
				if (!outfile) {
					fprintf (stderr, "Could not open output file '%s'.\n", optarg);
					return EXIT_FAILURE;
				}
				break;

			case 'q':
				single_qual_threshold = atoi (optarg);
				if (single_qual_threshold < 0) {
					fprintf (stderr, "Quality threshold must be >= 0\n");
					return EXIT_FAILURE;
				}
				break;

			case 'l':
				single_length_threshold = atoi (optarg);
				if (single_length_threshold < 0) {
					fprintf (stderr, "Length threshold must be >= 0\n");
					return EXIT_FAILURE;
				}
				break;

			case 'd':
				debug = 1;
				break;

			case_GETOPT_HELP_CHAR(single_usage)
			case_GETOPT_VERSION_CHAR(PROGRAM_NAME, VERSION, AUTHORS);

			case '?':
				single_usage (EXIT_FAILURE);
				break;

			default:
				single_usage (EXIT_FAILURE);
				break;
		}
	}


	if (!se || !outfile || qualtype == -1) {
		single_usage (EXIT_FAILURE);
	}


	fqrec = kseq_init (se);

	while ((l = kseq_read (fqrec)) >= 0) {

		p1cut = sliding_window (fqrec, qualtype, single_length_threshold, single_qual_threshold);

		if (p1cut >= 0) {
			fprintf (outfile, "@%s\n", fqrec->name.s);
			fprintf (outfile, "%.*s\n", p1cut, fqrec->seq.s);
			fprintf (outfile, "+%s\n", fqrec->name.s);
			fprintf (outfile, "%.*s\n", p1cut, fqrec->qual.s);
		}
	}

	kseq_destroy (fqrec);
	gzclose (se);
	fclose (outfile);
	return 0;
}
