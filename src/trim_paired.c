#include <assert.h>
#include <ctype.h>
#include <stdlib.h>
#include <zlib.h>
#include <stdio.h>
#include <getopt.h>
#include "sickle.h"
#include "kseq.h"

__KS_GETC(gzread, BUFFER_SIZE)
__KS_GETUNTIL(gzread, BUFFER_SIZE)
__KSEQ_READ

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
	{"no-fiveprime", optional_argument, 0, 'x'},
	{"discard-n", optional_argument, 0, 'n'},
	{"quiet", optional_argument, 0, 'z'},
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
-t, --qual-type, Type of quality values (illumina (CASAVA 1.3 to 1.7), sanger (which is CASAVA >= 1.8)) (required)\n", PROGRAM_NAME);

	fprintf (stderr, "-o, --output-pe1, Output trimmed fastq file 1 (required)\n\
-p, --output-pe2, Output trimmed fastq file 2 (required)\n\
-s, --output-single, Output trimmed singles fastq file (required)\n\
-q, --qual-threshold, Threshold for trimming based on average quality in a window. Default 20.\n\
-l, --length-threshold, Threshold to keep a read based on length after trimming. Default 20.\n\
-x, --no-fiveprime, Don't do five prime trimming.\n\
-n, --discard-n, Discard sequences with any Ns in them.\n");


	fprintf (stderr, "--quiet, do not output trimming info\n\
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
	cutsites* p1cut;
	cutsites* p2cut;
	char *outfn1=NULL;
	char *outfn2=NULL;
	char *sfn=NULL;
	char *infn1=NULL;
	char *infn2=NULL;
	int kept_p=0;
	int discard_p=0;
	int kept_s1=0;
	int kept_s2=0;
	int discard_s1=0;
	int discard_s2=0;
	int quiet=0;
	int no_fiveprime=0;
	int discard_n=0;

	while (1) {
		int option_index = 0;
		optc = getopt_long (argc, argv, "df:r:t:o:p:s:q:l:xn", paired_long_options, &option_index);

		if (optc == -1) break;

		switch (optc) {
			if (paired_long_options[option_index].flag != 0) break;

			case 'f':
				infn1 = (char*) malloc (strlen (optarg) + 1);
				strcpy (infn1, optarg);
				break;

			case 'r':
				infn2 = (char*) malloc (strlen (optarg) + 1);
				strcpy (infn2, optarg);
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
				outfn1 = (char*) malloc (strlen (optarg) + 1);
				strcpy (outfn1, optarg);
				break;

			case 'p':
				outfn2 = (char*) malloc (strlen (optarg) + 1);
				strcpy (outfn2, optarg);
				break;

			case 's':
				sfn = (char*) malloc (strlen (optarg) + 1);
				strcpy (sfn, optarg);
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

			case 'x':
				no_fiveprime = 1;
				break;

			case 'n':
				discard_n = 1;
				break;

			case 'z':
				quiet=1;
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


	if (qualtype == -1 || !infn1 || !infn2 || !outfn1 || !outfn2 || !sfn) {
		paired_usage (EXIT_FAILURE);
	}

	if (!strcmp (infn1, infn2) || !strcmp (infn1, outfn1) || !strcmp (infn1, outfn2) ||
		!strcmp (infn1, sfn) || !strcmp (infn2, outfn1) || !strcmp (infn2, outfn2) ||
		!strcmp (infn2, sfn) || !strcmp (outfn1, outfn2) || !strcmp (outfn1, sfn) ||
		!strcmp (outfn2, sfn)) {

		fprintf (stderr, "Error: Duplicate input and/or output file names.\n");
		return EXIT_FAILURE;
	}

	pe1 = gzopen (infn1, "r");
	if (!pe1) {
		fprintf (stderr, "Could not open input file '%s'.\n", infn1);
		return EXIT_FAILURE;
	}

	pe2 = gzopen (infn2, "r");
	if (!pe2) {
		fprintf (stderr, "Could not open input file '%s'.\n", infn2);
		return EXIT_FAILURE;
	}

	outfile1 = fopen (outfn1, "w");
	if (!outfile1) {
		fprintf (stderr, "Could not open output file '%s'.\n", outfn1);
		return EXIT_FAILURE;
	}

	outfile2 = fopen (outfn2, "w");
	if (!outfile2) {
		fprintf (stderr, "Could not open output file '%s'.\n", outfn2);
		return EXIT_FAILURE;
	}

	single = fopen (sfn, "w");
	if (!single) {
		fprintf (stderr, "Could not open output file '%s'.\n", sfn);
		return EXIT_FAILURE;
	}


	fqrec1 = kseq_init (pe1);
	fqrec2 = kseq_init (pe2);

	while ((l1 = kseq_read (fqrec1)) >= 0) {

		l2 = kseq_read (fqrec2);
		if (l2 < 0) {
			fprintf (stderr, "Error: PE file 2 is shorter than PE file 1. Disregarding rest of PE file 1.\n");
			break;
		}

		p1cut = sliding_window (fqrec1, qualtype, paired_length_threshold, paired_qual_threshold, no_fiveprime, discard_n);
		p2cut = sliding_window (fqrec2, qualtype, paired_length_threshold, paired_qual_threshold, no_fiveprime, discard_n);

		/* The sequence and quality print statements below print out the sequence string starting from the 5' cut */
		/* and then only print out to the 3' cut, however, we need to adjust the 3' cut */
		/* by subtracting the 5' cut because the 3' cut was calculated on the original sequence */

		/* if both sequences passed quality and length filters, then output both records */
		if (p1cut->three_prime_cut >= 0 && p2cut->three_prime_cut >= 0) {
			fprintf (outfile1, "@%s", fqrec1->name.s);
			if (fqrec1->comment.l) fprintf (outfile1, " %s\n", fqrec1->comment.s);
			else fprintf (outfile1, "\n");
			fprintf (outfile1, "%.*s\n", p1cut->three_prime_cut - p1cut->five_prime_cut, fqrec1->seq.s + p1cut->five_prime_cut);

			fprintf (outfile1, "+%s", fqrec1->name.s);
			if (fqrec1->comment.l) fprintf (outfile1, " %s\n", fqrec1->comment.s);
			else fprintf (outfile1, "\n");
			fprintf (outfile1, "%.*s\n", p1cut->three_prime_cut - p1cut->five_prime_cut, fqrec1->qual.s + p1cut->five_prime_cut);


			fprintf (outfile2, "@%s", fqrec2->name.s);
			if (fqrec2->comment.l) fprintf (outfile2, " %s\n", fqrec2->comment.s);
			else fprintf (outfile2, "\n");
			fprintf (outfile2, "%.*s\n", p2cut->three_prime_cut - p2cut->five_prime_cut, fqrec2->seq.s + p2cut->five_prime_cut);

			fprintf (outfile2, "+%s", fqrec2->name.s);
			if (fqrec2->comment.l) fprintf (outfile2, " %s\n", fqrec2->comment.s);
			else fprintf (outfile2, "\n");
			fprintf (outfile2, "%.*s\n", p2cut->three_prime_cut - p2cut->five_prime_cut, fqrec2->qual.s + p2cut->five_prime_cut);

			kept_p += 2;
		}

		/* if only one sequence passed filter, then put its record in singles and discard the other */
		else if (p1cut->three_prime_cut >= 0 && p2cut->three_prime_cut < 0) {
			fprintf (single, "@%s", fqrec1->name.s);
			if (fqrec1->comment.l) fprintf (single, " %s\n", fqrec1->comment.s);
			else fprintf (single, "\n");
			fprintf (single, "%.*s\n", p1cut->three_prime_cut - p1cut->five_prime_cut, fqrec1->seq.s + p1cut->five_prime_cut);

			fprintf (single, "+%s", fqrec1->name.s);
			if (fqrec1->comment.l) fprintf (single, " %s\n", fqrec1->comment.s);
			else fprintf (single, "\n");
			fprintf (single, "%.*s\n", p1cut->three_prime_cut - p1cut->five_prime_cut, fqrec1->qual.s + p1cut->five_prime_cut);

			kept_s1++;
			discard_s2++;
		}

		else if (p1cut->three_prime_cut < 0 && p2cut->three_prime_cut >= 0) {
			fprintf (single, "@%s", fqrec2->name.s);
			if (fqrec2->comment.l) fprintf (single, " %s\n", fqrec2->comment.s);
			else fprintf (single, "\n");
			fprintf (single, "%.*s\n", p2cut->three_prime_cut - p2cut->five_prime_cut, fqrec2->seq.s + p2cut->five_prime_cut);

			fprintf (single, "+%s", fqrec2->name.s);
			if (fqrec2->comment.l) fprintf (single, " %s\n", fqrec2->comment.s);
			else fprintf (single, "\n");
			fprintf (single, "%.*s\n", p2cut->three_prime_cut - p2cut->five_prime_cut, fqrec2->qual.s + p2cut->five_prime_cut);

			kept_s2++;
			discard_s1++;
		}

		else discard_p += 2;

		free(p1cut);
		free(p2cut);
	}

	if (l1 < 0) {
		l2 = kseq_read (fqrec2);
		if (l2 >= 0) {
			fprintf (stderr, "Error: PE file 1 is shorter than PE file 2. Disregarding rest of PE file 2.\n");
		}
	}

	if (!quiet) fprintf (stdout, "\nFastQ paired records kept: %d (%d pairs)\nFastQ single records kept: %d (from PE1: %d, from PE2: %d)\nFastQ paired records discarded: %d (%d pairs)\nFastQ single records discarded: %d (from PE1: %d, from PE2: %d)\n\n", kept_p, (kept_p/2), (kept_s1+kept_s2), kept_s1, kept_s2, discard_p, (discard_p/2), (discard_s1+discard_s2), discard_s1, discard_s2);

	kseq_destroy (fqrec1);
	kseq_destroy (fqrec2);
	gzclose (pe1);
	gzclose (pe2);
	fclose (outfile1);
	fclose (outfile2);
	fclose (single);
	return EXIT_SUCCESS;
}
