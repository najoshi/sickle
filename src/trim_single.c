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

int single_qual_threshold = 20;
int single_length_threshold = 20;

static struct option single_long_options[] = { { "fastq-file",
required_argument, 0, 'f' }, { "output-file", required_argument, 0, 'o' }, {
		"qual-type",
		required_argument, 0, 't' }, { "qual-threshold",
optional_argument, 0, 'q' }, { "length-threshold",
optional_argument, 0, 'l' }, { "prefix-length",
optional_argument, 0, 'c' }, { "no-fiveprime",
optional_argument, 0, 'x' }, { "discard-n", optional_argument, 0, 'n' }, {
		"quiet", optional_argument, 0, 'z' }, {
GETOPT_HELP_OPTION_DECL }, { GETOPT_VERSION_OPTION_DECL }, {
NULL, 0, NULL, 0 } };

void single_usage(int status) {

	fprintf(stderr,
			"\nUsage: %s se -f <fastq sequence file> -t <quality type> -o <trimmed fastq file>\n\
\n\
Options:\n\
-f, --fastq-file, Input fastq file (required)\n\
-t, --qual-type, Type of quality values (solexa (CASAVA < 1.3), illumina (CASAVA 1.3 to 1.7), sanger (which is CASAVA >= 1.8)) (required)\n",
			PROGRAM_NAME);

	fprintf(stderr,
			"-o, --output-file, Output trimmed fastq file (required)\n\
-q, --qual-threshold, Threshold for trimming based on average quality in a window. Default 20.\n\
-l, --length-threshold, Threshold to keep a read based on length after trimming. Default 20.\n\
-x, --no-fiveprime, Don't do five prime trimming.\n\
-n, --discard-n, Discard sequences with any Ns in them.\n");

	fprintf(stderr,
			"-c, --prefix-length, Prefix length for color-space sequence processing. Default 0\n\
-z, --quiet, Don't print out any trimming information\n\
--help, display this help and exit\n\
--version, output version information and exit\n\n");

	exit(status);
}

int single_main(int argc, char *argv[]) {

	gzFile se = NULL;
	kseq_t *fqrec;
	int l;
	FILE *outfile = NULL;
	int debug = 0;
	int optc;
	extern char *optarg;
	int qualtype = -1;
	cutsites* p1cut;
	char *outfn = NULL;
	int csprefix = 0;
	char *infn = NULL;
	int kept = 0;
	int discard = 0;
	int quiet = 0;
	int no_fiveprime = 0;
	int discard_n = 0;

	while (1) {
		int option_index = 0;
		optc = getopt_long(argc, argv, "df:t:o:q:l:c:zxn", single_long_options,
				&option_index);

		if (optc == -1)
			break;

		switch (optc) {
		if (single_long_options[option_index].flag != 0)
			break;

	case 'f':
		infn = (char*) malloc(strlen(optarg) + 1);
		strcpy(infn, optarg);
		break;
	case 'c':
		csprefix = atoi(optarg);
		if (csprefix < 0) {
			fprintf(stderr, "Prefix length must be >= 0\n");
			return EXIT_FAILURE;
		}
		break;

	case 't':
		if (!strcmp(optarg, "illumina"))
			qualtype = ILLUMINA;
		else if (!strcmp(optarg, "solexa"))
			qualtype = SOLEXA;
		else if (!strcmp(optarg, "sanger"))
			qualtype = SANGER;
		else {
			fprintf(stderr, "Error: Quality type '%s' is not a valid type.\n",
					optarg);
			return EXIT_FAILURE;
		}
		break;

	case 'o':
		outfn = (char*) malloc(strlen(optarg) + 1);
		strcpy(outfn, optarg);
		break;

	case 'q':
		single_qual_threshold = atoi(optarg);
		if (single_qual_threshold < 0) {
			fprintf(stderr, "Quality threshold must be >= 0\n");
			return EXIT_FAILURE;
		}
		break;

	case 'l':
		single_length_threshold = atoi(optarg);
		if (single_length_threshold < 0) {
			fprintf(stderr, "Length threshold must be >= 0\n");
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
		quiet = 1;
		break;

	case 'd':
		debug = 1;
		break;

	case_GETOPT_HELP_CHAR(single_usage)
	case_GETOPT_VERSION_CHAR(PROGRAM_NAME, VERSION, AUTHORS)
		;

	case '?':
		single_usage(EXIT_FAILURE);
		break;

	default:
		single_usage(EXIT_FAILURE);
		break;
		}
	}

	if (qualtype == -1 || !infn || !outfn) {
		single_usage(EXIT_FAILURE);
	}

	if (!strcmp(infn, outfn)) {
		fprintf(stderr, "Error: Input file is same as output file.\n");
		return EXIT_FAILURE;
	}

	se = gzopen(infn, "r");
	if (!se) {
		fprintf(stderr, "Could not open input file '%s'.\n", infn);
		return EXIT_FAILURE;
	}

	outfile = fopen(outfn, "w");
	if (!outfile) {
		fprintf(stderr, "Could not open output file '%s'.\n", outfn);
		return EXIT_FAILURE;
	}

	fqrec = kseq_init(se);

	while ((l = kseq_read(fqrec)) >= 0) {

		p1cut = sliding_window(fqrec, qualtype, single_length_threshold,
				single_qual_threshold, no_fiveprime, discard_n);

		if (debug)
			fprintf(stderr, "P1cut: \n%d, %d\n", p1cut->five_prime_cut,
					p1cut->three_prime_cut);
		/* if sequence quality and length pass filter then output record, else discard */
		if (p1cut->three_prime_cut >= 0) {
			fprintf(outfile, "@%s", fqrec->name.s);
			if (fqrec->comment.l)
				fprintf(outfile, " %s\n", fqrec->comment.s);
			else
				fprintf(outfile, "\n");

			/* This print statement prints out the sequence string starting from the 5' cut */
			/* and then only prints out to the 3' cut, however, we need to adjust the 3' cut */
			/* by subtracting the 5' cut because the 3' cut was calculated on the original sequence */
			if (csprefix == 0) {/*No color sequence prefix*/
				fprintf(outfile, "%.*s\n",
						p1cut->three_prime_cut - p1cut->five_prime_cut,
						fqrec->seq.s + p1cut->five_prime_cut);

				fprintf(outfile, "+\n");
				fprintf(outfile, "%.*s\n",
						p1cut->three_prime_cut - p1cut->five_prime_cut,
						fqrec->qual.s + p1cut->five_prime_cut);

			} else {
				fprintf(outfile, "%.*s\n", p1cut->three_prime_cut + csprefix,
						fqrec->seq.s);

				fprintf(outfile, "+\n");
				fprintf(outfile, "%.*s\n", p1cut->three_prime_cut,
						fqrec->qual.s);

			}

			kept++;
		}

		else {
			discard++;
			/*if ( debug ) fprintf (stderr, "\n%s\n%s\n", fqrec->seq.s, fqrec->qual.s);*/
		}

		free(p1cut);
	}

	if (!quiet)
		fprintf(stdout,
				"\nFastQ records kept: %d\nFastQ records discarded: %d\n\n",
				kept, discard);

	kseq_destroy(fqrec);
	gzclose(se);
	fclose(outfile);
	return 0;
}
