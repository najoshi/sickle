#include <assert.h>
#include <ctype.h>
#include <stdlib.h>
#include <zlib.h>
#include <stdio.h>
#include <getopt.h>
#include <pthread.h>
#include <string.h>
#include "sickle.h"

int single_qual_threshold = 20;
int single_length_threshold = 20;
int num_threads_single = 1;

static struct option single_long_options[] = {
  {"fastq-file", required_argument, 0, 'f'},
  {"output-file", required_argument, 0, 'o'},
  {"qual-type", required_argument, 0, 't'},
  {"qual-threshold", optional_argument, 0, 'q'},
  {"length-threshold", optional_argument, 0, 'l'},
  {"no-fiveprime", optional_argument, 0, 'x'},
  {"discard-n", optional_argument, 0, 'n'},
  {"quiet", optional_argument, 0, 'z'},
  {"num-threads", optional_argument, 0, 'u'},
  {GETOPT_HELP_OPTION_DECL},
  {GETOPT_VERSION_OPTION_DECL},
  {NULL, 0, NULL, 0}
};

void single_usage (int status) {
  
  fprintf (stderr, "\nUsage: %s se -f <fastq sequence file> -t <quality type> -o <trimmed fastq file>\n\
\n\
Options:\n\
-f, --fastq-file, Input fastq file (required)\n\
-t, --qual-type, Type of quality values (solexa (CASAVA < 1.3), illumina (CASAVA 1.3 to 1.7), sanger (which is CASAVA >= 1.8)) (required)\n", PROGRAM_NAME);

  fprintf (stderr, "-o, --output-file, Output trimmed fastq file (required)\n\
-q, --qual-threshold, Threshold for trimming based on average quality in a window. Default 20.\n\
-l, --length-threshold, Threshold to keep a read based on length after trimming. Default 20.\n\
-x, --no-fiveprime, Don't do five prime trimming.\n\
-n, --discard-n, Discard sequences with any Ns in them.\n\
-u, --num-threads, number of threads to use (default 1)\n");

fprintf (stderr, "--quiet, Don't print out any trimming information\n\
--help, display this help and exit\n\
--version, output version information and exit\n\n");

  exit (status);
}

int single_main (int argc, char *argv[]) {

	FILE *se=NULL;
	FILE *outfile=NULL;
	int debug=0;
	int optc;
	extern char *optarg;
	int qualtype=-1;
	char *outfn=NULL;
	char *infn=NULL;
	int kept=0;
	int discard=0;
	int quiet=0;
	int no_fiveprime=0;
	int discard_n=0;
    int i=0;
    pthread_t *threads=NULL;
    int waitnum = 0;

	while (1) {
		int option_index = 0;
		optc = getopt_long (argc, argv, "df:t:o:q:l:zxn",single_long_options, &option_index);

		if (optc == -1) break;

		switch (optc) {
			if (single_long_options[option_index].flag != 0) break;

			case 'f':
				infn = (char*) malloc (strlen (optarg) + 1);
				strcpy (infn, optarg);
				break;

			case 't':
				if (!strcmp (optarg, "illumina")) qualtype = ILLUMINA;
				else if (!strcmp (optarg, "solexa")) qualtype = SOLEXA;
				else if (!strcmp (optarg, "sanger")) qualtype = SANGER;
				else {
					fprintf (stderr, "Error: Quality type '%s' is not a valid type.\n", optarg);
					return EXIT_FAILURE;
				}
				break;

			case 'o':
				outfn = (char*) malloc (strlen (optarg) + 1);
				strcpy (outfn, optarg);
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

			case 'x':
				no_fiveprime = 1;
				break;

			case 'n':
				discard_n = 1;
				break;

            case 'u':
                num_threads_single = atoi (optarg);
                if (num_threads_single <= 0) {
                    fprintf (stderr, "Number of threads must be > 0\n");
                    return EXIT_FAILURE;
                }

			case 'z':
				quiet=1;
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


	if (qualtype == -1 || !infn || !outfn) {
		single_usage (EXIT_FAILURE);
	}

	if (!strcmp (infn, outfn)) {
		fprintf (stderr, "Error: Input file is same as output file.\n");
		return EXIT_FAILURE;
	}

	se = fopen (infn, "r");
	if (!se) {
		fprintf (stderr, "Could not open input file '%s'.\n", infn);
		return EXIT_FAILURE;
	}

	outfile = fopen (outfn, "w");
	if (!outfile) {
		fprintf (stderr, "Could not open output file '%s'.\n", outfn);
		return EXIT_FAILURE;
	}


    threads = malloc (sizeof(pthread_t) * num_threads_single);
    global_swd = malloc (sizeof(global_swd) * num_threads_single);

    for (i=0; i<num_threads_single; i++) {
        global_swd[i].fqrec_r1 = get_fq_record (se);

        global_swd[i].qualtype = qualtype;
        global_swd[i].length_threshold = single_length_threshold;
        global_swd[i].qual_threshold = single_qual_threshold;
        global_swd[i].no_fiveprime = no_fiveprime;
        global_swd[i].discard_n = discard_n;

        pthread_create (&threads[i], NULL, &sw_call, &i);
    }
	
	while (1) {

        pthread_join (threads[waitnum], NULL);

		/* if sequence quality and length pass filter then output record, else discard */
		if (global_swd[waitnum].three_prime_cut_r1 >= 0) {
			fprintf (outfile, "%s", global_swd[waitnum].fqrec_r1->h1);

			/* This print statement prints out the sequence string starting from the 5' cut */
			/* and then only prints out to the 3' cut, however, we need to adjust the 3' cut */
			/* by subtracting the 5' cut because the 3' cut was calculated on the original sequence */
			fprintf (outfile, "%.*s\n", global_swd[waitnum].three_prime_cut_r1 - global_swd[waitnum].five_prime_cut_r1, global_swd[waitnum].fqrec_r1->seq + global_swd[waitnum].five_prime_cut_r1);

			fprintf (outfile, "+\n");
			fprintf (outfile, "%.*s\n", global_swd[waitnum].three_prime_cut_r1 - global_swd[waitnum].five_prime_cut_r1, global_swd[waitnum].fqrec_r1->qual + global_swd[waitnum].five_prime_cut_r1);

			kept++;
		}

		else discard++;

        free (global_swd[waitnum].fqrec_r1);
        global_swd[waitnum].fqrec_r1 = get_fq_record (se);

        if (global_swd[waitnum].fqrec_r1 == NULL) {break;}
        pthread_create (&threads[waitnum], NULL, &sw_call, &waitnum);

        waitnum = (waitnum + 1) % num_threads_single;
	}

	if (!quiet) fprintf (stdout, "\nFastQ records kept: %d\nFastQ records discarded: %d\n\n", kept, discard);

	fclose (se);
	fclose (outfile);
	return 0;
}
