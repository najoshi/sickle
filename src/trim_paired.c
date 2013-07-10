#include <assert.h>
#include <ctype.h>
#include <stdlib.h>
#include <zlib.h>
#include <stdio.h>
#include <getopt.h>
#include <unistd.h>
#include <pthread.h>
#include <string.h>
#include "sickle.h"

int paired_qual_threshold = 20;
int paired_length_threshold = 20;
int num_threads_paired = 1;
pthread_t *threads=NULL;

swdata *global_swd;
int waitnum;

static struct option paired_long_options[] = {
	{"qual-type", required_argument, 0, 't'},
	{"pe-file1", required_argument, 0, 'f'},
	{"pe-file2", required_argument, 0, 'r'},
	{"output-pe1", required_argument, 0, 'o'},
	{"output-pe2", required_argument, 0, 'p'},
	{"output-single", required_argument, 0, 's'},
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

void paired_usage (int status) {

  fprintf(stderr, "\nIf your have separate files for forward and reverse reads...\n");
	fprintf (stderr, "\nUsage: %s pe -f <paired-end fastq file 1> -r <paired-end fastq file 2> -t <quality type> -o <trimmed pe file 1> -p <trimmed pe file 2> -s <trimmed singles file>\n",PROGRAM_NAME);
  fprintf(stderr, "\nIf your have one file with interleaved forward and reverse reads...\n");
	fprintf (stderr, "\nUsage: %s pe -c <combined input file> -t <quality type> -m <combined trimmed output> -s <trimmed singles file>\n",PROGRAM_NAME);
  fprintf (stderr, "\n Options:\n\
-t, --qual-type, Type of quality values (solexa (CASAVA < 1.3), illumina (CASAVA 1.3 to 1.7), sanger (which is CASAVA >= 1.8)) (required)\n\
-f, --pe-file1, Input paired-end fastq file 1 (optional, must have same number of records as pe2)\n\
-r, --pe-file2, Input paired-end fastq file 2 (optional, must have same number of records as pe1)\n\
-c, --pe-combo, Combined input paired-end fastq (optional)\n");

	fprintf (stderr, "-o, --output-pe1, Output trimmed fastq file 1 (optional)\n\
-p, --output-pe2, Output trimmed fastq file 2 (optional)\n\
-m, --output-combo, Output combined paired-end fastq file (optional)\n\
-s, --output-single, Output trimmed singles fastq file (required)\n");
	fprintf (stderr, "-q, --qual-threshold, Threshold for trimming based on average quality in a window. Default 20.\n\
-l, --length-threshold, Threshold to keep a read based on length after trimming. Default 20.\n\
-x, --no-fiveprime, Don't do five prime trimming.\n\
-n, --discard-n, Discard sequences with any Ns in them.\n");


	fprintf (stderr, "-u, --num-threads, number of threads to use (default 1)\n\
--quiet, do not output trimming info\n\
--help, display this help and exit\n\
--version, output version information and exit\n\n");

	exit (status);
}

int paired_main (int argc, char *argv[]) {

	FILE *pe1=NULL;        /* forward input file handle */
	FILE *pe2=NULL;        /* reverse input file handle */
	gzFile pec=NULL;        /* combined input file handle */
	FILE *outfile1=NULL;    /* forward output file handle */
	FILE *outfile2=NULL;    /* reverse output file handle */
	FILE *combo=NULL;       /* combined output file handle */
	FILE *single=NULL;      /* single output file handle */
	int debug=0;
	int optc;
	extern char *optarg;
	int qualtype=-1;
	char *outfn1=NULL;      /* forward file out name */
	char *outfn2=NULL;      /* reverse file out name */
	char *outfnc=NULL;      /* combined file out name */
	char *sfn=NULL;         /* single/combined file out name */
	char *infn1=NULL;       /* forward input filename */
	char *infn2=NULL;       /* reverse input filename */
	char *infnc=NULL;       /* combined input filename */
	int kept_p=0;
	int discard_p=0;
	int kept_s1=0;
	int kept_s2=0;
	int discard_s1=0;
	int discard_s2=0;
	int quiet=0;
	int no_fiveprime=0;
	int discard_n=0;
    int i=0;
    int *tids;
    waitnum = 0;


	while (1) {
		int option_index = 0;
		optc = getopt_long (argc, argv, "df:r:t:o:p:s:q:l:u:xn", paired_long_options, &option_index);

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

			case 'u':
				num_threads_paired = atoi (optarg);
				if (num_threads_paired <= 0) {
					fprintf (stderr, "Number of threads must be > 0\n");
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

  /* required: qualtype */
	if (qualtype == -1) {
        fprintf (stderr, "Error: Quality type (-t) is required.\n");
        paired_usage (EXIT_FAILURE);
	}

  /* required: singlton filename */
	if (!sfn) {
        fprintf (stderr, "Error: Single read filename (-s) is required.\n");
        paired_usage (EXIT_FAILURE);
	}

  /* make sure minimum input filenames are specified*/
  if (!infn1 && !infn2 && !infnc) {
    fprintf (stderr, "Error: Must have either -f and -r OR -c arguments.\n");
    paired_usage(EXIT_FAILURE);
  }

  if (infnc) { /* using combined input file */
    fprintf(stderr,"COMBO\n");
    /* check for duplicate file names */
    if (!strcmp (infnc, outfnc) || !strcmp (infnc, sfn) || !strcmp (outfnc, sfn)) {
      fprintf (stderr, "Error: Duplicate filename between combo input, combo output, and/or single output file names.\n");
      return EXIT_FAILURE;
    }
    /* make sure minimum input filenames are specified*/
    if (!infn1 && !infn2 && !infnc) {
      fprintf (stderr, "Error: Must have either -f and -r OR -c arguments.\n");
      paired_usage(EXIT_FAILURE);
    }
    if (infn1 || infn2) {
      fprintf (stderr, "Cannot have both the -c option for combined input and either -f or -r\n");
      return EXIT_FAILURE;
    }
    if (outfn1 || outfn2) {
      fprintf (stderr, "Cannot have -o or -p when using combined output combo (-m) option\n");
      return EXIT_FAILURE;
    }

    /* get combined output file */
    combo = fopen (outfnc, "w");
    if (!combo) {
      fprintf (stderr, "Could not open combo output file '%s'.\n", outfnc);
      return EXIT_FAILURE;
    }

    pec = gzopen(infnc,"r");
    if (!pec) {
      fprintf (stderr, "Could not open combined input file '%s'.\n", infnc);
      return EXIT_FAILURE;
    }
  } else { /* using forward and reverse input files */
    if ( (infn1 && !infn2) || (infn2 && !infn1)) {
      fprintf (stderr, "Error: Must have both -f and -r arguments if using forward and reverse files.\n");
      return EXIT_FAILURE;
    }

    if (!outfn1 && !outfn2) {
      fprintf (stderr, "Error: Must have -o and -p if using separate forward and reverse files.\n");
      return EXIT_FAILURE;
    }

    if ( (outfn1 && !outfn2) || (outfn2 && !outfn1)) {
      fprintf (stderr, "Error: Must have both -o and -p arguments if using forward and reverse files.\n");
      return EXIT_FAILURE;
    }

    if (!strcmp (infn1, infn2) || !strcmp (infn1, outfn1) || !strcmp (infn1, outfn2) ||
      !strcmp (infn1, sfn) || !strcmp (infn2, outfn1) || !strcmp (infn2, outfn2) ||
      !strcmp (infn2, sfn) || !strcmp (outfn1, outfn2) || !strcmp (outfn1, sfn) ||
      !strcmp (outfn2, sfn)) {

      fprintf (stderr, "Error: Duplicate input and/or output file names.\n");
      return EXIT_FAILURE;
    }
 
    pe1 = fopen (infn1, "r");
    if (!pe1) {
      fprintf (stderr, "Could not open input file '%s'.\n", infn1);
      return EXIT_FAILURE;
    }

    pe2 = fopen (infn2, "r");
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
  }

  /* get singles output file handle*/
  single = fopen (sfn, "w");
  if (!single) {
    fprintf (stderr, "Could not open single output file '%s'.\n", sfn);
    return EXIT_FAILURE;
  }

/* printf ("Got here\n"); */

  threads = malloc (sizeof(pthread_t) * num_threads_paired);
  global_swd = malloc (sizeof(swdata) * num_threads_paired);
  tids = malloc (sizeof(int) * num_threads_paired);

/* printf ("Got here2\n"); */

  for (i=0; i<num_threads_paired; i++) {

/* printf ("Got here3\n"); */

      global_swd[i].fqrec_r1 = get_fq_record (pe1);
      global_swd[i].fqrec_r2 = get_fq_record (pe2);

/* printf ("Got here3.5: %s\n", global_swd[i].fqrec_r1->seq); */

      global_swd[i].qualtype = qualtype;
      global_swd[i].length_threshold = paired_length_threshold;
      global_swd[i].qual_threshold = paired_qual_threshold;
      global_swd[i].no_fiveprime = no_fiveprime;
      global_swd[i].discard_n = discard_n;

/* printf ("Got here4: %d\n", global_swd[i].length_threshold); */
  }

  for (i=0; i<num_threads_paired; i++) {
/* printf ("about to call create: %d\n", i); */
      tids[i] = i;
      pthread_create (&threads[i], NULL, &sw_call, &tids[i]);
  }

/* printf ("Got here5\n"); */

  while (1) {

/* printf ("ABOUT TO WAIT on thread %d\n", waitnum); */

    pthread_join (threads[waitnum], NULL);

    /* The sequence and quality print statements below print out the sequence string starting from the 5' cut */
    /* and then only print out to the 3' cut, however, we need to adjust the 3' cut */
    /* by subtracting the 5' cut because the 3' cut was calculated on the original sequence */

    /* if both sequences passed quality and length filters, then output both records */
    if (global_swd[waitnum].three_prime_cut_r1 >= 0 && global_swd[waitnum].three_prime_cut_r2 >= 0) {

        fprintf (outfile1, "%s", global_swd[waitnum].fqrec_r1->h1);
        fprintf (outfile1, "%.*s\n", global_swd[waitnum].three_prime_cut_r1 - global_swd[waitnum].five_prime_cut_r1, global_swd[waitnum].fqrec_r1->seq + global_swd[waitnum].five_prime_cut_r1);
        fprintf (outfile1, "+\n");
        fprintf (outfile1, "%.*s\n", global_swd[waitnum].three_prime_cut_r1 - global_swd[waitnum].five_prime_cut_r1, global_swd[waitnum].fqrec_r1->qual + global_swd[waitnum].five_prime_cut_r1);

        fprintf (outfile2, "%s", global_swd[waitnum].fqrec_r2->h1);
        fprintf (outfile2, "%.*s\n", global_swd[waitnum].three_prime_cut_r2 - global_swd[waitnum].five_prime_cut_r2, global_swd[waitnum].fqrec_r2->seq + global_swd[waitnum].five_prime_cut_r2);
        fprintf (outfile2, "+\n");
        fprintf (outfile2, "%.*s\n", global_swd[waitnum].three_prime_cut_r2 - global_swd[waitnum].five_prime_cut_r2, global_swd[waitnum].fqrec_r2->qual + global_swd[waitnum].five_prime_cut_r2);
/*
        fprintf (outfile1, "@%s", fqrec1->name.s);
        if (fqrec1->comment.l) fprintf (outfile1, " %s\n", fqrec1->comment.s);
        else fprintf (outfile1, "\n");
        fprintf (outfile1, "%.*s\n", p1cut->three_prime_cut - p1cut->five_prime_cut, fqrec1->seq.s + p1cut->five_prime_cut);
        fprintf (outfile1, "+\n");
        fprintf (outfile1, "%.*s\n", p1cut->three_prime_cut - p1cut->five_prime_cut, fqrec1->qual.s + p1cut->five_prime_cut);
        fprintf (outfile2, "@%s", fqrec2->name.s);
        if (fqrec2->comment.l) fprintf (outfile2, " %s\n", fqrec2->comment.s);
        else fprintf (outfile2, "\n");
        fprintf (outfile2, "%.*s\n", p2cut->three_prime_cut - p2cut->five_prime_cut, fqrec2->seq.s + p2cut->five_prime_cut);
        fprintf (outfile2, "+\n");
        fprintf (outfile2, "%.*s\n", p2cut->three_prime_cut - p2cut->five_prime_cut, fqrec2->qual.s + p2cut->five_prime_cut);
*/
        kept_p += 2;
    }

    /* if only one sequence passed filter, then put its record in singles and discard the other */
    else if (global_swd[waitnum].three_prime_cut_r1 >= 0 && global_swd[waitnum].three_prime_cut_r2 < 0) {

        fprintf (single, "%s", global_swd[waitnum].fqrec_r1->h1);
        fprintf (single, "%.*s\n", global_swd[waitnum].three_prime_cut_r1 - global_swd[waitnum].five_prime_cut_r1, global_swd[waitnum].fqrec_r1->seq + global_swd[waitnum].five_prime_cut_r1);
        fprintf (single, "+\n");
        fprintf (single, "%.*s\n", global_swd[waitnum].three_prime_cut_r1 - global_swd[waitnum].five_prime_cut_r1, global_swd[waitnum].fqrec_r1->qual + global_swd[waitnum].five_prime_cut_r1);
    }

    else if (global_swd[waitnum].three_prime_cut_r1 < 0 && global_swd[waitnum].three_prime_cut_r2 >= 0) {

        fprintf (single, "%s", global_swd[waitnum].fqrec_r2->h1);
        fprintf (single, "%.*s\n", global_swd[waitnum].three_prime_cut_r2 - global_swd[waitnum].five_prime_cut_r2, global_swd[waitnum].fqrec_r2->seq + global_swd[waitnum].five_prime_cut_r2);
        fprintf (single, "+\n");
        fprintf (single, "%.*s\n", global_swd[waitnum].three_prime_cut_r2 - global_swd[waitnum].five_prime_cut_r2, global_swd[waitnum].fqrec_r2->qual + global_swd[waitnum].five_prime_cut_r2);
    }

    else discard_p += 2;

    free (global_swd[waitnum].fqrec_r1);
    free (global_swd[waitnum].fqrec_r2);

    global_swd[waitnum].fqrec_r1 = get_fq_record (pe1);
    global_swd[waitnum].fqrec_r2 = get_fq_record (pe2);

    if (global_swd[waitnum].fqrec_r1 == NULL) {break;}

/* printf ("about to call create again with waitnum %d\n", waitnum); */
    tids[waitnum] = waitnum;
    pthread_create (&threads[waitnum], NULL, &sw_call, &tids[waitnum]);

    waitnum = (waitnum + 1) % num_threads_paired;

  } /* end of while(1) */

	if (!quiet) {
    fprintf (stdout, "\nFastQ paired records kept: %d (%d pairs)\n",kept_p, (kept_p/2));
    if (pec) {
      fprintf (stdout, "FastQ single records kept: %d\n", (kept_s1+kept_s2));
    } else {
      fprintf (stdout, "FastQ single records kept: %d (from PE1: %d, from PE2: %d)\n", (kept_s1+kept_s2), kept_s1, kept_s2);
    }

    fprintf (stdout, "FastQ paired records discarded: %d (%d pairs)\n", discard_p, (discard_p/2));

    if (pec) {
      fprintf (stdout, "FastQ single records discarded: %d\n\n", (discard_s1+discard_s2));
    } else {
      fprintf (stdout, "FastQ single records discarded: %d (from PE1: %d, from PE2: %d)\n\n", (discard_s1+discard_s2), discard_s1, discard_s2);
    }
  }

    fclose (single);
    fclose (pe1);
    fclose (pe2);
    fclose (outfile1);
    fclose (outfile2);
	return EXIT_SUCCESS;
} /* end of paired_main() */
