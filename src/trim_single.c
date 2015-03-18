#include <assert.h>
#include <ctype.h>
#include <stdlib.h>
#include <zlib.h>
#include <stdio.h>
#include <getopt.h>
#include "sickle.h"
#include "kseq.h"
#include "print_record.h"

__KS_GETC(gzread, BUFFER_SIZE)
__KS_GETUNTIL(gzread, BUFFER_SIZE)
__KSEQ_READ

int single_qual_threshold = 20;
int single_length_threshold = 20;

static struct option single_long_options[] = {
    {"polyA-trimming", no_argument, 0, 'a'},
    {"fastq-file", required_argument, 0, 'f'},
    {"output-file", required_argument, 0, 'o'},
    {"qual-type", required_argument, 0, 't'},
    {"qual-threshold", required_argument, 0, 'q'},
    {"length-threshold", required_argument, 0, 'l'},
    {"no-fiveprime", no_argument, 0, 'x'},
    {"discard-n", no_argument, 0, 'n'},
    {"gzip-output", no_argument, 0, 'g'},
    {"quiet", no_argument, 0, 'z'},
    {GETOPT_HELP_OPTION_DECL},
    {GETOPT_VERSION_OPTION_DECL},
    {NULL, 0, NULL, 0}
};

void single_usage(int status, char *msg) {

    fprintf(stderr, "\nUsage: %s se [options] -f <fastq sequence file> -t <quality type> -o <trimmed fastq file>\n"
"\n"
"Options:\n"
"-f, --fastq-file, Input fastq file (required)\n"
"-t, --qual-type, Type of quality values (solexa (CASAVA < 1.3), illumina (CASAVA 1.3 to 1.7), sanger (which is CASAVA >= 1.8)) (required)\n"
"-o, --output-file, Output trimmed fastq file (required)\n", PROGRAM_NAME);

    fprintf(stderr, "-q, --qual-threshold, Threshold for trimming based on average quality in a window. Default 20.\n"
"-l, --length-threshold, Threshold to keep a read based on length after trimming. Default 20.\n"
"-x, --no-fiveprime, Don't do five prime trimming.\n"
"-n, --trunc-n, Truncate sequences at position of first N.\n, --polyA-min, Minimum length of Poly-A or Poly-T tail to trim\n"
"-A, --polyA-min, Minimum length of Poly-A or Poly-T tail to trim. Default 10\n");

    fprintf(stderr, "-E, --polyA-error, Maximum amount of errors allowed in Poly-A or Poly-T tail in when trimming. Default 3\n"
"-g, --gzip-output, Output gzipped files.\n"
"--quiet, Don't print out any trimming information\n"
"--help, display this help and exit\n"
"--version, output version information and exit\n\n");

    if (msg) fprintf(stderr, "%s\n\n", msg);
    exit(status);
}

int single_main(int argc, char *argv[]) {

    gzFile se = NULL;
    kseq_t *fqrec;
    int l;
    FILE *outfile = NULL;
    gzFile outfile_gzip = NULL;
    int debug = 0;
    int optc;
    extern char *optarg;
    int qualtype = -1;
    cutsites *p1cut;
    char *outfn = NULL;
    char *infn = NULL;
    int kept = 0;
    int discard = 0;
    int quiet = 0;
    int no_fiveprime = 0;
    int trunc_n = 0;
    int gzip_output = 0;
    int total=0;
    int poly_trimming=0;
    int cut_poly=0;
    int r1_poly=0;
    int five_prime_removed=0;
    int three_prime_removed=0;
    int r1_check = 0;
    int polyA_error = 3;
    int polyA_min = 10;

    while (1) {
        int option_index = 0;
        optc = getopt_long(argc, argv, "adf:t:o:q:l:zxngA:E:", single_long_options, &option_index);

        if (optc == -1)
            break;

        switch (optc) {
            if (single_long_options[option_index].flag != 0)
                break;
	case 'a':
	    poly_trimming=1;
	    break;		
        case 'A':
            polyA_min = atoi(optarg);
            if (polyA_min < 0) {
                fprintf(stderr, "PolyA Min tail must be larger than 0");
            }
            break;
        case 'E':
            polyA_error = atoi(optarg);
            if (polyA_error < 0) {
                fprintf(stderr, "PolyA Min tail must be larger than 0");
            }
            break;

        case 'f':
            infn = (char *) malloc(strlen(optarg) + 1);
            strcpy(infn, optarg);
            break;

        case 't':
            if (!strcmp(optarg, "illumina"))
                qualtype = ILLUMINA;
            else if (!strcmp(optarg, "solexa"))
                qualtype = SOLEXA;
            else if (!strcmp(optarg, "sanger"))
                qualtype = SANGER;
            else {
                fprintf(stderr, "Error: Quality type '%s' is not a valid type.\n", optarg);
                return EXIT_FAILURE;
            }
            break;

        case 'o':
            outfn = (char *) malloc(strlen(optarg) + 1);
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
            trunc_n = 1;
            break;

        case 'g':
            gzip_output = 1;
            break;

        case 'z':
            quiet = 1;
            break;

        case 'd':
            debug = 1;
            break;

        case_GETOPT_HELP_CHAR(single_usage)
        case_GETOPT_VERSION_CHAR(PROGRAM_NAME, VERSION, AUTHORS);

        case '?':
            single_usage(EXIT_FAILURE, NULL);
            break;

        default:
            single_usage(EXIT_FAILURE, NULL);
            break;
        }
    }


    if (qualtype == -1 || !infn || !outfn) {
        single_usage(EXIT_FAILURE, "****Error: Must have quality type, input file, and output file.");
    }

    if (!strcmp(infn, outfn)) {
        fprintf(stderr, "****Error: Input file is same as output file.\n\n");
        return EXIT_FAILURE;
    }

    se = gzopen(infn, "r");

    if (!se) {
        fprintf(stderr, "****Error: Could not open input file '%s'.\n\n", infn);
        return EXIT_FAILURE;
    }

    if (!gzip_output) {
	if (strcmp(outfn, "stdout") != 0 ) {
 	       outfile = fopen(outfn, "w");
        } else {
		outfile = stdout;
	}

	if (!outfile) {
            fprintf(stderr, "****Error: Could not open output file '%s'.\n\n", outfn);
            return EXIT_FAILURE;
        }
    } else {
        outfile_gzip = gzopen(outfn, "w");
        if (!outfile_gzip) {
            fprintf(stderr, "****Error: Could not open output file '%s'.\n\n", outfn);
            return EXIT_FAILURE;
        }
    }


    fqrec = kseq_init(se);

    while ((l = kseq_read(fqrec)) >= 0) {

        p1cut = sliding_window(fqrec, qualtype, single_length_threshold, single_qual_threshold, no_fiveprime, trunc_n, debug);
        total++;

	r1_check = 0;

	if (poly_trimming) {
 		cut_poly =  fqrec->seq.l - trim_poly_a(fqrec->seq.s, polyA_min, polyA_error, fqrec->seq.l, 0) + 1;
                p1cut->three_prime_cut = ((p1cut->three_prime_cut > cut_poly) ? cut_poly : p1cut->three_prime_cut);
		if (cut_poly == p1cut->three_prime_cut) {
			r1_check = 1;
			r1_poly++;
		}
                cut_poly = trim_poly_t(fqrec->seq.s, polyA_min, polyA_error, 0, fqrec->seq.l);
                p1cut->five_prime_cut = ((p1cut->five_prime_cut < cut_poly) ? cut_poly : p1cut->five_prime_cut);
		if (cut_poly == p1cut->five_prime_cut && r1_check != 1) {
			r1_check = 1;
			r1_poly++;
		}
	
	}


	if (p1cut->three_prime_cut - p1cut->five_prime_cut <= single_length_threshold) {
                p1cut->three_prime_cut = -1;
                p1cut->five_prime_cut = -1;
        } 

	if (p1cut->three_prime_cut > 0) {
                three_prime_removed += (fqrec->seq.l - p1cut->three_prime_cut);
        }

        if (p1cut->five_prime_cut > 0) {
                five_prime_removed += p1cut->five_prime_cut;
        }
	



        if (debug) printf("P1cut: %d,%d\n", p1cut->five_prime_cut, p1cut->three_prime_cut);

        /* if sequence quality and length pass filter then output record, else discard */
        if (p1cut->three_prime_cut >= 0 && (p1cut->three_prime_cut - p1cut->five_prime_cut) >= single_length_threshold) {
            if (!gzip_output) {
                /* This print statement prints out the sequence string starting from the 5' cut */
                /* and then only prints out to the 3' cut, however, we need to adjust the 3' cut */
                /* by subtracting the 5' cut because the 3' cut was calculated on the original sequence */

                print_record (outfile, fqrec, p1cut);
            } else {
                print_record_gzip (outfile_gzip, fqrec, p1cut);
            }

            kept++;
        }

        else discard++;

        free(p1cut);
    }

    if (!quiet) {
	fprintf(stderr, "SE input file: %s\n\nTotal FastQ records: %d\nFastQ records kept: %d\nFastQ records discarded: %d\n\n", infn, total, kept, discard);
	if (poly_trimming) {
		fprintf(stderr, "Poly AT tails: %d\n", r1_poly);
	}
	fprintf(stderr, "Base pairs left removed: %d\nBase pairs right removed: %d\n", five_prime_removed, three_prime_removed);
    } 
    kseq_destroy(fqrec);
    gzclose(se);
    if (!gzip_output) fclose(outfile);
    else gzclose(outfile_gzip);

    return EXIT_SUCCESS;
}
