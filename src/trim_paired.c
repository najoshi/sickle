#include <assert.h>
#include <ctype.h>
#include <stdlib.h>
#include <zlib.h>
#include <stdio.h>
#include <getopt.h>
#include <unistd.h>
#include "sickle.h"
#include "kseq.h"
#include "print_record.h"

__KS_GETC(gzread, BUFFER_SIZE)
__KS_GETUNTIL(gzread, BUFFER_SIZE)
__KSEQ_READ

int paired_qual_threshold = 20;
int paired_length_threshold = 20;

static struct option paired_long_options[] = {
    {"qual-type", required_argument, 0, 't'},
    {"pe-file1", required_argument, 0, 'f'},
    {"pe-file2", required_argument, 0, 'r'},
    {"pe-combo", required_argument, 0, 'c'},
    {"output-pe1", required_argument, 0, 'o'},
    {"output-pe2", required_argument, 0, 'p'},
    {"output-single", required_argument, 0, 's'},
    {"output-combo", required_argument, 0, 'm'},
    {"qual-threshold", required_argument, 0, 'q'},
    {"length-threshold", required_argument, 0, 'l'},
    {"no-fiveprime", no_argument, 0, 'x'},
    {"truncate-n", no_argument, 0, 'n'},
    {"gzip-output", no_argument, 0, 'g'},
    {"output-combo-all", required_argument, 0, 'M'},
    {"quiet", no_argument, 0, 'z'},
    {GETOPT_HELP_OPTION_DECL},
    {GETOPT_VERSION_OPTION_DECL},
    {NULL, 0, NULL, 0}
};

void paired_usage (int status, char *msg) {

    fprintf(stderr, "\nIf you have separate files for forward and reverse reads:\n");
    fprintf(stderr, "Usage: %s pe [options] -f <paired-end forward fastq file> -r <paired-end reverse fastq file> -t <quality type> -o <trimmed PE forward file> -p <trimmed PE reverse file> -s <trimmed singles file>\n\n", PROGRAM_NAME);
    fprintf(stderr, "If you have one file with interleaved forward and reverse reads:\n");
    fprintf(stderr, "Usage: %s pe [options] -c <interleaved input file> -t <quality type> -m <interleaved trimmed paired-end output> -s <trimmed singles file>\n\n\
If you have one file with interleaved reads as input and you want ONLY one interleaved file as output:\n\
Usage: %s pe [options] -c <interleaved input file> -t <quality type> -M <interleaved trimmed output>\n\n", PROGRAM_NAME, PROGRAM_NAME);
    fprintf(stderr, "Options:\n\
Paired-end separated reads\n\
--------------------------\n\
-f, --pe-file1, Input paired-end forward fastq file (Input files must have same number of records)\n\
-r, --pe-file2, Input paired-end reverse fastq file\n\
-o, --output-pe1, Output trimmed forward fastq file\n\
-p, --output-pe2, Output trimmed reverse fastq file. Must use -s option.\n\n\
Paired-end interleaved reads\n\
----------------------------\n");
    fprintf(stderr,"-c, --pe-combo, Combined (interleaved) input paired-end fastq\n\
-m, --output-combo, Output combined (interleaved) paired-end fastq file. Must use -s option.\n\
-M, --output-combo-all, Output combined (interleaved) paired-end fastq file with any discarded read written to output file as a single N. Cannot be used with the -s option.\n\n\
Global options\n\
--------------\n\
-t, --qual-type, Type of quality values (solexa (CASAVA < 1.3), illumina (CASAVA 1.3 to 1.7), sanger (which is CASAVA >= 1.8)) (required)\n");
    fprintf(stderr, "-s, --output-single, Output trimmed singles fastq file\n\
-q, --qual-threshold, Threshold for trimming based on average quality in a window. Default 20.\n\
-l, --length-threshold, Threshold to keep a read based on length after trimming. Default 20.\n\
-x, --no-fiveprime, Don't do five prime trimming.\n\
-n, --truncate-n, Truncate sequences at position of first N.\n");


    fprintf(stderr, "-g, --gzip-output, Output gzipped files.\n--quiet, do not output trimming info\n\
--help, display this help and exit\n\
--version, output version information and exit\n\n");

    if (msg) fprintf(stderr, "%s\n\n", msg);
    exit(status);
}


int paired_main(int argc, char *argv[]) {

    gzFile pe1 = NULL;          /* forward input file handle */
    gzFile pe2 = NULL;          /* reverse input file handle */
    gzFile pec = NULL;          /* combined input file handle */
    kseq_t *fqrec1 = NULL;
    kseq_t *fqrec2 = NULL;
    int l1, l2;
    FILE *outfile1 = NULL;      /* forward output file handle */
    FILE *outfile2 = NULL;      /* reverse output file handle */
    FILE *combo = NULL;         /* combined output file handle */
    FILE *single = NULL;        /* single output file handle */
    gzFile outfile1_gzip = NULL;
    gzFile outfile2_gzip = NULL;
    gzFile combo_gzip = NULL;
    gzFile single_gzip = NULL;
    int debug = 0;
    int optc;
    extern char *optarg;
    int qualtype = -1;
    cutsites *p1cut;
    cutsites *p2cut;
    char *outfn1 = NULL;        /* forward file out name */
    char *outfn2 = NULL;        /* reverse file out name */
    char *outfnc = NULL;        /* combined file out name */
    char *sfn = NULL;           /* single/combined file out name */
    char *infn1 = NULL;         /* forward input filename */
    char *infn2 = NULL;         /* reverse input filename */
    char *infnc = NULL;         /* combined input filename */
    int kept_p = 0;
    int discard_p = 0;
    int kept_s1 = 0;
    int kept_s2 = 0;
    int discard_s1 = 0;
    int discard_s2 = 0;
    int quiet = 0;
    int no_fiveprime = 0;
    int trunc_n = 0;
    int gzip_output = 0;
    int combo_all=0;
    int combo_s=0;
    int total=0;

    while (1) {
        int option_index = 0;
        optc = getopt_long(argc, argv, "df:r:c:t:o:p:m:M:s:q:l:xng", paired_long_options, &option_index);

        if (optc == -1)
            break;

        switch (optc) {
            if (paired_long_options[option_index].flag != 0)
                break;

        case 'f':
            infn1 = (char *) malloc(strlen(optarg) + 1);
            strcpy(infn1, optarg);
            break;

        case 'r':
            infn2 = (char *) malloc(strlen(optarg) + 1);
            strcpy(infn2, optarg);
            break;

        case 'c':
            infnc = (char *) malloc(strlen(optarg) + 1);
            strcpy(infnc, optarg);
            break;

        case 't':
            if (!strcmp(optarg, "illumina")) qualtype = ILLUMINA;
            else if (!strcmp(optarg, "solexa")) qualtype = SOLEXA;
            else if (!strcmp(optarg, "sanger")) qualtype = SANGER;
            else {
                fprintf(stderr, "Error: Quality type '%s' is not a valid type.\n", optarg);
                return EXIT_FAILURE;
            }
            break;

        case 'o':
            outfn1 = (char *) malloc(strlen(optarg) + 1);
            strcpy(outfn1, optarg);
            break;

        case 'p':
            outfn2 = (char *) malloc(strlen(optarg) + 1);
            strcpy(outfn2, optarg);
            break;

        case 'm':
            outfnc = (char *) malloc(strlen(optarg) + 1);
            strcpy(outfnc, optarg);
            combo_s = 1;
            break;

        case 'M':
            outfnc = (char *) malloc(strlen(optarg) + 1);
            strcpy(outfnc, optarg);
            combo_all = 1;
            break;

        case 's':
            sfn = (char *) malloc(strlen(optarg) + 1);
            strcpy(sfn, optarg);
            break;

        case 'q':
            paired_qual_threshold = atoi(optarg);
            if (paired_qual_threshold < 0) {
                fprintf(stderr, "Quality threshold must be >= 0\n");
                return EXIT_FAILURE;
            }
            break;

        case 'l':
            paired_length_threshold = atoi(optarg);
            if (paired_length_threshold < 0) {
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

        case_GETOPT_HELP_CHAR(paired_usage);
        case_GETOPT_VERSION_CHAR(PROGRAM_NAME, VERSION, AUTHORS);

        case '?':
            paired_usage(EXIT_FAILURE, NULL);
            break;

        default:
            paired_usage(EXIT_FAILURE, NULL);
            break;
        }
    }

    /* required: qualtype */
    if (qualtype == -1) {
        paired_usage(EXIT_FAILURE, "****Error: Quality type is required.");
    }

    /* make sure minimum input filenames are specified */
    if (!infn1 && !infnc) {
        paired_usage(EXIT_FAILURE, "****Error: Must have either -f OR -c argument.");
    }

    if (infnc) {      /* using combined input file */

        if (infn1 || infn2 || outfn1 || outfn2) {
            paired_usage(EXIT_FAILURE, "****Error: Cannot have -f, -r, -o, or -p options with -c.");
        }

        if ((combo_all && combo_s) || (!combo_all && !combo_s)) {
            paired_usage(EXIT_FAILURE, "****Error: Must have only one of either -m or -M options with -c.");
        }

        if ((combo_s && !sfn) || (combo_all && sfn)) {
            paired_usage(EXIT_FAILURE, "****Error: -m option must have -s option, and -M option cannot have -s option.");
        }

        /* check for duplicate file names */
        if (!strcmp(infnc, outfnc) || (combo_s && (!strcmp(infnc, sfn) || !strcmp(outfnc, sfn)))) {
            fprintf(stderr, "****Error: Duplicate filename between combo input, combo output, and/or single output file names.\n\n");
            return EXIT_FAILURE;
        }

        /* get combined output file */
        if (!gzip_output) {
            combo = fopen(outfnc, "w");
            if (!combo) {
                fprintf(stderr, "****Error: Could not open combo output file '%s'.\n\n", outfnc);
                return EXIT_FAILURE;
            }
        } else {
            combo_gzip = gzopen(outfnc, "w");
            if (!combo_gzip) {
                fprintf(stderr, "****Error: Could not open combo output file '%s'.\n\n", outfnc);
                return EXIT_FAILURE;
            }
        }

        pec = gzopen(infnc, "r");
        if (!pec) {
            fprintf(stderr, "****Error: Could not open combined input file '%s'.\n\n", infnc);
            return EXIT_FAILURE;
        }

    } else {     /* using forward and reverse input files */

        if (infn1 && (!infn2 || !outfn1 || !outfn2 || !sfn)) {
            paired_usage(EXIT_FAILURE, "****Error: Using the -f option means you must have the -r, -o, -p, and -s options.");
        }

        if (infn1 && (infnc || combo_all || combo_s)) {
            paired_usage(EXIT_FAILURE, "****Error: The -f option cannot be used in combination with -c, -m, or -M.");
        }

        if (!strcmp(infn1, infn2) || !strcmp(infn1, outfn1) || !strcmp(infn1, outfn2) ||
            !strcmp(infn1, sfn) || !strcmp(infn2, outfn1) || !strcmp(infn2, outfn2) || 
            !strcmp(infn2, sfn) || !strcmp(outfn1, outfn2) || !strcmp(outfn1, sfn) || !strcmp(outfn2, sfn)) {

            fprintf(stderr, "****Error: Duplicate input and/or output file names.\n\n");
            return EXIT_FAILURE;
        }

        pe1 = gzopen(infn1, "r");
        if (!pe1) {
            fprintf(stderr, "****Error: Could not open input file '%s'.\n\n", infn1);
            return EXIT_FAILURE;
        }

        pe2 = gzopen(infn2, "r");
        if (!pe2) {
            fprintf(stderr, "****Error: Could not open input file '%s'.\n\n", infn2);
            return EXIT_FAILURE;
        }

        if (!gzip_output) {
            outfile1 = fopen(outfn1, "w");
            if (!outfile1) {
                fprintf(stderr, "****Error: Could not open output file '%s'.\n\n", outfn1);
                return EXIT_FAILURE;
            }

            outfile2 = fopen(outfn2, "w");
            if (!outfile2) {
                fprintf(stderr, "****Error: Could not open output file '%s'.\n\n", outfn2);
                return EXIT_FAILURE;
            }
        } else {
            outfile1_gzip = gzopen(outfn1, "w");
            if (!outfile1_gzip) {
                fprintf(stderr, "****Error: Could not open output file '%s'.\n\n", outfn1);
                return EXIT_FAILURE;
            }

            outfile2_gzip = gzopen(outfn2, "w");
            if (!outfile2_gzip) {
                fprintf(stderr, "****Error: Could not open output file '%s'.\n\n", outfn2);
                return EXIT_FAILURE;
            }

        }
    }

    /* get singles output file handle */
    if (sfn && !combo_all) {
        if (!gzip_output) {
            single = fopen(sfn, "w");
            if (!single) {
                fprintf(stderr, "****Error: Could not open single output file '%s'.\n\n", sfn);
                return EXIT_FAILURE;
            }
        } else {
            single_gzip = gzopen(sfn, "w");
            if (!single_gzip) {
                fprintf(stderr, "****Error: Could not open single output file '%s'.\n\n", sfn);
                return EXIT_FAILURE;
            }
        }
    }

    if (pec) {
        fqrec1 = kseq_init(pec);
        fqrec2 = (kseq_t *) malloc(sizeof(kseq_t));
        fqrec2->f = fqrec1->f;
    } else {
        fqrec1 = kseq_init(pe1);
        fqrec2 = kseq_init(pe2);
    }

    while ((l1 = kseq_read(fqrec1)) >= 0) {

        l2 = kseq_read(fqrec2);
        if (l2 < 0) {
            fprintf(stderr, "Warning: PE file 2 is shorter than PE file 1. Disregarding rest of PE file 1.\n");
            break;
        }

        p1cut = sliding_window(fqrec1, qualtype, paired_length_threshold, paired_qual_threshold, no_fiveprime, trunc_n, debug);
        p2cut = sliding_window(fqrec2, qualtype, paired_length_threshold, paired_qual_threshold, no_fiveprime, trunc_n, debug);
        total += 2;

        if (debug) printf("p1cut: %d,%d\n", p1cut->five_prime_cut, p1cut->three_prime_cut);
        if (debug) printf("p2cut: %d,%d\n", p2cut->five_prime_cut, p2cut->three_prime_cut);

        /* The sequence and quality print statements below print out the sequence string starting from the 5' cut */
        /* and then only print out to the 3' cut, however, we need to adjust the 3' cut */
        /* by subtracting the 5' cut because the 3' cut was calculated on the original sequence */

        /* if both sequences passed quality and length filters, then output both records */
        if (p1cut->three_prime_cut >= 0 && p2cut->three_prime_cut >= 0) {
            if (!gzip_output) {
                if (pec) {
                    print_record (combo, fqrec1, p1cut);
                    print_record (combo, fqrec2, p2cut);
                } else {
                    print_record (outfile1, fqrec1, p1cut);
                    print_record (outfile2, fqrec2, p2cut);
                }
            } else {
                if (pec) {
                    print_record_gzip (combo_gzip, fqrec1, p1cut);
                    print_record_gzip (combo_gzip, fqrec2, p2cut);
                } else {
                    print_record_gzip (outfile1_gzip, fqrec1, p1cut);
                    print_record_gzip (outfile2_gzip, fqrec2, p2cut);
                }
            }

            kept_p += 2;
        }

        /* if only one sequence passed filter, then put its record in singles and discard the other */
        /* or put an "N" record in if that option was chosen. */
        else if (p1cut->three_prime_cut >= 0 && p2cut->three_prime_cut < 0) {
            if (!gzip_output) {
                if (combo_all) {
                    print_record (combo, fqrec1, p1cut);
                    print_record_N (combo, fqrec2, qualtype);
                } else {
                    print_record (single, fqrec1, p1cut);
                }
            } else {
                if (combo_all) {
                    print_record_gzip (combo_gzip, fqrec1, p1cut);
                    print_record_N_gzip (combo_gzip, fqrec2, qualtype);
                } else {
                    print_record_gzip (single_gzip, fqrec1, p1cut);
                }
            }

            kept_s1++;
            discard_s2++;
        }

        else if (p1cut->three_prime_cut < 0 && p2cut->three_prime_cut >= 0) {
            if (!gzip_output) {
                if (combo_all) {
                    print_record_N (combo, fqrec1, qualtype);
                    print_record (combo, fqrec2, p2cut);
                } else {
                    print_record (single, fqrec2, p2cut);
                }
            } else {
                if (combo_all) {
                    print_record_N_gzip (combo_gzip, fqrec1, qualtype);
                    print_record_gzip (combo_gzip, fqrec2, p2cut);
                } else {
                    print_record_gzip (single_gzip, fqrec2, p2cut);
                }
            }

            kept_s2++;
            discard_s1++;

        } else {

            /* If both records are to be discarded, but the -M option */
            /* is being used, then output two "N" records */
            if (combo_all) {
                if (!gzip_output) {
                    print_record_N (combo, fqrec1, qualtype);
                    print_record_N (combo, fqrec2, qualtype);
                } else {
                    print_record_N_gzip (combo_gzip, fqrec1, qualtype);
                    print_record_N_gzip (combo_gzip, fqrec2, qualtype);
                }
            }

            discard_p += 2;
        }

        free(p1cut);
        free(p2cut);
    }             /* end of while ((l1 = kseq_read (fqrec1)) >= 0) */

    if (l1 < 0) {
        l2 = kseq_read(fqrec2);
        if (l2 >= 0) {
            fprintf(stderr, "Warning: PE file 1 is shorter than PE file 2. Disregarding rest of PE file 2.\n");
        }
    }

    if (!quiet) {
        if (infn1 && infn2) fprintf(stdout, "\nPE forward file: %s\nPE reverse file: %s\n", infn1, infn2);
        if (infnc) fprintf(stdout, "\nPE interleaved file: %s\n", infnc);
        fprintf(stdout, "\nTotal input FastQ records: %d (%d pairs)\n", total, (total / 2));
        fprintf(stdout, "\nFastQ paired records kept: %d (%d pairs)\n", kept_p, (kept_p / 2));
        if (pec) fprintf(stdout, "FastQ single records kept: %d\n", (kept_s1 + kept_s2));
        else fprintf(stdout, "FastQ single records kept: %d (from PE1: %d, from PE2: %d)\n", (kept_s1 + kept_s2), kept_s1, kept_s2);

        fprintf(stdout, "FastQ paired records discarded: %d (%d pairs)\n", discard_p, (discard_p / 2));

        if (pec) fprintf(stdout, "FastQ single records discarded: %d\n\n", (discard_s1 + discard_s2));
        else fprintf(stdout, "FastQ single records discarded: %d (from PE1: %d, from PE2: %d)\n\n", (discard_s1 + discard_s2), discard_s1, discard_s2);
    }

    kseq_destroy(fqrec1);
    if (pec) free(fqrec2);
    else kseq_destroy(fqrec2);

    if (sfn && !combo_all) {
        if (!gzip_output) fclose(single);
        else gzclose(single_gzip);
    }

    if (pec) {
        gzclose(pec);
        if (!gzip_output) fclose(combo);
        else gzclose(combo_gzip);
    } else {
        gzclose(pe1);
        gzclose(pe2);
        if (!gzip_output) {
            fclose(outfile1);
            fclose(outfile2);
        } else {
            gzclose(outfile1_gzip);
            gzclose(outfile2_gzip);
        }
    }

    return EXIT_SUCCESS;
}                               /* end of paired_main() */
