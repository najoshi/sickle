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

static const char *single_short_options = "df:t:o:q:l:w:zxng";
static struct option single_long_options[] = {
    { "fastq-file",         required_argument, NULL, 'f' },
    { "output-file",        required_argument, NULL, 'o' },
    { "qual-type",          required_argument, NULL, 't' },
    { "qual-threshold",     required_argument, NULL, 'q' },
    { "length-threshold",   required_argument, NULL, 'l' },
    { "window",             required_argument, NULL, 'w' },
    { "no-fiveprime",       no_argument,       NULL, 'x' },
    { "truncate-n",         no_argument,       NULL, 'n' },
    { "gzip-output",        no_argument,       NULL, 'g' },
    { "quiet",              no_argument,       NULL, 'z' },
    { GETOPT_HELP_OPTION_DECL },
    { GETOPT_VERSION_OPTION_DECL },
    { NULL, 0, NULL, 0 }
};

void single_usage(int status, char *msg) {
    static const char *usage_format =
        "\n"
        "Usage: %1$s se [options] -f <fastq sequence file>\n"
        "       %2$*3$s           -t <quality type>\n"
        "       %2$*3$s           -o <trimmed fastq file>\n"
        "\n"
        "Options:\n"
        "\n"
        "-f FILE, --fastq-file FILE   Input fastq file (required)\n"
        "-o FILE, --output-file FILE  Output trimmed fastq file (required)\n"
        "-t TYPE, --qual-type TYPE    Type of quality values, one of:\n"
        "                                 solexa (CASAVA < 1.3)\n"
        "                                 illumina (CASAVA 1.3 to 1.7)\n"
        "                                 sanger (which is CASAVA >= 1.8)\n"
        "                             (required)\n"
        "-q #, --qual-threshold #     Threshold for trimming based on average quality\n"
        "                             in a window. Default %3$d.\n"
        "-l #, --length-threshold #   Threshold to keep a read based on length after\n"
        "                             trimming. Default %4$d.\n"
        "-w #, --window #             Fixed window size to use.  Default is a dynamic\n"
        "                             window size of 0.1 of the read length.\n"
        "-x, --no-fiveprime           Don't do five prime trimming.\n"
        "-n, --truncate-n             Truncate sequences at position of first N.\n"
        "-g, --gzip-output            Output gzipped files.\n"
        "--quiet                      Do not output trimming info\n"
        "--help                       Display this help and exit\n"
        "--version                    Output version information and exit\n"
    ;

    if (msg) fprintf( STDERR_OR_OUT(status), "%s\n", msg );

    fprintf(
        STDERR_OR_OUT(status),
        usage_format,
        PROGRAM_NAME,
        "", strlen(PROGRAM_NAME) + 3, // +3 makes it possible to keep text visually aligned in the format string itself
        single_qual_threshold,
        single_length_threshold
    );

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
    int window_size = 0;

    while (1) {
        int option_index = 0;
        optc = getopt_long(argc, argv, single_short_options, single_long_options, &option_index);

        if (optc == -1)
            break;

        switch (optc) {
            if (single_long_options[option_index].flag != 0)
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

        case 'w':
            window_size = atoi(optarg);
            if (window_size < 1) {
                fprintf(stderr, "Fixed window size must be >= 1\n");
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
        outfile = fopen(outfn, "w");
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

        p1cut = sliding_window(fqrec, qualtype, single_length_threshold, single_qual_threshold, window_size, no_fiveprime, trunc_n, debug);
        total++;

        if (debug) printf("P1cut: %d,%d\n", p1cut->five_prime_cut, p1cut->three_prime_cut);

        /* if sequence quality and length pass filter then output record, else discard */
        if (p1cut->three_prime_cut >= 0) {
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

    if (!quiet) fprintf(stdout, "\nSE input file: %s\n\nTotal FastQ records: %d\nFastQ records kept: %d\nFastQ records discarded: %d\n\n", infn, total, kept, discard);

    kseq_destroy(fqrec);
    gzclose(se);
    if (!gzip_output) fclose(outfile);
    else gzclose(outfile_gzip);

    return EXIT_SUCCESS;
}
