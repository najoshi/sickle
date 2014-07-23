#ifndef SICKLE_H
#define SICKLE_H

#include <limits.h>
#include <zlib.h>
#include "kseq.h"


/* KSEQ_INIT() cannot be called here, because we only need the types
   defined. Calling KSEQ_INIT() would also define functions, leading
   to an unused function warning with GCC. So, the basic typedefs
   kseq.h has are included here, and each file that reads needs:

   __KS_GETC(gzread, BUFFER_SIZE)
   __KS_GETUNTIL(gzread, BUFFER_SIZE)
   __KSEQ_READ

*/

#define BUFFER_SIZE 4096
__KS_TYPE(gzFile)
__KS_BASIC(gzFile, BUFFER_SIZE)
__KSEQ_TYPE(gzFile)
__KSEQ_BASIC(gzFile)

#ifndef PROGRAM_NAME
#define PROGRAM_NAME "sickle"
#endif

#ifndef AUTHORS
#define AUTHORS "Nikhil Joshi, UC Davis Bioinformatics Core\n"
#endif

#ifndef VERSION
#define VERSION 0.0
#endif

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
#define case_GETOPT_HELP_CHAR(Usage_call) \
case GETOPT_HELP_CHAR: \
Usage_call(EXIT_SUCCESS, NULL); \
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

typedef enum {
  PHRED,
  SANGER,
  SOLEXA,
  ILLUMINA
} quality_type;

static const char typenames[4][10] = {
	{"Phred"},
	{"Sanger"},
	{"Solexa"},
	{"Illumina"}
};

#define Q_OFFSET 0
#define Q_MIN 1
#define Q_MAX 2

static const int quality_constants[4][3] = {
  /* offset, min, max */
  {0, 4, 60}, /* PHRED */
  {33, 33, 126}, /* SANGER */
  {64, 58, 112}, /* SOLEXA; this is an approx; the transform is non-linear */
  {64, 64, 110} /* ILLUMINA */
};

typedef struct __cutsites_ {
    int five_prime_cut;
	int three_prime_cut;
} cutsites;


/* Function Prototypes */
int single_main (int argc, char *argv[]);
int paired_main (int argc, char *argv[]);
cutsites* sliding_window (kseq_t *fqrec, int qualtype, int length_threshold, int qual_threshold, int no_fiveprime, int trunc_n, int debug);

#endif /*SICKLE_H*/
