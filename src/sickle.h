#ifndef SICKLE_H
#define SICKLE_H

#include <limits.h>
int single_main(int, char**);
int paired_main(int, char**);

#ifndef PROGRAM_NAME
#define PROGRAM_NAME "sickle"
#endif

#ifndef AUTHORS
#define AUTHORS "Nikhil Joshi, UC Davis Bioinformatics Core\n"
#endif

#ifndef VERSION
#define VERSION 0.0
#endif

enum {ILLUMINA_TYPE, PHRED_TYPE, SANGER_TYPE};

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
Usage_call(EXIT_SUCCESS); \
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


#endif /*SICKLE_H*/
