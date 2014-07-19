#include <assert.h>
#include <ctype.h>
#include <stdlib.h>
#include <limits.h>
#include <zlib.h>
#include <stdio.h>
#include <getopt.h>
#include <string.h>
#include "sickle.h"

void main_usage (int status) {

	fprintf (stdout, "\nUsage: %s <command> [options]\n\
\n\
Command:\n\
pe\tpaired-end sequence trimming\n\
se\tsingle-end sequence trimming\n\
\n\
--help, display this help and exit\n\
--version, output version information and exit\n\n", PROGRAM_NAME);

	exit (status);
}

int main (int argc, char *argv[]) {
	int retval=0;

	if (argc < 2 || (strcmp (argv[1],"pe") != 0 && strcmp (argv[1],"se") != 0 && strcmp (argv[1],"--version") != 0 && strcmp (argv[1],"--help") != 0)) {
		main_usage (EXIT_FAILURE);
	}

	if (strcmp (argv[1],"--version") == 0) {
		fprintf(stdout, "%s version %0.2f\nCopyright (c) 2011 The Regents of University of California, Davis Campus.\n%s is free software and comes with ABSOLUTELY NO WARRANTY.\nDistributed under the MIT License.\n\nWritten by %s\n", PROGRAM_NAME, VERSION, PROGRAM_NAME, AUTHORS);

		exit (EXIT_SUCCESS);

	}

	else if (strcmp (argv[1],"--help") == 0) {
		main_usage (EXIT_SUCCESS);
	}

	else if (strcmp (argv[1],"pe") == 0) {
		retval = paired_main (argc, argv);
		return (retval);
	}

	else if (strcmp (argv[1],"se") == 0) {
		retval = single_main (argc, argv);
		return (retval);
	}

	return 0;
}
