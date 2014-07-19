#ifndef PRINT_RECORD_H
#define PRINT_RECORD_H

#include <stdio.h>
#include <zlib.h>
#include "kseq.h"

void print_record (FILE *fp, kseq_t *fqr, cutsites *cs);
void print_record_gzip (gzFile fp, kseq_t *fqr, cutsites *cs);
void print_record_N (FILE *fp, kseq_t *fqr, int qualtype);
void print_record_N_gzip (gzFile fp, kseq_t *fqr, int qualtype);

#endif /* PRINT_RECORD_H */
