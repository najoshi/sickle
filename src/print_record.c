#include <assert.h>
#include <ctype.h>
#include <stdlib.h>
#include <zlib.h>
#include <stdio.h>
#include <unistd.h>
#include "sickle.h"
#include "kseq.h"


void print_record (FILE *fp, kseq_t *fqr, cutsites *cs) {
    fprintf(fp, "@%s", fqr->name.s);
    if (fqr->comment.l) fprintf(fp, " %s\n", fqr->comment.s);
    else fprintf(fp, "\n");
    fprintf(fp, "%.*s\n", cs->three_prime_cut - cs->five_prime_cut, fqr->seq.s + cs->five_prime_cut);
    fprintf(fp, "+\n");
    fprintf(fp, "%.*s\n", cs->three_prime_cut - cs->five_prime_cut, fqr->qual.s + cs->five_prime_cut);
}

void print_record_gzip (gzFile fp, kseq_t *fqr, cutsites *cs) {
    gzprintf(fp, "@%s", fqr->name.s);
    if (fqr->comment.l) gzprintf(fp, " %s\n", fqr->comment.s);
    else gzprintf(fp, "\n");
    gzprintf(fp, "%.*s\n", cs->three_prime_cut - cs->five_prime_cut, fqr->seq.s + cs->five_prime_cut);
    gzprintf(fp, "+\n");
    gzprintf(fp, "%.*s\n", cs->three_prime_cut - cs->five_prime_cut, fqr->qual.s + cs->five_prime_cut);
}

void print_record_N (FILE *fp, kseq_t *fqr, int qualtype) {
    fprintf(fp, "@%s", fqr->name.s);
    if (fqr->comment.l) fprintf(fp, " %s\n", fqr->comment.s);
    else fprintf(fp, "\n");
    fprintf(fp, "N\n+\n%c\n", quality_constants[qualtype][Q_MIN]);
}

void print_record_N_gzip (gzFile fp, kseq_t *fqr, int qualtype) {
    gzprintf(fp, "@%s", fqr->name.s);
    if (fqr->comment.l) gzprintf(fp, " %s\n", fqr->comment.s);
    else gzprintf(fp, "\n");
    gzprintf(fp, "N\n+\n%c\n", quality_constants[qualtype][Q_MIN]);
}

