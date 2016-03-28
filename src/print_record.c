#include <assert.h>
#include <ctype.h>
#include <stdlib.h>
#include <zlib.h>
#include <stdio.h>
#include <unistd.h>
#include "sickle.h"
#include "kseq.h"

void print_record_tab (FILE *fp, kseq_t *fqr1, kseq_t *fqr2, cutsites *cs1, cutsites *cs2) {

    fprintf(fp, "@%s\t", fqr1->name.s);
    fprintf(fp, "%.*s\t", cs1->three_prime_cut - cs1->five_prime_cut, fqr1->seq.s + cs1->five_prime_cut);
    fprintf(fp, "%.*s\t", cs1->three_prime_cut - cs1->five_prime_cut, fqr1->qual.s + cs1->five_prime_cut);
    fprintf(fp, "%.*s\t", cs2->three_prime_cut - cs2->five_prime_cut, fqr2->seq.s + cs2->five_prime_cut);
    fprintf(fp, "%.*s\n", cs2->three_prime_cut - cs2->five_prime_cut, fqr2->qual.s + cs2->five_prime_cut);
}


void print_record (FILE *fp, kseq_t *fqr, cutsites *cs) {
    fprintf(fp, "@%s", fqr->name.s);
    if (fqr->comment.l) fprintf(fp, " %s\n", fqr->comment.s);
    else fprintf(fp, "\n");
    fprintf(fp, "%.*s\n", cs->three_prime_cut - cs->five_prime_cut, fqr->seq.s + cs->five_prime_cut);
    fprintf(fp, "+\n");
    fprintf(fp, "%.*s\n", cs->three_prime_cut - cs->five_prime_cut, fqr->qual.s + cs->five_prime_cut);
}
void print_record_tab_s (FILE *fp, kseq_t *fqr, cutsites *cs) {
    fprintf(fp, "@%s\t", fqr->name.s);
    fprintf(fp, "%.*s\t", cs->three_prime_cut - cs->five_prime_cut, fqr->seq.s + cs->five_prime_cut);
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

