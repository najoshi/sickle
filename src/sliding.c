#include <assert.h>
#include <ctype.h>
#include <stdlib.h>
#include <zlib.h>
#include <stdio.h>
#include <getopt.h>
#include <pthread.h>
#include <string.h>
#include "sickle.h"

int get_quality_num (char qualchar, int qualtype, fqrec *fqrec, int pos) {
  /* 
     Return the adjusted quality, depending on quality type.

     Note that this uses the array in sickle.h, which *approximates*
     the SOLEXA (pre-1.3 pipeline) qualities as linear. This is
     inaccurate with low-quality bases.
  */

  int qual_value = (int) qualchar;

  if (qual_value < quality_constants[qualtype][Q_MIN] || qual_value > quality_constants[qualtype][Q_MAX]) {
	fprintf (stderr, "ERROR: Quality value (%d) does not fall within correct range for %s encoding.\n", qual_value, typenames[qualtype]);
	fprintf (stderr, "Range for %s encoding: %d-%d\n", typenames[qualtype], quality_constants[qualtype][Q_MIN], quality_constants[qualtype][Q_MAX]);
	fprintf (stderr, "FastQ record: %s\n", fqrec->h1);
	fprintf (stderr, "Quality string: %s\n", fqrec->qual);
	fprintf (stderr, "Quality char: '%c'\n", qualchar);
	fprintf (stderr, "Quality position: %d\n", pos+1);
	exit(1);
  }

  return (qual_value - quality_constants[qualtype][Q_OFFSET]);
}

fqrec* get_fq_record (FILE *fp) {
    char *h1 = NULL;
    size_t h1size = 0;
    ssize_t h1len;

    char *seq = NULL;
    size_t seqsize = 0;
    ssize_t seqlen;

    char *h2 = NULL;
    size_t h2size = 0;
    ssize_t h2len;

    char *qual = NULL;
    size_t qsize = 0;
    ssize_t qlen;

    fqrec* retval;

/* printf ("inside get_fq_record\n"); */

    h1len = getline (&h1, &h1size, fp);
    if (h1len == -1) {return NULL;}

/* printf ("got first line\n"); */

    seqlen = getline (&seq, &seqsize, fp);
    h2len = getline (&h2, &h2size, fp);
    qlen = getline (&qual, &qsize, fp);

    if (seq[seqlen - 1] == '\n') {
        seq[seqlen - 1] = '\0';
        seqlen--;
    }

    if (qual[qlen - 1] == '\n') {
        qual[qlen - 1] = '\0';
        qlen--;
    }

    retval = malloc (sizeof(fqrec));
    retval->h1 = h1;
    retval->seq = seq;
    retval->h2 = NULL;
    retval->qual = qual;
    retval->seqlen = seqlen;
    retval->qlen = qlen;

/* printf ("got here inside get fq record: %s\n", seq); */

    return (retval);
}

void* sw_call (void *tn) {
    int *thread_num;
    thread_num = (int *) tn;

/* printf ("got inside sw_call: thread_num is %d, waitnum is %d\n", *thread_num, waitnum);
printf ("Inside sw_call: %s, %d\n", global_swd[*thread_num].fqrec_r1->seq, global_swd[*thread_num].length_threshold);
*/
    sliding_window (global_swd[*thread_num].fqrec_r1, global_swd[*thread_num].qualtype, global_swd[*thread_num].length_threshold, global_swd[*thread_num].qual_threshold, global_swd[*thread_num].no_fiveprime, global_swd[*thread_num].discard_n, &(global_swd[*thread_num].five_prime_cut_r1), &(global_swd[*thread_num].three_prime_cut_r1));
    sliding_window (global_swd[*thread_num].fqrec_r2, global_swd[*thread_num].qualtype, global_swd[*thread_num].length_threshold, global_swd[*thread_num].qual_threshold, global_swd[*thread_num].no_fiveprime, global_swd[*thread_num].discard_n, &(global_swd[*thread_num].five_prime_cut_r2), &(global_swd[*thread_num].three_prime_cut_r2));

    pthread_exit(0);
}

void sliding_window (fqrec *fqrec, int qualtype, int length_threshold, int qual_threshold, int no_fiveprime, int discard_n, int *five_prime_cut, int *three_prime_cut) {

	int window_size = (int) (0.1 * fqrec->seqlen);
	int i,j;
	int window_start=0;
	int window_total=0;
    int found_five_prime = 0;
    double window_avg;
	*three_prime_cut = fqrec->seqlen;
	*five_prime_cut = 0;
/*
printf ("got here in sliding window1: %s\n", fqrec->seq);
if (fqrec == NULL) {printf ("fqrec is NULL!!!\n");}
*/
	/* If the sequence contains an "N" then discard if the option has been selected */
	/* Also discard if the length of the sequence is less than the length threshold */
	if ((discard_n && (strstr(fqrec->seq, "N") || strstr(fqrec->seq, "n"))) || 
			(fqrec->seqlen < length_threshold)) {
		*three_prime_cut = -1;
		*five_prime_cut = -1;
	}

/* printf ("got here in sliding window2\n"); */

	/* if the seq length is less then 10bp, */
	/* then make the window size the length of the seq */
	if (window_size == 0) window_size = fqrec->seqlen;

	for (i=0; i<window_size; i++) {
		window_total += get_quality_num (fqrec->qual[i], qualtype, fqrec, i);
	}

	for (i=0; i <= fqrec->qlen - window_size; i++) {

		window_avg = (double)window_total / (double)window_size;

		/* Finding the 5' cutoff */
		/* Find when the average quality in the window goes above the threshold starting from the 5' end */
		if (no_fiveprime == 0 && found_five_prime == 0 && window_avg >= qual_threshold) {

			/* at what point in the window does the quality go above the threshold? */
			for (j=window_start; j<window_start+window_size; j++) {
				if (get_quality_num (fqrec->qual[j], qualtype, fqrec, j) >= qual_threshold) {
					*five_prime_cut = j;
					break;
				}
			}

			found_five_prime = 1;
		}

		/* Finding the 3' cutoff */
		/* if the average quality in the window is less than the threshold */
		/* or if the window is the last window in the read */
		if ((window_avg < qual_threshold) || 
			(window_start+window_size > fqrec->qlen)) {

			/* at what point in the window does the quality dip below the threshold? */
			for (j=window_start; j<window_start+window_size; j++) {
				if (get_quality_num (fqrec->qual[j], qualtype, fqrec, j) < qual_threshold) {
					*three_prime_cut = j;

					/* if cutting length is less than threshold then return -1 for both */
					/* to indicate that the read should be discarded */
					if (*three_prime_cut - *five_prime_cut < length_threshold) {
						*three_prime_cut = -1;
						*five_prime_cut = -1;
					}
					break;
				}
			}

			break;
		}

		/* instead of sliding the window, subtract the first qual and add the next qual */
		window_total -= get_quality_num (fqrec->qual[window_start], qualtype, fqrec, window_start);
		if (window_start+window_size < fqrec->qlen) {
			window_total += get_quality_num (fqrec->qual[window_start+window_size], qualtype, fqrec, window_start+window_size);
		}
		window_start++;
	}
}
