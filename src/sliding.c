#include <assert.h>
#include <ctype.h>
#include <stdlib.h>
#include <zlib.h>
#include <stdio.h>
#include <getopt.h>
#include "sickle.h"
#include "kseq.h"

int get_quality_num (char qualchar, int qualtype) {
  /* 
     Return the adjusted quality, depending on quality type.

     Note that this uses the array in sickle.h, which *approximates*
     the SOLEXA (pre-1.6 pipeline) qualities as linear. This is
     inaccurate with low-quality bases.
  */
  return((int) qualchar - quality_constants[qualtype][Q_OFFSET]);
}


int sliding_window (kseq_t *fqrec, int qualtype, int length_threshold, int qual_threshold) {

	int window_size = (int) (0.1 * fqrec->seq.l);
	int i,j;
	int window_start=0;
	int window_total=0;
	int p1cut = fqrec->seq.l;

	/* if the seq length is less then 10bp, */
	/* then make the window size the length of the seq */
	if (window_size == 0) window_size = fqrec->seq.l;

	for (i=0; i<window_size; i++) {
		window_total += get_quality_num (fqrec->qual.s[i], qualtype);
	}

	for (i=0; i<fqrec->qual.l; i++) {

		/* if the average quality in the window is less than the threshold */
		/* or if the window is the last window in the read */
		if (((double)window_total / (double)window_size < qual_threshold) || 
			(window_start+window_size > fqrec->qual.l)) {

			/* at what point in the window does the quality dip below the threshold? */
			for (j=window_start; j<window_start+window_size; j++) {
				if (get_quality_num (fqrec->qual.s[j], qualtype) < qual_threshold) {
					p1cut = j;

					/* if cutting length is less than threshold then return -1 */
					/* to indicate that the read should be discarded */
					if (p1cut < length_threshold) {p1cut = -1;}
					break;
				}
			}

			break;
		}

		/* instead of sliding the window, subtract the first qual and add the next qual */
		window_total -= get_quality_num (fqrec->qual.s[window_start], qualtype);
		window_total += get_quality_num (fqrec->qual.s[window_start+window_size], qualtype);
		window_start++;
	}

	return p1cut;
}
