#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <stdint.h>
#include <inttypes.h>
#include "sickle.h"

#if defined(__GNUC__) && defined(__SSE2__)
#	define WITH_SSE2
#endif

#ifdef WITH_SSE2
#	include <xmmintrin.h>
#	include <emmintrin.h>
#endif

#ifdef WITH_SSE2

cutsites* vectorized_sliding_window(kseq_t *fqrec, int qualtype, int lenth_threshold, int qual_threshold, int no_fiveprime, int trunc_n, int debug) {

	cutsites* returnCuts = (cutsites*)malloc(sizeof(cutsites));
	
	if (no_fiveprime) {
		returnCuts->five_prime_cut = check5prime(fqrec->qual.s, qual_threshold + 32, 0, fqrec->seq.l);
	}
	
	returnCuts->three_prime_cut = check3prime(fqrec->qual.s, qual_threshold + 32, 0);
		
	return returnCuts;
}


void print128_num(__m128i test) {
	int64_t *val = (int64_t*) &test;
	printf("%.16llx %.16llx\n", val[1], val[0]);
}

uint32_t sum_array(__m128i v) {
 
	//Creates 128 bit number full of zeros 
    	const __m128i vk0 = _mm_set1_epi8(0);
	//Creates 128 bit number HEX -> 000100010001000100010001
	const __m128i vk1 = _mm_set1_epi16(1);

	
	//sets bits vsum to zero
    	__m128i vsum = _mm_set1_epi32(0);


	//Unpacks the low of the bits  - seperates lowest signicant 8 bits by 00
     	__m128i vl = _mm_unpacklo_epi8(v, vk0);
	//Unpacks the top of the bits  - seperates highest signicant 8 bits by 00
     	__m128i vh = _mm_unpackhi_epi8(v, vk0);


	//Adds low bits
     	vsum = _mm_add_epi32(vsum, _mm_madd_epi16(vl, vk1));
	//Adds high bits
     	vsum = _mm_add_epi32(vsum, _mm_madd_epi16(vh, vk1));
	
    	vsum = _mm_add_epi32(vsum, _mm_srli_si128(vsum, 8));
    	vsum = _mm_add_epi32(vsum, _mm_srli_si128(vsum, 4));

	uint32_t sum = _mm_cvtsi128_si32(vsum);	
	return sum;

	
}

int check3prime(char *qual, uint32_t minQual, uint32_t sumMinQual) {

	
	__m128i loading = _mm_set1_epi32(0);

	uint32_t t = 0;
	uint32_t min = minQual * 16;
	uint32_t i = 0;

	//Loads 16 characters (0 - 15) up in the vectorization function
	//Checks makes sure they are not out of bounds
	
	//Checks the first set (the case nothing has to be trimmed on the 3prime end
	loading = _mm_loadu_si128((const void *)qual);
	t = sum_array(loading);
	
	while (t < min) {
		if (qual[15] != '\n' && qual[15] !='\0') {
			loading = _mm_loadu_si128((const void *)qual);
			t = sum_array(loading);
			i++;
			qual++;
		} else {
			return -1;
		}
	}
	//Do not need to check \0 or newline since we checked it in the vectorization
	//Additionally, this was easier to implement outside the vector	
	if  (qual[i] != '\0' && qual[i] == '\n') { return -1; }

	while (qual[i] <= minQual && qual[i] != '\0' && qual[i] == '\n') {
		i++;	
	}

	return i;
}



int check5prime(char *qual, uint32_t minQual, uint32_t sumMinQual, uint32_t qualLen) {

	qual += qualLen - 17;
	__m128i loading = _mm_set1_epi32(0);

	uint32_t t = 0;
	uint32_t min = minQual * 16;
	int i = qualLen;

	//Loads 16 characters (0 - 15) up in the vectorization function
	//Checks makes sure they are not out of bounds
	
	//Checks the first set (the case nothing has to be trimmed on the 3prime end
	loading = _mm_loadu_si128((const void *)qual);
	t = sum_array(loading);

	while (t < min) {
		loading = _mm_loadu_si128((const void *)qual);
		t = sum_array(loading);
		i--;
		qual--;
	}

	//Do not need to check \0 or newline since we checked it in the vectorization
	//Additionally, this was easier to implement outside the vector	
	uint32_t tmpi = 15;;
	while (qual[tmpi] < minQual) {
		tmpi--;
		i--;	
	}

		
	return i;

}

#endif
