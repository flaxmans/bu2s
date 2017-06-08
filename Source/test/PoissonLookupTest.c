// This code is deprecated but was used to test out the Poisson lookup functions
// as written in PoissonLookup.c

// when it was used, it was in that file, and compiled with the following command:
// gcc -lm -I/usr/local/include/gsl -L/usr/local/lib -lgsl -lgslcblas -O3 -o PoissonLookup PoissonLookup.c

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <gsl_rng.h>            // gnu scientific library //
#include <gsl_randist.h>        // gnu scientific library //
#include <gsl_cdf.h>		// needed for making the Poisson lookup table
#include "bu2s.h"
#include "MT/dSFMT.h"


//////////////////////////////
// begin section of TEST CODE:
//
//#define METHOD 1 // "1" for lookup table inverse CDF method; "2" for using the Poisson rng from the gsl
//#define TOTAL_MAP_LENGTH 1000.0 // this is here so we can compile without linking to bu2s.c for now
//const gsl_rng_type *rngType;    /* generator type */
//gsl_rng *rngState;              /* rng instance */
//
//int main(void) {
//	double *poissonTable;
//	int i, maxNum;
//	long int j, bigNum = 100000000;
//	double expectedNumRecomb;
//	FILE *foo;
//
//	gsl_rng_env_setup(); // set up the environment variables for the RNG
//	rngType = gsl_rng_mt19937; // Mersenne Twister// for default you can say rngType = gsl_rng_default;
//	rngState = gsl_rng_alloc(rngType);
//	gsl_rng_set(rngState, 2); // 1 is the seed
//
//
//	if ( METHOD == 1 ) {
//		maxNum = makePoissonLookup();
//		//	for ( i = 0; i <= maxNum; i++ )
//		//		printf("%i\t%E\n", i, poissonTable[i]);
//		foo = fopen("foo.txt", "w"); //
//		for ( j = 0; j < bigNum; j++ ) {
//			i = lookupPoissonValue(poissonTable);
//			if ( j < 100000 )
//				fprintf(foo, "%i\n", i);
//		}
//		fclose(foo);
//	}
//	if ( METHOD == 2 ) {
//		expectedNumRecomb = TOTAL_MAP_LENGTH / CENTIMORGANS;
//		foo = fopen("bar.txt", "w");
//		for ( j = 0; j < bigNum; j++ ) {
//			i = gsl_ran_poisson (rngState, expectedNumRecomb);
//			if ( j < 100000 )
//				fprintf(foo, "%i\n", i);
//		}
//		fclose(foo);
//	}
//
//
//	return 0;
//}
// end of this section of test code
///////////////////////////////

