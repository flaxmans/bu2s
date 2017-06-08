// functions for creating a table to lookup Poisson random numbers
// for the number of recombination events in the BU2S simulation
// goal is to remove calls of log() that were in makeZygoteChromosomes()

// you can compile this file on its own with the following command:
// gcc -lm -I/usr/local/include/gsl -L/usr/local/lib -lgsl -lgslcblas -O3 -o PoissonLookup PoissonLookup.c

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
//#include <gsl_rng.h>            // gnu scientific library //
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


// the "makePoissonLookup()" function: make the lookup table; this would be called ONCE
// during initialization phase of the program:
unsigned int makePoissonLookup(void) {
	// assume poissonTable is a pointer passed in with nothing assigned
	// yet; it will be the basis of the lookup
	// it truncates the Poisson at 3x the expected value
	double expectedNumRecomb;
	unsigned int maxNum, i;
	double truncationFactor = 3.0; // how many times the mean we want to truncate
	
	expectedNumRecomb = TOTAL_MAP_LENGTH / CENTIMORGANS;
	// global TOTAL_MAP_LENGTH used here defined in bu2s.c and declared as extern in bu2s.h
	// CENTIMORGANS defined in bu2s.h = 100.0
	
	// as a cheat, let's allow up to three times that number as the max
	maxNum = (int) ((truncationFactor * expectedNumRecomb) + 0.5); // round up
	poissonTable = (double *) malloc( (maxNum + 1) * sizeof(double) );
	// the "+ 1" allows for an array that has positions corresponding to 0, 1, 2, ..., naxNum
	// (rather than only up to maxNum - 1).
	
	// use the Poisson probability distribution function to set the cutoffs for the lookup table
	for ( i = 0; i < maxNum; i++ )
		poissonTable[i] = gsl_cdf_poisson_P(i, expectedNumRecomb);
	
	// make a safety to ensure we don't run off the end of the table:
	poissonTable[maxNum] = 1.1; // lookup value will never be > 1.0
	
	return(maxNum);
}



int lookupPoissonValue (void) {
	// this function would be called every time a total number of recombinations events was needed
	// this is an inverse CDF function using a truncated Poisson.
	int i = 0;
	double testVal;
	testVal = randU(); // generate value for lookup
	
	while ( testVal > (*(poissonTable + i)) )
		i++;
	
	return(i);
}






