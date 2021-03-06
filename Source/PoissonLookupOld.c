// functions for creating a table to lookup Poisson random numbers
// for the number of recombination events in the BU2S simulation
// goal is to remove calls of log() that were in makeZygoteChromosomes()

#include <stdlib.h>
#include <gsl_randist.h>        // gnu scientific library //
#include <gsl_cdf.h>			// needed for making the Poisson lookup table
#include "bu2s.h"				// header file for bu2s program

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


void getCrossoverLocations(int totalCOcount, double *crossoverLocations) {
	
	int i;
	
	if ( totalCOcount > 0 ) {
		//dsfmt_fill_array_open_open( &dsfmt, crossoverLocations, totalCOcount); // get a vector of random numbers
		for ( i = 0; i < totalCOcount; i++ ) {
			crossoverLocations[i] = randU() * TOTAL_MAP_LENGTH; // multiplication for scale
		}
		// two possible sorting algorithms; insertion sort may be marginally faster based on limited testing
		//qsort(crossoverLocations, totalCOcount, sizeof(double), compare_doubles);
		myInsertSort(crossoverLocations, totalCOcount);
		crossoverLocations[totalCOcount] = TOTAL_MAP_LENGTH + 1.0;
	}
	else
		crossoverLocations[0] = TOTAL_MAP_LENGTH + 1.0;
}


void myInsertSort(double * array, int numElements) {
	
	// code modified from: http://www.programmingsimplified.com/c/source-code/c-program-insertion-sort
	// accessed 6/9/17
	int c, d, d1;
	double foo;
	
	for (c = 1; c < numElements; c++) {
		
		d = c;
		d1 = d - 1;
		
		while ( d > 0 && *(array + d) < *(array + d1) ) {
			foo = *(array + d);
			*(array + d) = *(array + d1);
			*(array + d1) = foo;
			d--;
			d1--;
		}
		
	}
}



