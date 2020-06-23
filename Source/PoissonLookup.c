// functions for creating a table to lookup Poisson random numbers
// for the number of recombination events in the BU2S simulation
// goal is to remove calls of log() that were in makeZygoteChromosomes()

#include <stdlib.h>
#include <gsl/gsl_randist.h>        // gnu scientific library //
#include <gsl/gsl_cdf.h>			// needed for making the Poisson lookup table
#include "bu2s.h"				// header file for bu2s program

#define SORT_THRESHOLD 60
void QuickInsert(double *a, int l, int r);

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
		if (totalCOcount < SORT_THRESHOLD)
			myInsertSort(crossoverLocations, 0, totalCOcount-1);
		else
			QuickInsert(crossoverLocations, 0, totalCOcount-1);

		crossoverLocations[totalCOcount] = TOTAL_MAP_LENGTH + 1.0;
	}
	else
		crossoverLocations[0] = TOTAL_MAP_LENGTH + 1.0;
}


void
myInsertSort(double *a, int l, int r)
{
    int i, j;
    double v;

    for (i = l+1; i <= r; i++) {
        v = a[i];

        for (j = i; j > l && a[j-1] > v; j--)
            a[j] = a[j-1];

        a[j] = v;
    }
}

#define swap(x, y) { double t = y; y = x; x = t; }


//
// QuickSort that reverts to insertion sort when size is below a threshold.
// Optimized for bu2s: sorting random numbers that aren't already in order.
// Allows a little more efficiency to not handle the cases where QuickSort
// goes to O(N^2).
//
void
QuickInsert(double *a, int l, int r)
{
    int m;
    double v;

    if ((r - l) < SORT_THRESHOLD) {     // Insertion sort is faster for small arrays.
        myInsertSort(a, l, r);
        return;
    }

    //
    // Find median of first three elements. Aiming to get a number close to the
    // midpoint to get the best performance.
    //
    // Note: We're guaranteed at this point that there are at least N elements in the array
    // because of the preceding code that switches to an insertion sort for small arrays.
    //
    if (a[l] > a[l+1])
        swap(a[l], a[l+1]);
    if (a[l] > a[l+2])
        swap(a[l], a[l+2]);
    if (a[l+1] > a[l+2])
        swap(a[l+1], a[l+2]);

    //
    // At this point, we know a[l] <= a[l+1] <= a[l+2]
    // median is a[l+1]
    //
    m = l+1;
    v = a[m];

    for (int i = l+3; i <= r; i++) {
        if (a[i] < v) {
            m++;                // can't put this in the swap() call because it's a macro
            swap(a[i], a[m]);
        }
    }

    swap(a[l+1], a[m]);
    QuickInsert(a, l, m-1);
    QuickInsert(a, m+1, r);
}
