// header file for bu2s.c and associated source files
// started May 24, 2017

// basic parameters that would NEVER be modified on the command line
// placed here to avoid magic numbers (and eventually clean up the bu2s.c file a bit)
#define CENTIMORGANS 100.0 // expected map distance between consecutive recombination events

#include "MT/dSFMT.h"			// Mersenne Twister library //
// code for using Mersenne Twister RNG
dsfmt_t dsfmt;
extern dsfmt_t dsfmt; // rng state needs to be known everywhere
#define seedRand(s) dsfmt_init_gen_rand(&dsfmt, s)
#define	randU() dsfmt_genrand_close_open(&dsfmt)
#define	randI() (unsigned)dsfmt_genrand_uint32(&dsfmt)

// globals defined in bu2s.c but needed in other files
extern double TOTAL_MAP_LENGTH;
extern double *poissonTable;

// function prototypes from other .c files
// these from PoissonLookup.c:
unsigned int makePoissonLookup(void);
int lookupPoissonValue(void);
void getCrossoverLocations(int totalCOcount, double *crossoverLocations);
void myInsertSort(double *a, int l, int r);

// from bu2s.c
int compare_doubles(const void *a, const void *b);
