// header file for bu2s.c and associated source files
// started May 24, 2017

// basic parameters that would NEVER be modified on the command line
// placed here to avoid magic numbers (and eventually clean up the bu2s.c file a bit)
#define CENTIMORGANS 100.0 // expected map distance between consecutive recombination events

// globals defined in bu2s.c but needed in other files
extern double TOTAL_MAP_LENGTH;
extern gsl_rng *rngState;

// function prototypes from other .c files
// these from PoissonLookup.c:
double * makePoissonLookup(int *maxNumAddress);
int lookupPoissonValue(double *poissonTable);
