// header file for bu2s.c and associated source files
// started May 24, 2017

// basic parameters that would NEVER be modified on the command line
// placed here to avoid magic numbers (and eventually clean up the bu2s.c file a bit)
#define CENTIMORGANS 100.0 // expected map distance between consecutive recombination events
#define RECORD_LD_VALUES_DEFAULT 0
#define LD_LOCI_SUBSAMPLE 1 // LD only recorded for pairs involving every LD_LOCI_SUBSAMPLEth locus
#define LD_LOWER_BOUND 0.001 // LD only recorded when fabs(DD) > LD_LOWER_BOUND and...


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
extern int nPATCHES;
extern long int nLOCI;
extern long int *locusID;
extern double START_THRESH_FOR_LD;
extern double END_THRESH_FOR_LD;
extern int *variable_loci;
extern int *IS_SELECTED_LOCUS;
extern int *chromosomeMembership;
extern short int *genotypes;
extern double *MAP;
extern long int totalGenerationsElapsed;
extern int N;
extern double *allele_frequencies;
extern double *S_MAX1;



// function prototypes from other .c files:
// these from PoissonLookup.c:
unsigned int makePoissonLookup(void);
int lookupPoissonValue(void);
void getCrossoverLocations(int totalCOcount, double *crossoverLocations);
void myInsertSort(double *a, int l, int r);

// from LDcalculations.c
void calculateLDaverageVals( long int pairCount, int *categoryCounts, int *pairCategories, double *individualLDvals, double *averageLDvalues );
void calculateLDwithinDemes( int gatherLDvalues, double *alleleFrequenciesByPatch, int *patchAbundancePtr );
void calculateLDonePairInDemes( int i, int j, double *afbppti, double *afbpptj, double *LDstorageArray, int *patchAbundancePtr );
void openLDdataFiles( int gatherLDvalues );
void closeLDdataFiles( int gatherLDvalues );


// from bu2s.c
int compare_doubles(const void *a, const void *b);
