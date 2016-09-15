/*
 *  bu2s.c
 *
 *
 *  Created by Samuel Melvin Flaxman on 4/10/12.
 *  Copyright 2012-2016 Samuel Melvin Flaxman. All rights reserved.
 *
 */

const char *version = "bu2s_3.6.1";

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <getopt.h>
#include <time.h>
#include <sys/time.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>

#include "MT/dSFMT.h"

int loci_count_0 = 0;
int loci_count_1 = 0;
int loci_count_2 = 0;
int loci_count_many = 0;

// code for using Mersenne Twister RNG
dsfmt_t dsfmt;
#define seedRand(s) dsfmt_init_gen_rand(&dsfmt, s)
#define	randU() dsfmt_genrand_close_open(&dsfmt)
#define	randI() (unsigned)dsfmt_genrand_uint32(&dsfmt)


// constants that have to be set here in source code
#define D 1 // dimensionality of habitat; set to 1 or 2
#define PATCHES 2 // number of patches across the habitat in a single dimension; see also nPATCHES below, which is the actual number = PATCHES^D
#define USE_MUTATIONS_FROM_FILES 0 // whether mutation sequence will be defined by a file or at random
#define ADDITIVE_FITNESS 0
#define MULTIPLICATIVE_FITNESS 1

//#define KEEP_REINTRODUCING 0
//#define nMUTATIONS_TO_ESTABLISH 10
//#define nMAX_TRIES_TO_ESTABLISH 5000


// parameters that can be changed with command line options, though some options have not yet been built....
#define TWO_DEME_DEFAULT 1 // if this is 1, then D above should be 1 and PATCHES should be 2, and space will be treated as discrete rather than as continuous
#define DETERMINISTIC_DEFAULT 1
#define nGENERATIONS_MAX_DEFAULT 1 // number of generations between mutations
#define nMUTATIONS_DEFAULT 15000000
#define MUTATIONS_PER_GENERATION_DEFAULT 1
#define INITIAL_CONDITIONS_DEFAULT 0
#define INITIAL_POPULATION_SIZE_DEFAULT 4000
#define MEAN_S_DEFAULT 0.01 // mean value of selection coefficients drawn from -exponential distribution
#define SYMMETRIC_MUTATIONS_DEFAULT 1
#define MAP_TYPE_DEFAULT 1
#define GAMETE_PRODUCTION_MODE_DEFAULT 0 // 0 for DH and GH, 1 for GH only, 2 for "bag of genes"
#define FIXED_N_DEFAULT 1 // whether N will be fixed or variable (latter would be density dependent regulation)
#define K_DEFAULT 10000.0 // default carrying capacity per patch; only used if FIXED_N = 0
#define H_DEFAULT 0.5 // dominance coefficient
#define MOSAIC_DEFAULT 0
#define SD_MOVE_DEFAULT 0.1 // if TWO_DEME is used, this IS the gross migration rate
#define OIRL_DEFAULT 1 // offspring in random locations or not
#define nCHROMOSOMES_DEFAULT 10 // only used if MAP_TYPE == 1
#define TOTAL_MAP_LENGTH_DEFAULT 1000.0 // only used if MAP_TYPE == 1
#define MUTATION_DISTRIBUTION_DEFAULT 0 // 0 for exponential, 1 for fattened tail, 2 for early "flat" (uniform) distribution, 3 for constant S
#define FATTEN_TAIL_PROPORTION_DEFAULT 0.0001 // proportion to flatten or fatten
#define FATTEN_TAIL_MAX_DEFAULT 1.0 // max value to use when flattening or fattening
#define NUMBER_BIG_TO_ESTABLISH_DEFAULT 0
#define START_WITH_BIG_IN_PLACE_DEFAULT 0 // boolean for starting with a few mutations having S = BIG_START_S_VAL to try to see if they will provide seeds for "islands"
#define BIG_START_S_VAL_DEFAULT 0.5 // if the previous bool is set to 1, then this is used
// #define SECONDARY_CNTCT_MUTATIONS_DEFAULT -1 // set to -1 for primary contact; set to > 0 for secondary contact after a certain number of mutations.
#define END_PERIOD_ALLOPATRY_DEFAULT -1 // replaces SECONDARY_CNTCT_MUTATIONS_DEFAULT
#define PERIOD_ALLOPATRY 1000
#define START_PERIOD_ALLOPATRY_DEFAULT -1
#define DMI_MEAN_EFFECT_DEFAULT 0.0 // proportional reduction in fitness arising from negative epistatic interaction
#define DMI_PROB_DEFAULT 0.0 // proportion of pairs of loci at which a negative epistatic interaction arises.
// set this probability to < 0 for no epistatis
#define PROBABILITY_DERIVED_REVERSED_DEFAULT 0.5
#define POSITIVE_EPI_PROB_DEFAULT 0.0 // probability that a pair of loci can have a positive epistatic interaction
#define MEAN_POSITIVE_EPI_COEFF_DEFAULT 0.0 // factor to enhance pairwise derived loci by = 1 + this

#define DEME0_CONSTANT_S_DEFAULT 0.01
#define DEME1_CONSTANT_S_DEFAULT 0.01

#define FRACTION_SELECTED_LOCI_DEFAULT 0.090909 // fraction of mutations subject to divergent selection; rest are neutral

// parameters that affect data recording only
#define nTIME_SAMPLES_DEFAULT 100
#define RECORD_LD_VALUES_DEFAULT 0
double START_THRESH_FOR_LD = 0.87;
double END_THRESH_FOR_LD = 0.06;
// #define LD_SAMPLING_FREQUENCY 25 // time subsampling of LD
#define LD_LOCI_SUBSAMPLE 1 // LD only recorded for pairs involving every LD_LOCI_SUBSAMPLEth locus
#define AL_FREQ_SUBSAMPLE 1
#define LD_LowerBound 0.001 // LD only recorded when fabs(DD) > LD_LowerBound and...
//#define LD_UpperBound 1.0 // LD only recorded when fabs(Delta) < LD_UpperBound
double RI_THRESH = 0.0001; // threshold for ending simulations
#define RECORD_FIT_TS_DEFAULT 1 // data recording Boolean
int nRECORD_FIT = 100;
#define FST_MIN_RECORDING_THRESH -1.0 // minimum value of FST, below which, values for a locus will NOT be recorded in data files
#define MIN_TS_SAMP_FREQ 100
#define MAX_TS_SAMP_FREQ 1000
#define RECORDING_TIMES_IN_FILE_DEFAULT 0
#define PURE_NEUTRAL_MODEL_DEFAULT 0
#define PURE_ALLOPATRY_MODEL_DEFAULT 0

// globals whose values can be alterd on command line
int DETERMINISTIC = DETERMINISTIC_DEFAULT;
long int nGENERATIONS_MAX = nGENERATIONS_MAX_DEFAULT;
int nMUTATIONS = nMUTATIONS_DEFAULT;
int INITIAL_CONDITIONS = INITIAL_CONDITIONS_DEFAULT;
double MEAN_S = MEAN_S_DEFAULT;
_Bool SYMMETRIC_MUTATIONS = SYMMETRIC_MUTATIONS_DEFAULT; // whether 0 and 1 alleles will have symmetric, opposite effects
int MAP_TYPE = MAP_TYPE_DEFAULT;
int nCHROMOSOMES = nCHROMOSOMES_DEFAULT;   // total number of chromosomes
int GAMETE_PRODUCTION_MODE = GAMETE_PRODUCTION_MODE_DEFAULT;
_Bool FIXED_N = FIXED_N_DEFAULT;
int INITIAL_POPULATION_SIZE = INITIAL_POPULATION_SIZE_DEFAULT;
double K = K_DEFAULT;
double HVAL = H_DEFAULT;
_Bool MOSAIC = MOSAIC_DEFAULT;
double SD_MOVE = SD_MOVE_DEFAULT;
int nTIME_SAMPLES = nTIME_SAMPLES_DEFAULT;
long int TS_SAMPLING_FREQUENCY;
_Bool OFFSPRING_IN_RANDOM_LOCATIONS = OIRL_DEFAULT;
double TOTAL_MAP_LENGTH = TOTAL_MAP_LENGTH_DEFAULT;
int MUTATION_DISTRIBUTION = MUTATION_DISTRIBUTION_DEFAULT;
double FATTEN_TAIL_PROPORTION = FATTEN_TAIL_PROPORTION_DEFAULT;  // only used if MUTATION_DISTRIBUTION == 1
double FATTEN_TAIL_MAX = FATTEN_TAIL_MAX_DEFAULT; // only used if MUTATION_DISTRIBUTION == 1
int NUMBER_BIG_TO_ESTABLISH = NUMBER_BIG_TO_ESTABLISH_DEFAULT;
_Bool TWO_DEME = TWO_DEME_DEFAULT;
_Bool RECORD_LD_VALUES = RECORD_LD_VALUES_DEFAULT;
// int SECONDARY_CNTCT_MUTATIONS = SECONDARY_CNTCT_MUTATIONS_DEFAULT;
long int END_PERIOD_ALLOPATRY = END_PERIOD_ALLOPATRY_DEFAULT;
// double EM_THRESH_FOR_LD; // threshold for recording of LD values
_Bool BeginRecordingLD = 0;
_Bool START_WITH_BIG_IN_PLACE = START_WITH_BIG_IN_PLACE_DEFAULT; // boolean for starting with a few mutations having S = BIG_START_S_VAL to try to see if they will provide seeds for "islands"
double BIG_START_S_VAL = BIG_START_S_VAL_DEFAULT;
double DMI_MEAN_EFFECT = DMI_MEAN_EFFECT_DEFAULT;
double DMI_PROBABILITY = DMI_PROB_DEFAULT;
int EARLY_TEST_KILL = 0;
double PROBABILITY_DERIVED_REVERSED = PROBABILITY_DERIVED_REVERSED_DEFAULT;
double POSITIVE_EPI_PROBABILITY = POSITIVE_EPI_PROB_DEFAULT;
double MEAN_POSITIVE_EPI_COEFF = MEAN_POSITIVE_EPI_COEFF_DEFAULT;
_Bool RECORD_FIT_TS = RECORD_FIT_TS_DEFAULT;
long int START_PERIOD_ALLOPATRY = START_PERIOD_ALLOPATRY_DEFAULT;
int MUTATIONS_PER_GENERATION = MUTATIONS_PER_GENERATION_DEFAULT;
double DEME0_CONSTANT_S = DEME0_CONSTANT_S_DEFAULT;
double DEME1_CONSTANT_S = DEME1_CONSTANT_S_DEFAULT;
double FRACTION_SELECTED_LOCI = FRACTION_SELECTED_LOCI_DEFAULT;
_Bool PURE_NEUTRAL_MODEL = PURE_NEUTRAL_MODEL_DEFAULT;
_Bool PURE_ALLOPATRY_MODEL = PURE_ALLOPATRY_MODEL_DEFAULT;



// pointers for malloc() calls
double *x_locations;  // array of length POPULAION_SIZE storing a coordinate
#if (D == 2)
double *y_locations;  // array of length POPULAION_SIZE storing a coordinate
#endif
short int *genotypes;	// array of length ( POPULATION_SIZE * nLOCI * 2 ) storing genotypes of all individual in the population
int *patch_locations;  // array of length POULATION_SIZE storing integer patch number to which each individual belongs
double *MUTATION_ORDER;  // array of lenght nMUTATIONS that stores the indexes of loci where "1" mutations arise, and the order in which this happens
int *MUTATION_TYPE_SEQUENCE; // array of length nMUTATIONS that stores a "1" if the mutation was selected and a zero otherwise
int *MUTATION_LOCATIONS; // array of length nMUTATIONS that stores the index of the PATCH were the ith mutation arises
int *LOCI_PER_CHROMOSOME; // array of length nCHROMOSOMES that stores the number of loci on each of the nCHROMOSOMES
double *MAP_LENGTHS;  // array of length nCHROMOSOMES that stores the map length of each chromosome
int *previousPatches; // array to help keep records of effective migration rate
double *S_MAX1_SEQUENCE;
double *S_MAX0_SEQUENCE;
double *H_SEQUENCE;
long int *locusID;
double *allele_frequencies; // frequency of "1" allele at each locus across the whole population
double *MAP; // genetic map of all chromosomes; each zero is the start of a new chromosome
double *S_MAX1; // values of selection coefficients for 1 alleles
double *S_MAX0; // values of selection coefficients for 0 alleles; may be symmetric with S_MAX1
int *IS_SELECTED_LOCUS;
int *variable_loci; // record which loci have more than one allele at a given time
short int *fixed_allele; // will only be used for loci that are NOT variable
int *chromosomeMembership;
double *H;
int *is_reference_locus;
int totalMigrants;
double *fit0;
double *fit1;
double *fitH;
short int *epistasisMatrix;
_Bool *is_reversed_locus;
_Bool *REVERSAL_SEQUENCE;
_Bool *any_epis_for_this_locus;
double *epi_coeff_matrix;


// other globals
long int nLOCI = 1;
int N;              // total population size
int nSELECTED_LOCI = 0;
int nNEUTRAL_LOCI = 0;
long int m = 0;     // mutation counter
long int t;
const double TWOPI = 2.0 * M_PI;
long int numberOfMutationTries = 0;
int nVariableLoci = 0;
long int totalGenerationsElapsed = 0;
_Bool RI_REACHED = 0; // threshold for ending a run based upon very low effective migration rates as set by RI_THRESH
long int newestLocus;
long int totalDMIs = 0;
long int totalPositiveEpis = 0;
_Bool CONSIDER_EPISTASIS;
long int nPositiveEpisReversed = 0;
long int nPositiveEpisRegular = 0;
long int nRegularLoci = 0;
long int nReversedLoci = 0;
long int lastTimeResNearImm = -1;
long int lastTimeResNearRand = -1;
long int firstTimeResNearMaxFit = -1;
_Bool approachingSpeciationThreshold = 0;
_Bool keepIntroducingMutations = 1;
_Bool NO_MUTATIONS_DURING_ALLOPATRY = 0;
//int timeSampleCounter = 0;
_Bool RECORDING_TIMES_IN_FILE = RECORDING_TIMES_IN_FILE_DEFAULT;
int nRecordingTimes = -1;
long int *vectorOfRecordingTimes;
long int nextRecordingTime = -1;
int recordingTimesCompleted = 0;

FILE *effMigRates, *FSTtimeSeries, *grossMigRates, *mutationLog, *nVarTS, *DXYTS;
FILE *allelesLost, *logOfRemovedLoci, *LDfpt, *epistasisTS, *fixationLog;
FILE *fitTSdeme0, *fitTSdeme1, *AFtimeSeries, *selectedFrequencies, *neutralFrequencies;
FILE *LDselSitesAvg, *LDselSitesSame, *LDselSitesDiff, *AFSTS, *selAFSTS, *neutAFSTS;
FILE *LDneutSitesAvg, *LDneutSitesSame, *LDneutSitesDiff, *effPopSizeData;

#if (D == 1)
int n_in_each_patch[PATCHES];
int nPATCHES = PATCHES;
int BEST_ALLELE_IN_PATCH[PATCHES];
double epi_patch_multipliers[PATCHES];
int migrationCount[PATCHES];
#elif (D == 2)
int n_in_each_patch[(PATCHES * PATCHES)];
int nPATCHES = PATCHES * PATCHES;
int BEST_ALLELE_IN_PATCH[(PATCHES * PATCHES)];
double epi_patch_multipliers[(PATCHES * PATCHES)];
int migrationCount[(PATCHES * PATCHES)];
#endif


// core biology functions
void addALocus(double mapLoc, int locType);
void bagOfGenes(short int *ogtpt, int *noff, int newN);
int calculateAlleleFrequencies(void);
void calculateFitnesses(double *f, double *fsum);
void calculateLDpair(int l1, int l2, double dist, double *Dsum, double *DprimeSum, double *DeltaSum, int gatherLDvalues);
double calculateLDpairOneOff(int l1, int l2);
void calculateLD(int gatherLDvalues);
void calculateLDneutralSitesOnly(int gatherLDvalues);
void calculateLDselectedSitesOnly(int gatherLDvalues);
void makeZygoteChromosomes(int parent, short int *ogtpt);
void move(void);
void nextMutation(void);
void removeALocus(int locusNumber);
void reproduce(int gatherLDvalues);
void setUpFitnesses(void);
void setUpMap(void);
void setUpMutationSequence(void);


// utility functions
void adjustFitnesses(void);
void allocateGlobals(void);
double boxMuller(double mu, double sd);
void calcDXY(int nSamples); // Dxy using actual random samples
void calcDXY2(double *alleleFrequenciesByPatch); // expected DXY using allele frequencies
void calcExpectedME(double *fitpt, double *fitsumpt, int gatherLDvalues);
double calcMeanMagnitude(double *valuesArray, long int n);
double calcRandomHWFitness(int patchNum);
double calculateMaxPossibleFitness(int patchNum);
double calculateMinPossibleFitness(int patchNum);
double calcVariance(double *valuesArray, long int n);
void chooseDemeSample(int *sampleArray, int demeNumber, int nSamples);
//void collectDetailedTimeSample(void);
int compare_doubles(const void *a, const void *b);
int findNearestSelectedNeighbor(int focalLocus);
void growEpistasisMatrix(void);
void initializePopulation(void);
void loadCustomPopulationFromFiles(void);
void loadMapFromFile(void);
int maxInt(int i1, int i2);
int minInt(int i1, int i2);
void openDataFilesForRecording( int gatherLDvalues );
long int Poisson(double mm);
void printParameters(long int estRunLength);
int scalarPatchNumber(double x, double y);
int seed_gen(void);
long int setTSfreq(void);
void shrinkEpistasisMatrix(int locusToRemove);
void sortPopulation(void);
void sortPopulation2(void);
void testPrints(void);
void usage(char *s);
void warmUpRNG(void);
void printTime(long);

//
// pooled malloc() utility functions
//

#define PALLOC_TRACK_LIMIT 100	// number of allocation sizes to track
#define PALLOC_MINSIZE	(32 * 1024 * 1024) // minimum size for pooled allocation

struct palloc_pool_track_t {
    int	tag;
    size_t	size;
    int	requested;
    int	allocated;
    int	reallocated;
} palloc_pool_track[PALLOC_TRACK_LIMIT];

int palloc_track_count = 0;

typedef struct palloc_hdr {
    size_t	size;	// size of allocation (not counting palloc header)
    int	tag;	// tag supplied by caller
    void	*next;	// pointer to next free buffer on freelist
} palloc_hdr;

//
// Size of the header for allocated memory managed by the palloc system,
// rounded up to a nice boundary.
//
#define PALLOC_HDR_SIZE ((((int)sizeof(palloc_hdr) + 128) / 128) * 128)

//
// set PALLOC environment variable to:
//	1	prints stats about allocations
//	2	prints stats + each call to palloc()
//	3	stats + palloc() plus pfree()
//	4	stats + palloc() plus pfree() + reclaim info
//	5	stats + palloc() plus pfree() + reclaim info + freelist
//
int palloc_debug = 0;			// from PALLOC env var

palloc_hdr *palloc_free_list = NULL;	// linked list of free allocations
int palloc_count = 0;			// count of calls to palloc
int palloc_reclaimed = 0;		// count of reclaimed allocations
int palloc_realloced = 0;		// count of realloced allocations
int palloc_active = 0;			// number of active allocations
int palloc_max_active = 0;		// max number of active allocations
int palloc_free = 0;			// number of free allocations
int palloc_max_free = 0;		// max number of free allocations

void *palloc(int pool, size_t size);
void pfree(void *ptr);
void palloc_stats();

#define	POOL_PREVPATCH 1           // for previousPatches allocations
#define	POOL_GENOTYPES 2           // for genotypes allocations
#define	POOL_PATCH_LOCATIONS 3     // for patch_locations allocations
#define	POOL_X_LOCATIONS 4         // for x_locations allocations
#define	POOL_Y_LOCATIONS 5         // for y_locations allocations
#define	POOL_FITNESSES 6           // for fitness array allocations
#define   POOL_EPISTASIS 7           // for epistasisMatrix allocations
#define   POOL_EPI_COEFFS 8          // for epi_coeff_matrix allocations

int
main(int argc, char *argv[])
{
    
    // read in optional command line arguments ...
    int ch, foo = 0, nSamples = 200;
    char *progname = argv[0];
    _Bool keepTrying = 1, mWasIncremented;
    int gatherLDvalues = 0;
    struct timeval tod;
    unsigned long startTime, elapsedTime;
    int i, ii, jj, lostByDrift, lossCount;
    short int *sipt;
    long int dumli, count, estRunLength;
    _Bool keepGoing = 1, equilibriumReached = 0; // while loop variables
    long int mutationsThisTime, mtt;
    
    while ((ch = getopt(argc, argv, "A:a:B:b:C:c:D:d:E:e:F:f:G:g:H:h:I:i:K:L:l:M:m:N:n:O:P:p:R:rS:s:T:t:U:u:V:W:w:Z:?")) != -1) {
        switch (ch) {
            case 'A':
                END_PERIOD_ALLOPATRY = atoi(optarg);
                break;
            case 'a':
                START_PERIOD_ALLOPATRY = atoi(optarg);
                break;
            case 'B':
                NUMBER_BIG_TO_ESTABLISH = atoi(optarg);
                break;
            case 'b':
                START_WITH_BIG_IN_PLACE = atoi(optarg);
                break;
            case 'C':
                nCHROMOSOMES = atoi(optarg);
                break;
            case 'c':
                RECORD_FIT_TS = atoi(optarg);
                break;
            case 'D':
                DETERMINISTIC = atoi(optarg);
                break;
            case 'd':
                SD_MOVE = strtod(optarg, (char **)NULL);
                break;
            case 'E':
                DMI_MEAN_EFFECT = strtod(optarg, (char **)NULL);
                break;
            case 'e':
                TWO_DEME = atoi(optarg);
                break;
            case 'F':
                FIXED_N = atoi(optarg);
                break;
            case 'f':
                MUTATION_DISTRIBUTION = atoi(optarg);
                break;
            case 'G':
                GAMETE_PRODUCTION_MODE = atoi(optarg);
                break;
            case 'g':
                gatherLDvalues = atoi(optarg);
                break;
            case 'H':
                HVAL = strtod(optarg, (char **)NULL);
                break;
            case 'h':
                PURE_ALLOPATRY_MODEL = atoi(optarg);
                break;
            case 'I':
                DMI_PROBABILITY = strtod(optarg, (char **)NULL);
                break;
            case 'i':
                POSITIVE_EPI_PROBABILITY = strtod(optarg, (char **)NULL);
                break;
            case 'K':
                K = strtod(optarg, (char **)NULL);
                break;
            case 'L':
                OFFSPRING_IN_RANDOM_LOCATIONS = atoi(optarg);
                break;
            case 'l':
                TOTAL_MAP_LENGTH = strtod(optarg, (char **)NULL);
                break;
            case 'M':
                MAP_TYPE = atoi(optarg);
                break;
            case 'm':
                nMUTATIONS = atoi(optarg);
                break;
            case 'N':
                INITIAL_POPULATION_SIZE = atoi(optarg);
                break;
            case 'n':
                PURE_NEUTRAL_MODEL = atoi(optarg);
                break;
            case 'O':
                MOSAIC = atoi(optarg);
                break;
            case 'P':
                FATTEN_TAIL_PROPORTION = strtod(optarg, (char **)NULL);
                break;
            case 'p':
                MEAN_POSITIVE_EPI_COEFF = strtod(optarg, (char **)NULL);
                break;
            case 'R':
                RI_THRESH = strtod(optarg, (char **)NULL);
                break;
            case 'r':
                RECORDING_TIMES_IN_FILE = 1;
                break;
            case 'S':
                MEAN_S = strtod(optarg, (char **)NULL);
                break;
            case 's':
                SYMMETRIC_MUTATIONS = atoi(optarg);
                break;
            case 'T':
                nTIME_SAMPLES = atoi(optarg);
                break;
            case 't':
                nGENERATIONS_MAX = (long int) atoi(optarg);
                break;
            case 'U':
                MUTATIONS_PER_GENERATION = atoi(optarg);
                break;
            case 'u':
                FRACTION_SELECTED_LOCI = strtod(optarg, (char **)NULL);
                break;
            case 'V':
                RECORD_LD_VALUES = (_Bool) atoi(optarg);
                break;
            case 'W':
                DEME0_CONSTANT_S = strtod(optarg, (char **)NULL);
                break;
            case 'w':
                DEME1_CONSTANT_S = strtod(optarg, (char **)NULL);
                break;
            case 'Z':
                BIG_START_S_VAL = strtod(optarg, (char **)NULL);
                break;
            case '?':
            default:
                usage(progname);
                exit(-1);
        }
    }
    
    if ( gatherLDvalues >= 3 )
        RECORD_LD_VALUES = 1;
    
    estRunLength = setTSfreq();
    if ( PURE_NEUTRAL_MODEL ) {
        FRACTION_SELECTED_LOCI = 0.0;
        RECORD_FIT_TS = 0;
    }
    
    RI_THRESH = fmin( RI_THRESH, (INITIAL_POPULATION_SIZE * SD_MOVE / 1000) );
    
    if ( (PROBABILITY_DERIVED_REVERSED > 0.0) && (PROBABILITY_DERIVED_REVERSED < 1.0) ) {
        DMI_PROBABILITY = DMI_PROBABILITY / ( 2.0 * PROBABILITY_DERIVED_REVERSED * (1.0 - PROBABILITY_DERIVED_REVERSED) );
        POSITIVE_EPI_PROBABILITY = POSITIVE_EPI_PROBABILITY / ( (PROBABILITY_DERIVED_REVERSED * PROBABILITY_DERIVED_REVERSED) + ((1.0 - PROBABILITY_DERIVED_REVERSED) * (1.0 - PROBABILITY_DERIVED_REVERSED)) );
        // multiplied to account for the fact that many
        // of the pairs of loci are ineligible since
        // we require that one locus has a reversed
        // fitness scheme and one does not for DMI, or both are teh
        // same for + epistatic interaction
        // hence, prob{DMI} becomes prob / (2pq), and prob{+ epistatic} becomes prob / (p^2 + q^2)
    }
    
    if ( PROBABILITY_DERIVED_REVERSED > 0.0 && PROBABILITY_DERIVED_REVERSED < 1.0 && DMI_PROBABILITY > 0.0 )
        CONSIDER_EPISTASIS = 1;
    else if ( POSITIVE_EPI_PROBABILITY > 0.0 )
        CONSIDER_EPISTASIS = 1;
    else
        CONSIDER_EPISTASIS = 0;  // no need to consider epistatic interactions unless at least one of the above is true
    
    if ( START_WITH_BIG_IN_PLACE && (nPATCHES > 2) ) {
        fprintf(stderr,"\nHey shortsighted programmer: \n\t\tYour algorithm for starting big alleles\n\t\tat mutation-selection balance does not work\n\t\tfor this many patches!\n\n");
        exit(1);
    }
    
    // random number seed
    if (DETERMINISTIC) {
        int seed, rcount;
        FILE *fpt;
        fpt = fopen("RnumSeed.txt","r");
        if (fpt == NULL) {
            perror("Can't open RnumSeed.txt");
            exit(-1);
        }
        
        rcount = fscanf(fpt,"%i",&seed);
        if ( rcount ) {
            seedRand(seed);	// fixed random number seed
            fclose(fpt);
        }
        else {
            fprintf(stderr, "\n\n\tError! nothing read from file! Exiting!\n\n");
            exit(-1);
        }
        //fprintf(stderr, "\n\nSeed = %i\n\n",seed);
    }
    else {
        seedRand(seed_gen());	// get random number seed (system time)
    }
    warmUpRNG();			// necessary?  can't hurt.
    
    
    // get to work
    
    setUpMap();
    setUpMutationSequence();
    
    initializePopulation();
    setUpFitnesses();
    
    openDataFilesForRecording( gatherLDvalues );
    
    // record parameter values
    printParameters(estRunLength);
    
    
    if ( EARLY_TEST_KILL > 0 ) {
        if ( EARLY_TEST_KILL > 1 )
            testPrints();
        exit(0);
    }
    //
    // keep track of elapsed time for simulation.  gettimeofday()
    // used to watch wall time for the main part of the simulation,
    // ignoring setup and final data dump
    //
    gettimeofday(&tod, NULL);
    startTime = (long)tod.tv_sec * 1000000 + (long)tod.tv_usec;
    
    // now for the actual loop that does all the business
    
    for ( m = 0; m < nMUTATIONS; m++ ) {
        keepGoing = 1;
        t = 0;
        equilibriumReached = 0;
        lossCount = 0;
        
        while (keepGoing) {
            previousPatches = (int*) palloc(POOL_PREVPATCH,  (sizeof(int) * N) );
            
            // migration
            move();
            
            // a sort is needed for other algorithms after migration happens
            if ( !PURE_ALLOPATRY_MODEL )
                sortPopulation2();
            
            // introduce the next mutation(s)
            mWasIncremented = 0;
            if ( t == 0 && keepIntroducingMutations ) {
                for ( mtt = 0; mtt < MUTATIONS_PER_GENERATION; mtt++ ) {
                    if ( m < nMUTATIONS ) {
                        nextMutation(); // introduce a mutation
                        numberOfMutationTries++;
                        m++;
                        mWasIncremented = 1;
                    }
                }
                if ( mWasIncremented )
                    m--; // decrement by one since it also gets incremented one unit above in the for loop.
            }
            
            // replace parents with offspring
            reproduce(gatherLDvalues);
            
            // get some stats about allele frequencies and their changes
            lostByDrift = calculateAlleleFrequencies();
            
            // some bookkeeping at regular intervals
            if ( (totalGenerationsElapsed % TS_SAMPLING_FREQUENCY == 0 && totalGenerationsElapsed > 0 ) || ( (totalGenerationsElapsed+1) == (nMUTATIONS * nGENERATIONS_MAX) ) || RI_REACHED || ( RECORDING_TIMES_IN_FILE && totalGenerationsElapsed == nextRecordingTime ) ) {
                
                
                
                fprintf(nVarTS,"%li %i %i %i %li %li %li %li %li %li",totalGenerationsElapsed, nVariableLoci, nSELECTED_LOCI, nNEUTRAL_LOCI, totalDMIs, totalPositiveEpis, nPositiveEpisRegular, nPositiveEpisReversed, nRegularLoci, nReversedLoci);
                
                /* // testcheck
                 if ( (nPositiveEpisRegular + nPositiveEpisReversed) != totalPositiveEpis ) {
                 fprintf( stderr, "\n Positive Epi numbers don't add up!\n\treg = %li, rev = %li, total = %li\nExiting!\n", nPositiveEpisRegular, nPositiveEpisReversed, totalPositiveEpis);
                 exit(1);
                 } */
                
                for ( ii = 0; ii < nCHROMOSOMES; ii++ )
                    fprintf(nVarTS, " %i", (LOCI_PER_CHROMOSOME[ii] - 2));  // subtract 2 for the reference loci
                fprintf(nVarTS, "\n");
                
                //calcDXY(nSamples);
                
                if ( CONSIDER_EPISTASIS ) {
                    for ( ii = 0; ii < (nLOCI-1); ii++ ) {
                        sipt = epistasisMatrix + (ii * nLOCI) + ii + 1;
                        for ( jj = ii + 1; jj < nLOCI; jj++ ) {
                            if ( (*sipt) ) {  // there is an incompatibility
                                fprintf( epistasisTS,"%li %li %E %i %i %li %E %i %i %i %E\n", totalGenerationsElapsed, locusID[ii], MAP[ii], chromosomeMembership[ii], is_reversed_locus[ii], locusID[jj], MAP[jj], chromosomeMembership[jj], is_reversed_locus[jj], *sipt, (*(epi_coeff_matrix + (ii * nLOCI) + jj)) );
                            }
                            sipt++;
                        }
                    }
                }
                
                if ( RECORDING_TIMES_IN_FILE ) {
                    recordingTimesCompleted++;
                    if ( recordingTimesCompleted < nRecordingTimes ) {
                        nextRecordingTime = vectorOfRecordingTimes[recordingTimesCompleted];
                    }
                    else {
                        // all done; break the loops
                        keepGoing = 0;
                        m = nMUTATIONS + 1;
                        printParameters(estRunLength);
                    }
                }
                
            }
            
            if ( RI_REACHED ) {
                fprintf(stderr,"\nRI_THRESH reached.  Exiting simulation during generation %li\n", totalGenerationsElapsed);
                printParameters(estRunLength);
                keepGoing = 0; // break the while loop
                m = nMUTATIONS + 1; // break the for loop
            }
            
            pfree(previousPatches);
            
            t++;
            totalGenerationsElapsed++;
            
            if ( START_PERIOD_ALLOPATRY > 0 ) {
                if ( totalGenerationsElapsed == START_PERIOD_ALLOPATRY ) {
                    END_PERIOD_ALLOPATRY = START_PERIOD_ALLOPATRY + PERIOD_ALLOPATRY;
                    if ( NO_MUTATIONS_DURING_ALLOPATRY )
                        keepIntroducingMutations = 0;
                }
            }
            if ( totalGenerationsElapsed == (END_PERIOD_ALLOPATRY))
                keepIntroducingMutations = 1;
            
            if ( equilibriumReached || t >= nGENERATIONS_MAX ) {
                keepGoing = 0;
            }
        }
        
        
    }
    
    
    //
    // Report run time and estimate of time need to do a full
    // run on target test machine.  This is based on a ratio of
    // the current number of timesteps and generations to a "full" run.
    //
    gettimeofday(&tod, NULL);
    elapsedTime = (long)tod.tv_sec*1000000 + (long)tod.tv_usec - startTime;
    printf("\nElapsed time:  ");
    printTime(elapsedTime);
    
    printf("Estimated time for full run:  ");
    printTime((long)elapsedTime *
              ((double)nGENERATIONS_MAX_DEFAULT/(double)nGENERATIONS_MAX) *
              ((double)nMUTATIONS_DEFAULT / (double)nMUTATIONS));
    
    fprintf(stderr, "\nEst. run length = %li, actual = %li, diff = %li, obs./exp ratio = %f\n", estRunLength, totalGenerationsElapsed, (estRunLength - totalGenerationsElapsed), ((double) totalGenerationsElapsed)/((double) estRunLength));
    
    testPrints();
    if ( !RI_REACHED )
        printParameters(estRunLength);
    
    fclose(effMigRates);
    fclose(FSTtimeSeries);
    fclose(AFtimeSeries);
    fclose(grossMigRates);
    fclose(mutationLog);
    fclose(nVarTS);
    fclose(logOfRemovedLoci);
    fclose(fixationLog);
    fclose(selectedFrequencies);
    fclose(neutralFrequencies);
    if ( CONSIDER_EPISTASIS ) {
        fclose(epistasisTS);
        pfree(epistasisMatrix);
        pfree(epi_coeff_matrix);
    }
    
    if ( RECORD_FIT_TS ) {
        fclose(fitTSdeme0);
        fclose(fitTSdeme1);
    }
    fclose(DXYTS);
    if ( gatherLDvalues >= 1 ) {
        fclose(LDselSitesAvg);
        fclose(LDneutSitesAvg);
    }
    if ( gatherLDvalues >= 3 ) {
        fclose(LDselSitesSame);
        fclose(LDselSitesDiff);
        fclose(LDneutSitesDiff);
        fclose(LDneutSitesSame);
        fclose(LDfpt);
    }
    fclose(effPopSizeData);
    fclose(AFSTS);
    fclose(selAFSTS);
    fclose(neutAFSTS);
    
    
    pfree(x_locations);
#if (D == 2)
    pfree(y_locations);
#endif
    pfree(genotypes);
    pfree(patch_locations);
    
    
    free(fit0);
    free(fit1);
    free(fitH);
    free(MUTATION_ORDER);
    free(MUTATION_TYPE_SEQUENCE);
    free(MUTATION_LOCATIONS);
    free(LOCI_PER_CHROMOSOME);
    free(MAP_LENGTHS);
    free(S_MAX1_SEQUENCE);
    free(S_MAX0_SEQUENCE);
    free(locusID);
    free(vectorOfRecordingTimes);
    
    palloc_stats();
    return 0;
}


/*  ***************************************************** */
/*  *************  CORE BIOLOGY FUNCTIONS  ************** */
/*  ***************************************************** */


void addALocus(double mapLoc, int locType)
{
    int i, j, focalChrom, locus = 0, ii;
    _Bool looking = 1;
    double positionOnFocalChrom = mapLoc, totLength = 0.0;
    
    
    // find where the new locus goes
    i = 0;
    while ( looking ) {
        if ( i >= nCHROMOSOMES ) {
            fprintf(stderr,"\n\nERROR!  Failed to find correct spot for new locus!\ni = %i, nCHROMOSOMES = %i, positionOnFocalChrom = %f, t = %li, mapLoc = %f, totLength = %f\n", i, nCHROMOSOMES, positionOnFocalChrom, totalGenerationsElapsed, mapLoc, totLength);
            exit(1);
        }
        
        totLength += MAP_LENGTHS[i];
        
        if ( mapLoc <= totLength ) {
            focalChrom = i;
            looking = 0;
            while ( (positionOnFocalChrom > MAP[locus]) && locus < nLOCI ) {
                locus++;
            }
            newestLocus = locus; // the point from which all arrays dealing with genotypes will have to shift
        }
        else {
            positionOnFocalChrom = positionOnFocalChrom - MAP_LENGTHS[i];
            if ( positionOnFocalChrom < 0 ) {
                fprintf(stderr,"\n\nERROR!  Add locus algorithm produced negative position = %E\n\n", positionOnFocalChrom);
                exit(1);
            }
            locus += LOCI_PER_CHROMOSOME[i];
        }
        i++;
    }
    if ( (MUTATION_DISTRIBUTION == 2) && (nVariableLoci < NUMBER_BIG_TO_ESTABLISH) && (!START_WITH_BIG_IN_PLACE) ) {
        S_MAX1_SEQUENCE[m] = randU() * FATTEN_TAIL_MAX;
        locType = 1; // override original layout
        MUTATION_TYPE_SEQUENCE[m] = 1; // override original layout
        if ( SYMMETRIC_MUTATIONS )
            S_MAX0_SEQUENCE[m] = S_MAX1_SEQUENCE[m];
        else
            S_MAX0_SEQUENCE[m] = randU() * FATTEN_TAIL_MAX;
    }
    else if ( MUTATION_DISTRIBUTION == 3 && locType ) {
        S_MAX1_SEQUENCE[m] = DEME1_CONSTANT_S;
        if ( SYMMETRIC_MUTATIONS )
            S_MAX0_SEQUENCE[m] = DEME1_CONSTANT_S;
        else
            S_MAX0_SEQUENCE[m] = DEME0_CONSTANT_S;
    }
    fprintf(mutationLog,"%li %li %E %i %E %E %E %E %i %i %i\n", totalGenerationsElapsed, m, mapLoc, focalChrom, positionOnFocalChrom, S_MAX0_SEQUENCE[m], S_MAX1_SEQUENCE[m], H_SEQUENCE[m], MUTATION_LOCATIONS[m], REVERSAL_SEQUENCE[m], locType);
    // test prints
    // the map
    /*
     fprintf(stderr,"\nThe old map:\n");
     for ( i = 0; i < nLOCI; i++ )
     fprintf(stderr,"%f ", MAP[i]);
     fprintf(stderr, "\nThe old number of LOCI_PER_CHROMOSOME:\n");
     for (i=0; i<nCHROMOSOMES; i++) {
     fprintf(stderr, "%i  ", LOCI_PER_CHROMOSOME[i]);
     }
     fprintf(stderr,"\n\nThe lengths of chromosomes in the map of total length %f:\n",TOTAL_MAP_LENGTH);
     for ( i = 0; i < nCHROMOSOMES; i++ )
     fprintf(stderr,"%f ", MAP_LENGTHS[i]);
     */
    // end test prints
    
    // adjust global variables accordingly by resizing and copying arrays
    nLOCI++;
    double *newAlleleFreq, *newMAP, *newSMAX1, *newSMAX0, *newH;
    int *newISSELLOCUS, *newChromMembership;
    int *newVariableLoci, *newIsRefLocus;
    short int *newFixedAllele, *newGtpt, *oldGtpt;
    short int *newGenotypes;
    int nIndexesBefore = newestLocus; // {0,1,2,3,...,(newestLocus - 1) }
    int nIndexesAfter;
    double *newFit0, *newFit1, *newFitH;
    long int *newLocusID;
    _Bool *newIsReversedLocus, *newAnyEpisAtLocus;
    
    nIndexesAfter = (nLOCI - 1) - newestLocus;
    
    newAlleleFreq = (double *) malloc( sizeof(double) * nLOCI );
    newMAP = (double *) malloc( sizeof(double) * nLOCI );
    newSMAX0 = (double *) malloc( sizeof(double) * nLOCI );
    newSMAX1 = (double *) malloc( sizeof(double) * nLOCI );
    newH = (double *) malloc( sizeof(double) * nLOCI );
    newISSELLOCUS = (int *) malloc( sizeof(int) * nLOCI );
    newChromMembership = (int *) malloc( sizeof(int) * nLOCI );
    newVariableLoci = (int *) malloc( sizeof(int) * nLOCI );
    newIsRefLocus = (int *) malloc( sizeof(int) * nLOCI );
    newFixedAllele = (short int *) malloc( sizeof(short int) * nLOCI);
    newGenotypes = (short int *) palloc( POOL_GENOTYPES, sizeof(short int) * N * nLOCI * 2 );
    newLocusID = (long int *) malloc( sizeof(long int) * nLOCI );
    newIsReversedLocus = (_Bool *) malloc( sizeof(_Bool) * nLOCI );
    if ( CONSIDER_EPISTASIS )
        newAnyEpisAtLocus = (_Bool *) malloc( sizeof(_Bool) * nLOCI );
    
    if ( nIndexesBefore > 0 ) {
        memcpy( newAlleleFreq, allele_frequencies, (sizeof(double) * nIndexesBefore) );
        memcpy( newMAP, MAP, (sizeof(double) * nIndexesBefore) );
        memcpy( newSMAX0, S_MAX0, (sizeof(double) * nIndexesBefore) );
        memcpy( newSMAX1, S_MAX1, (sizeof(double) * nIndexesBefore) );
        memcpy( newH, H, (sizeof(double) * nIndexesBefore) );
        memcpy( newISSELLOCUS, IS_SELECTED_LOCUS, (sizeof(int) * nIndexesBefore) );
        memcpy( newChromMembership, chromosomeMembership, (sizeof(int) * nIndexesBefore) );
        memcpy( newVariableLoci, variable_loci, (sizeof(int) * nIndexesBefore) );
        memcpy( newIsRefLocus, is_reference_locus, (sizeof(int) * nIndexesBefore) );
        memcpy( newFixedAllele, fixed_allele, (sizeof(short int) * nIndexesBefore) );
        memcpy( newLocusID, locusID, (sizeof(long int) * nIndexesBefore) );
        memcpy( newIsReversedLocus, is_reversed_locus, (sizeof(_Bool) * nIndexesBefore) );
        if ( CONSIDER_EPISTASIS )
            memcpy( newAnyEpisAtLocus, any_epis_for_this_locus, (sizeof(_Bool) * nIndexesBefore) );
    }
    
    newAlleleFreq[newestLocus] = 0.0;
    newMAP[newestLocus] = positionOnFocalChrom;
    newSMAX0[newestLocus] = S_MAX0_SEQUENCE[m];
    newSMAX1[newestLocus] = S_MAX1_SEQUENCE[m];
    newH[newestLocus] = H_SEQUENCE[m];
    newISSELLOCUS[newestLocus] = locType;
    if ( locType ) {
        nSELECTED_LOCI++;
    }
    else {
        nNEUTRAL_LOCI++;
    }
    newChromMembership[newestLocus] = focalChrom;
    newVariableLoci[newestLocus] = 1;
    newIsRefLocus[newestLocus] = 0; // added loci are NOT reference loci
    newFixedAllele[newestLocus] = -1;
    newLocusID[newestLocus] = m; // use mutation counter as the locus ID
    
    newIsReversedLocus[newestLocus] = REVERSAL_SEQUENCE[m];
    if ( (newIsReversedLocus[newestLocus]) )
        nReversedLoci++;
    else
        nRegularLoci++;
    
    if ( CONSIDER_EPISTASIS )
        newAnyEpisAtLocus[newestLocus] = 0;  // will be altered if needed in GrowEpistasisMatrix function
    
    if ( nIndexesAfter > 0 ) {
        memcpy( (newAlleleFreq + newestLocus + 1), (allele_frequencies + newestLocus), (sizeof(double) * nIndexesAfter) );
        memcpy( (newMAP + newestLocus + 1), (MAP + newestLocus), (sizeof(double) * nIndexesAfter) );
        memcpy( (newSMAX0 + newestLocus + 1), (S_MAX0 + newestLocus), (sizeof(double) * nIndexesAfter) );
        memcpy( (newSMAX1 + newestLocus + 1), (S_MAX1 + newestLocus), (sizeof(double) * nIndexesAfter) );
        memcpy( (newH + newestLocus + 1), (H + newestLocus), (sizeof(double) * nIndexesAfter) );
        memcpy( (newISSELLOCUS + newestLocus + 1), (IS_SELECTED_LOCUS + newestLocus), (sizeof(int) * nIndexesAfter) );
        memcpy( (newChromMembership + newestLocus + 1), (chromosomeMembership + newestLocus), (sizeof(int) * nIndexesAfter) );
        memcpy( (newVariableLoci + newestLocus + 1), (variable_loci + newestLocus), (sizeof(int) * nIndexesAfter) );
        memcpy( (newIsRefLocus + newestLocus + 1), (is_reference_locus + newestLocus), (sizeof(int) * nIndexesAfter) );
        memcpy( (newFixedAllele + newestLocus + 1), (fixed_allele + newestLocus), (sizeof(short int) * nIndexesAfter) );
        memcpy( (newLocusID + newestLocus + 1), (locusID + newestLocus), (sizeof(long int) * nIndexesAfter) );
        memcpy( (newIsReversedLocus + newestLocus + 1), (is_reversed_locus + newestLocus), (sizeof(_Bool) * nIndexesAfter) );
        if ( CONSIDER_EPISTASIS )
            memcpy( (newAnyEpisAtLocus + newestLocus + 1), (any_epis_for_this_locus + newestLocus), (sizeof(_Bool) * nIndexesAfter) );
    }
    
    // copy and adjust genotypes!
    newGtpt = newGenotypes;
    oldGtpt = genotypes;
    for ( i = 0; i < N; i++ ) {
        
        if ( nIndexesBefore > 0 )
            memcpy( newGtpt, oldGtpt, (sizeof(short int) * nIndexesBefore * 2) );
        
        *(newGtpt + (newestLocus*2)) = 0;
        *(newGtpt + (newestLocus*2) + 1) = 0;
        
        if ( nIndexesAfter > 0 )
            memcpy( (newGtpt + (2 * (newestLocus+1))), (oldGtpt + (2 * newestLocus)), (sizeof(short int) * 2 * nIndexesAfter));
        
        newGtpt += (2 * nLOCI);
        oldGtpt += (2 * (nLOCI - 1));
    }
    
    pfree(genotypes);
    genotypes = newGenotypes;
    free(allele_frequencies);
    allele_frequencies = newAlleleFreq;
    free(MAP);
    MAP = newMAP;
    free(S_MAX0);
    S_MAX0 = newSMAX0;
    free(S_MAX1);
    S_MAX1 = newSMAX1;
    free(H);
    H = newH;
    free(IS_SELECTED_LOCUS);
    IS_SELECTED_LOCUS = newISSELLOCUS;
    free(chromosomeMembership);
    chromosomeMembership = newChromMembership;
    free(variable_loci);
    variable_loci = newVariableLoci;
    free(is_reversed_locus);
    is_reversed_locus = newIsReversedLocus;
    
    free(is_reference_locus);
    is_reference_locus = newIsRefLocus;
    free(fixed_allele);
    fixed_allele = newFixedAllele;
    free(locusID);
    locusID = newLocusID;
    
    LOCI_PER_CHROMOSOME[focalChrom] = LOCI_PER_CHROMOSOME[focalChrom] + 1;
    newFit0 = (double *) malloc( (nLOCI * nPATCHES * sizeof(double)) );
    newFit1 = (double *) malloc( (nLOCI * nPATCHES * sizeof(double)) );
    newFitH = (double *) malloc( (nLOCI * nPATCHES * sizeof(double)) );
    
    // all the loci before the newest one
    if ( nIndexesBefore > 0 ) {
        memcpy( newFit0, fit0, (sizeof(double) * nIndexesBefore * nPATCHES) );
        memcpy( newFit1, fit1, (sizeof(double) * nIndexesBefore * nPATCHES) );
        memcpy( newFitH, fitH, (sizeof(double) * nIndexesBefore * nPATCHES) );
    }
    
    // all the loci after the newest one
    if ( nIndexesAfter > 0 ) {
        memcpy( (newFit0 + (nPATCHES * (newestLocus + 1))), (fit0 + (nPATCHES * newestLocus)), (sizeof(double) * nIndexesAfter * nPATCHES) );
        memcpy( (newFit1 + (nPATCHES * (newestLocus + 1))), (fit1 + (nPATCHES * newestLocus)), (sizeof(double) * nIndexesAfter * nPATCHES) );
        memcpy( (newFitH + (nPATCHES * (newestLocus + 1))), (fitH + (nPATCHES * newestLocus)), (sizeof(double) * nIndexesAfter * nPATCHES) );
    }
    
    // switch the array assignments
    free(fitH);
    fitH = newFitH;
    free(fit0);
    fit0 = newFit0;
    free(fit1);
    fit1 = newFit1;
    
    // now put in the correct fitness values for the new locus
    adjustFitnesses();
    
    if ( CONSIDER_EPISTASIS ) {
        free(any_epis_for_this_locus);
        any_epis_for_this_locus = newAnyEpisAtLocus;
        growEpistasisMatrix();
    }
    
    // */
    
    
    
    
    /*
     // test prints
     // the map
     fprintf(stderr,"\nThe NEW map:\n");
     for ( i = 0; i < nLOCI; i++ )
     fprintf(stderr,"%.2f  ", MAP[i]);
     fprintf(stderr, "\n\nThe NEW number of LOCI_PER_CHROMOSOME:\n");
     for (i=0; i<nCHROMOSOMES; i++) {
     fprintf(stderr, "%i  ", LOCI_PER_CHROMOSOME[i]);
     }
     
     fprintf(stderr, "\n\nThe NEW is_reference_locus:\n");
     for (i=0; i<nLOCI; i++) {
     fprintf(stderr, "%i ", is_reference_locus[i]);
     }
     
     //fprintf(stderr,"\n\nThe lengths of chromosomes in the map of total length %f:\n",TOTAL_MAP_LENGTH);
     //	for ( i = 0; i < nCHROMOSOMES; i++ )
     //		fprintf(stderr,"%f ", MAP_LENGTHS[i]);
     
     fprintf(stderr,"\n\nThe new location in the map:\t%f",mapLoc);
     fprintf(stderr,"\nwill be at position:\t\t%f",positionOnFocalChrom);
     fprintf(stderr,"\non chromosome\t\t\t%i",focalChrom);
     fprintf(stderr,"\nThe index of the locus will be:\t#%li (%li of %li)\n",newestLocus,newestLocus+1,nLOCI);
     
     fprintf(stderr,"\nLocus IDs:\n");
     for ( i = 0; i < nLOCI; i++ )
     fprintf(stderr,"%li ", locusID[i]);
     fprintf(stderr,"\n");
     
     fprintf(stderr,"\nChromosome Membership\n");
     for ( i = 0; i < nLOCI; i++ )
     fprintf(stderr,"%i ", chromosomeMembership[i]);
     fprintf(stderr,"\n");
     
     fprintf(stderr,"\nVariable Loci:\n");
     for ( i = 0; i < nLOCI; i++ )
     fprintf(stderr,"%i ", variable_loci[i]);
     fprintf(stderr,"\n");
     for ( i = 0; i < nLOCI; i++ )
     if ( variable_loci[i] )
     fprintf(stderr,"%i ", i);
     fprintf(stderr,"\n");
     fprintf(stderr,"nVariableLoci = %i\n",nVariableLoci);
     
     fprintf(stderr,"\nfit0\t\t\tfit1\t\t\tfitH\n");
     newFit0 = fit0;
     newFit1 = fit1;
     newFitH = fitH;
     for ( i = 0; i < nLOCI; i++ ) {
     for ( j = 0; j < nPATCHES; j++ ) {
     fprintf(stderr,"%.3f ", (*newFit0));
     newFit0++;
     }
     fprintf(stderr,"\t\t");
     for ( j = 0; j < nPATCHES; j++ ) {
     fprintf(stderr,"%.3f ", (*newFit1));
     newFit1++;
     }
     fprintf(stderr,"\t\t");
     for ( j = 0; j < nPATCHES; j++ ) {
     fprintf(stderr,"%.3f ", (*newFitH));
     newFitH++;
     }
     fprintf(stderr,"\n");
     }
     
     
     fprintf(stderr,"\n*********************************************\n");
     
     
     // end of test prints
     //  */
}




void bagOfGenes(short int *ogtpt, int *noff, int newN)
{
    int i, j, l, *ipt = noff, patch, nborn, locus;
    short int *opt = ogtpt, *stpt, *spt;
    double gtsum, patchAlleleFrequencies[nPATCHES][nLOCI], Ninv;
    double p, q, prHet, prHomo, weightedFitVal, wp, dum;
    long int fitArrayIndex;
    
    // calculate allele frequencies in each patch
    for ( i = 0; i < nPATCHES; i++ ) {
        for ( j = 0; j < nLOCI; j++ ) {
            patchAlleleFrequencies[i][j] = 0.0;
        }
    }
    for ( i = 0; i < N; i++ ) {
        patch = patch_locations[i];
        stpt = genotypes + ( 2 * i * nLOCI );
        for ( locus = 0; locus < nLOCI; locus++ ) {
            if ( variable_loci[locus] ) {
                spt = stpt + ( 2 * locus );
                gtsum = ((double) ((*spt) + (*(spt+1))));
                patchAlleleFrequencies[patch][locus] = patchAlleleFrequencies[patch][locus] + gtsum;
            }
        }
    }
    for ( i = 0; i < nPATCHES; i++ ) {
        Ninv = 0.5 / ((double) n_in_each_patch[i]);
        for ( locus = 0; locus < nLOCI; locus++ ) {
            if ( variable_loci[locus] )
                patchAlleleFrequencies[i][locus] = patchAlleleFrequencies[i][locus] * Ninv;
            //else   // this else clause is unncessary
            //	patchAlleleFrequencies[i][locus] = (double) fixed_allele[locus];
            /* if 1 is fixed, then the frequency is 1.0
             if 0 is fixed then then the frequency is 0.0 */
        }
    }
    
    // set all alleles to zero by default, then change only those necessary
    for ( i = 0; i < ( 2 * nLOCI * newN ); i++ ) {
        *opt = 0;
        opt++;
    }
    
    stpt = ogtpt; // points to first allele, first locus, first offspring, first patch
    
    
    // change some to ones
    for ( i = 0; i < nPATCHES; i++ ) {
        nborn = noff[i];
        for ( locus = 0; locus < nLOCI; locus++ ) {
            if ( variable_loci[locus] ) {
                fitArrayIndex = (locus * nPATCHES) + i;
                p = patchAlleleFrequencies[i][locus];
                opt = stpt + ( 2 * locus ); /* stpt points to first allele of first locus of first offspring in THIS PATCH
                                             opt points to first allele of focal locus of first offspring in this patch */
                spt = opt + 1;				/* spt points to second allele of focal locus of first offspring in this patch */
                if ( p > 0.0 && p < 1.0 ) {
                    q = 1.0 - p;
                    if ( PURE_NEUTRAL_MODEL ) {
                        prHomo = p * p;
                        prHet = prHomo + (2.0 * p * q);
                    }
                    else {
                        if ( is_reversed_locus[locus] )
                            weightedFitVal = ( p * q * (*(fitH + fitArrayIndex)) ) + ( p * p * (*(fit0 + fitArrayIndex)) );
                        else
                            weightedFitVal = ( p * q * (*(fitH + fitArrayIndex)) ) + ( p * p * (*(fit1 + fitArrayIndex)) ); /* H-W assumption;
                                                                                                                             no "2" on first term because only half
                                                                                                                             alleles count */
                        if ( is_reversed_locus[locus] )
                            wp = weightedFitVal / ( ( 2.0 * p * q * (*(fitH + fitArrayIndex)) ) + ( p * p * (*(fit0 + fitArrayIndex)) ) + ( q * q * (*(fit1 + fitArrayIndex)) ));
                        else
                            wp = weightedFitVal / ( ( 2.0 * p * q * (*(fitH + fitArrayIndex)) ) + ( p * p * (*(fit1 + fitArrayIndex)) ) + ( q * q * (*(fit0 + fitArrayIndex)) )); /* weighted probability of reproduction involving p */
                        prHomo = wp * wp;
                        prHet = prHomo + (2.0 * wp * (1.0 - wp)); // cumulative probability
                    }
                    for ( j = 0; j < nborn; j++ ) {
                        dum = randU();
                        if ( dum < prHomo ) {
                            *opt = 1;
                            *spt = 1;
                        }
                        else if ( dum < prHet ) {
                            if ( randU() < 0.5 )
                                *opt = 1;
                            else
                                *spt = 1; // have to have this randomization so that LD doesn't get artificially elevated!
                        }
                        opt += ( 2 * nLOCI ); // advance to first allele of same locus in next offspring
                        spt += ( 2 * nLOCI ); // advance to second allele of same locus in next offspring
                    }
                }
                else if ( p >= 1.0 ) {
                    for ( j = 0; j < nborn; j++ ) {
                        *opt = 1;
                        *spt = 1;
                        opt += ( 2 * nLOCI ); // advance to first allele of same locus in next offspring
                        spt += ( 2 * nLOCI ); // advance to second allele of same locus in next offspring
                    }
                }
            }
            else if ( fixed_allele[locus] ) { // not a variable locus; 1 is fixed
                opt = stpt + ( 2 * locus ); /* stpt points to first allele of first locus of first offspring in THIS PATCH
                                             opt points to first allele of focal locus of first offspring in this patch */
                spt = opt + 1;				/* spt points to second allele of focal locus of first offspring in this patch */
                for ( j = 0; j < nborn; j++ ) {
                    *opt = 1;
                    *spt = 1;
                    opt += ( 2 * nLOCI ); // advance to first allele of same locus in next offspring
                    spt += ( 2 * nLOCI ); // advance to second allele of same locus in next offspring
                }
            }
            
        }
        stpt += ( 2 * nLOCI * nborn );
    }
}



int calculateAlleleFrequencies(void)
{
    int i, j, patch, locus, locType;
    int allelesByPatch[nLOCI][nPATCHES], totalAlleles[nLOCI], gtsum, taln;
    double Ninv, patchWeights[nPATCHES];
    double alleleFrequenciesByPatch[nLOCI][nPATCHES], AFdiff;
    double HSsum, FST[nLOCI], HT, p, Nsinv[nPATCHES], nsnDist, nsnS, nsnLD;
    short int *spt, *stpt;
    int nRemovals = 0, lociToRemove[nLOCI], nsn;
    int dim0, tempAlleleCount;
    dim0 = (2 * N) + 1;
    int AFS[dim0], selectedAFS[dim0], neutralAFS[dim0];
    FILE *JAFSTS;
    char fname[80];
    
    for ( i = 0; i < dim0; i++ ) {
        AFS[i] = 0;
        selectedAFS[i] = 0;
        neutralAFS[i] = 0;
    }
    
    /* allele frequencies have been set to zero initially for all loci;
     we only actually  need to check frequencies at loci where mutations have occurred. */
    
    for ( i = 0; i < nLOCI; i++ ) {
        totalAlleles[i] = 0;
        FST[i] = 0.0;
        for ( j = 0; j < nPATCHES; j++ ) {
            allelesByPatch[i][j] = 0;
        }
    }
    
    Ninv = 1.0 / ((double) N);
    for ( i = 0; i < nPATCHES; i++ ) {
        patchWeights[i] = ((double) n_in_each_patch[i]) * Ninv;
        Nsinv[i] = 0.5 / ((double) n_in_each_patch[i]); // 0.5 since each has two alleles for each locus; reduces calc's below
    }
    
    for ( i = 0; i < N; i++ ) {
        patch = patch_locations[i];
        stpt = genotypes + ( 2 * i * nLOCI );
        for ( locus = 0; locus < nLOCI; locus++ ) {
            if ( variable_loci[locus] ) {
                spt = stpt + ( 2 * locus );
                gtsum = (*spt) + (*(spt+1));
                totalAlleles[locus] = totalAlleles[locus] + gtsum;
                allelesByPatch[locus][patch] = allelesByPatch[locus][patch] + gtsum;
            }
        }
    }
    
    if ( (totalGenerationsElapsed % TS_SAMPLING_FREQUENCY == 0 && totalGenerationsElapsed > 0 ) || ( (totalGenerationsElapsed+1) == (nMUTATIONS * nGENERATIONS_MAX) ) || RI_REACHED || ( RECORDING_TIMES_IN_FILE && totalGenerationsElapsed == nextRecordingTime ) ) {
        // site frequency spectrum:
        for ( locus = 0; locus < nLOCI; locus++ ) {
            if ( !is_reference_locus[locus] ) {
                tempAlleleCount = totalAlleles[locus];
                AFS[tempAlleleCount] = AFS[tempAlleleCount] + 1;
                if ( IS_SELECTED_LOCUS[locus] )
                    selectedAFS[tempAlleleCount] = selectedAFS[tempAlleleCount] + 1;
                else
                    neutralAFS[tempAlleleCount] = neutralAFS[tempAlleleCount] + 1;
            }
        }
        for ( i = 0; i < dim0; i++ ) {
            if ( AFS[i] > 0 )
                fprintf( AFSTS, "%li,%i,%i\n", totalGenerationsElapsed, i, AFS[i] );
            if ( selectedAFS[i] > 0 )
                fprintf( selAFSTS, "%li,%i,%i\n", totalGenerationsElapsed, i, selectedAFS[i] );
            if ( neutralAFS[i] > 0 )
                fprintf( neutAFSTS, "%li,%i,%i\n", totalGenerationsElapsed, i, neutralAFS[i] );
        }
        
        // raw data needed for joint site frequency spectrum
        sprintf(fname, "JAFSdata/JAFSdata%li.csv", totalGenerationsElapsed);
        JAFSTS = fopen(fname, "w");
        // header row:
        fprintf( JAFSTS, "locusID,IS_SELECTED_LOCUS" );
        for ( i = 0; i < nPATCHES; i++ )
            fprintf(JAFSTS, ",CountPatch%i", i);
        fprintf(JAFSTS, "\n");
        // data:
        for ( locus = 0; locus < nLOCI; locus++ ) {
            if ( !is_reference_locus[locus] ) {
                fprintf( JAFSTS, "%li,%i", locusID[locus], IS_SELECTED_LOCUS[locus] );
                for ( i = 0; i < nPATCHES; i++ )
                    fprintf( JAFSTS, ",%i", allelesByPatch[locus][i] );
                fprintf( JAFSTS, "\n" );
            }
        }
        fclose(JAFSTS);
        // number of individuals in each patch data for JAFS:
        // as .R file:
        sprintf(fname, "JAFSdata/JAFSnumPerPatch%li.R", totalGenerationsElapsed);
        JAFSTS = fopen(fname, "w");
        fprintf( JAFSTS, "numPerPatch%li <- c(", totalGenerationsElapsed);
        for ( i = 0; i < nPATCHES; i++ ) {
            if ( i == 0 )
                fprintf( JAFSTS, "%i", n_in_each_patch[i] );
            else
                fprintf( JAFSTS, ",%i", n_in_each_patch[i] );
        }
        fprintf( JAFSTS, ")\n" );
        fclose(JAFSTS);
        // as .m file:
        sprintf(fname, "JAFSdata/JAFSnumPerPatch%li.m", totalGenerationsElapsed);
        JAFSTS = fopen(fname, "w");
        fprintf( JAFSTS, "numPerPatch%li = [", totalGenerationsElapsed);
        for ( i = 0; i < nPATCHES; i++ ) {
            if ( i == 0 )
                fprintf( JAFSTS, "%i", n_in_each_patch[i] );
            else
                fprintf( JAFSTS, ",%i", n_in_each_patch[i] );
        }
        fprintf( JAFSTS, "];\n" );
        fclose(JAFSTS);
    }
    
    
    
    // now calculate allele frequencies and FST
    nVariableLoci = 0;
    for ( locus = 0; locus < nLOCI; locus++ ) {
        if ( variable_loci[locus] ) {
            p = ((double) totalAlleles[locus]) * Ninv * 0.5;
            allele_frequencies[locus] = p;
            if ( p <= 0.0 || p >= 1.0 ) { // inequalities here to account for possibility of imprecision in discrete math
                variable_loci[locus] = 0;
                if ( p >= 1.0 ) {
                    fixed_allele[locus] = 1;
                    if ( !is_reference_locus[locus] ) {
                        lociToRemove[nRemovals] = locus;
                        nRemovals++;
                    }
                    fprintf( fixationLog, "%li %li %E %i %E %i %i\n", totalGenerationsElapsed, locusID[locus], MAP[locus], chromosomeMembership[locus], S_MAX1[locus], is_reversed_locus[locus], IS_SELECTED_LOCUS[locus] );
                }
                else {
                    fixed_allele[locus] = 0;
                    if ( !is_reference_locus[locus] ) {
                        lociToRemove[nRemovals] = locus;
                        nRemovals++;
                    }
                }
            }
            else
                nVariableLoci++;
            
            if ( ((totalGenerationsElapsed % TS_SAMPLING_FREQUENCY) == 0 && totalGenerationsElapsed > 0 ) || ( (totalGenerationsElapsed+1) == (nMUTATIONS * nGENERATIONS_MAX) ) || ( RECORDING_TIMES_IN_FILE && totalGenerationsElapsed == nextRecordingTime ) ) {
                if ( p <= 0.0 || p >= 1.0 )
                    FST[locus] = 0.0;
                else {
                    HT = 2.0 * p * (1.0 - p);
                    HSsum = 0.0;
                    for ( j = 0; j < nPATCHES; j++ ) {
                        p = ((double) allelesByPatch[locus][j]) * Nsinv[j];
                        alleleFrequenciesByPatch[locus][j] = p;
                        // /* // test check
                        if ( p < 0.0 || p > 1.0 ) {
                            fprintf(stderr, "\nERROR!, p = %E\n",p);
                        }
                        // */ // end test check
                        HSsum += ( 2.0 * p * (1.0 - p) ) * patchWeights[j];
                    }
                    FST[locus] = ( HT - HSsum ) / HT;
                    // /* // test check
                    if ( FST[locus] > 1.0 || FST[locus] < -1.0E-13 ) {
                        fprintf(stderr, "\n\nFST ERROR!!\nFST[%i] = %E, HT = %E, HSsum = %E\n", locus, FST[locus],HT, HSsum);
                        exit(1);
                    }
                    if ( FST[locus] < 0 )
                        FST[locus] = 0;
                    // */ // end test check
                }
            }
        }
    }
    if ( ((totalGenerationsElapsed % TS_SAMPLING_FREQUENCY) == 0 && totalGenerationsElapsed > 0 ) || ( (totalGenerationsElapsed+1) == (nMUTATIONS * nGENERATIONS_MAX) ) || ( RECORDING_TIMES_IN_FILE && totalGenerationsElapsed == nextRecordingTime ) )
        calcDXY2( &alleleFrequenciesByPatch[0][0] );
    
    if ( ((totalGenerationsElapsed % TS_SAMPLING_FREQUENCY) == 0 && totalGenerationsElapsed > 0 ) || ( (totalGenerationsElapsed+1) == (nMUTATIONS * nGENERATIONS_MAX) ) || ( RECORDING_TIMES_IN_FILE && totalGenerationsElapsed == nextRecordingTime ) ) {
        for ( i = 0; i < nLOCI; i++ ) {
            if ( variable_loci[i] && ((locusID[i] % AL_FREQ_SUBSAMPLE) == 0) ) {
                if ( FST[i] > FST_MIN_RECORDING_THRESH ) {
                    locType = IS_SELECTED_LOCUS[i];
                    
                    // stats for nearest selected neighboring locus
                    if ( GAMETE_PRODUCTION_MODE > 0 || PURE_NEUTRAL_MODEL )
                        nsn = -1;
                    else
                        nsn = findNearestSelectedNeighbor(i);
                    
                    if ( nsn > 0 ) {
                        nsnDist = fabs( (MAP[i] - MAP[nsn]) );
                        nsnS = S_MAX1[nsn];
                        nsnLD = calculateLDpairOneOff(i, nsn);
                    }
                    else {
                        nsnDist = 50.0;
                        nsnS = 0.0;
                        nsnLD = 0.0;
                    }
                    
                    
                    fprintf(FSTtimeSeries,"%li %li %E %E %E %E %i %E %i\n", totalGenerationsElapsed, locusID[i], FST[i], allele_frequencies[i], S_MAX1[i], S_MAX0[i], chromosomeMembership[i], MAP[i], locType);
                    
                    AFdiff = alleleFrequenciesByPatch[i][0] - alleleFrequenciesByPatch[i][1];
                    
                    fprintf(AFtimeSeries, "%li %li", totalGenerationsElapsed, locusID[i]);
                    for ( j = 0; j < nPATCHES; j++ )
                        fprintf(AFtimeSeries, " %E", alleleFrequenciesByPatch[i][j]);
                    fprintf(AFtimeSeries, " %i %i %E %E\n", is_reversed_locus[i], locType, allele_frequencies[i], AFdiff);
                    
                    
                    if ( locType )
                        fprintf(selectedFrequencies,"%li %li %E %E %E %E %E %E %E\n", totalGenerationsElapsed, locusID[i], allele_frequencies[i], AFdiff, S_MAX1[i], FST[i], nsnDist, nsnS, nsnLD);
                    else
                        fprintf(neutralFrequencies,"%li %li %E %E %E %E %E\n", totalGenerationsElapsed, locusID[i], allele_frequencies[i], AFdiff, nsnDist, nsnS, nsnLD);
                }
            }
        }
    }
    
    taln = totalAlleles[newestLocus];
    
    if ( nRemovals > 0 ) {
        for ( i = (nRemovals-1); i >= 0; i-- ) { // go backwards so that indexes make sense even after teh first removal
            removeALocus(lociToRemove[i]);
        }
    }
    
    if ( taln < 1 )
        return 1; // latest allele introduced was lost
    else
        return 0;
}



void calculateFitnesses(double *f, double *fsum)
{
    int i, j, p, l, locus, locus2, nEpistaticInteractions;
    short int *gtpt, gt, gt1, gt2, epiInteraction;
    short int selectedLociGTs[nLOCI];
    double *dpt = f, fitVal, ec, fv1, fv2, resid, ecmult;
    long int fitArrayIndex[nLOCI], SELECTED_LOCI[nSELECTED_LOCI];
    _Bool locusTakenCareOf[nLOCI], invertFitness, useOppositeArray, timeToRecordSstar;
    
    if ( ((totalGenerationsElapsed % TS_SAMPLING_FREQUENCY) == 0 && totalGenerationsElapsed > 0 ) || ( (totalGenerationsElapsed+1) == (nMUTATIONS * nGENERATIONS_MAX) ) || ( RECORDING_TIMES_IN_FILE && totalGenerationsElapsed == nextRecordingTime ) )
        timeToRecordSstar = 1;
    else
        timeToRecordSstar = 0;
    
    int bar = 1;
    p = 0;
    for ( i = 0; i < nLOCI; i++ ) {
        fitArrayIndex[i] = 0;
        if ( IS_SELECTED_LOCUS[i] ) {
            SELECTED_LOCI[p] = i;
            p++;
        }
        selectedLociGTs[i] = 0;
        locusTakenCareOf[i] = 0;
    }
    
    for ( i = 0; i < N; i++ ) { // over all individuals
        fitVal = 1.0; // base fitness
        p = patch_locations[i]; // get patch number
        if ( !PURE_NEUTRAL_MODEL ) {
            for ( l = 0; l < nSELECTED_LOCI; l++ ) {
                locus = SELECTED_LOCI[l];
                gtpt = genotypes + ( i * 2 * nLOCI ) + ( 2 * locus );
                gt = (*gtpt) + (*(gtpt+1)); // genotype represented as a scalar; works as long as only possible alleles are 0 and 1
                selectedLociGTs[locus] = gt; // stored for calculations;
                fitArrayIndex[locus] = ( locus * nPATCHES ) + p;
                locusTakenCareOf[locus] = 0;
            }
            // now we have this individual's genotype at each locus represented as a scalar
#if ( MULTIPLICATIVE_FITNESS )
            if ( CONSIDER_EPISTASIS ) {
                for ( l = 0; l < nSELECTED_LOCI; l++ ) {
                    locus = SELECTED_LOCI[l];
                    gt1 = selectedLociGTs[locus];
                    if ( any_epis_for_this_locus[locus] ) {
                        if ( l < (nSELECTED_LOCI - 1) ) {
                            // we need to account for epistatic interactions
                            for ( j = l+1; j < nSELECTED_LOCI; j++ ) { // ensures we don't count an interaction twice
                                locus2 = SELECTED_LOCI[j];
                                gt2 = selectedLociGTs[locus2];
                                epiInteraction = *(epistasisMatrix + (locus * nLOCI) + locus2); // locus-th row, locus2-th column
                                if ( epiInteraction == -1 ) { // DMI
                                    if ( gt1 && gt2 ) { // this individual has at least one derived allele at each of the loci
                                        ec = *(epi_coeff_matrix + (nLOCI * locus) + locus2);
                                        fitVal *= (1.0 - ec);  // reduction of 100*ec% in fitness of individual
                                        locusTakenCareOf[locus] = 1;
                                        locusTakenCareOf[locus2] = 1;
                                    }
                                }
                                else if ( epiInteraction == 1 ) { // positive interaction
                                    if ( (gt1 > 0) && (gt2 > 0) ) { // this individual has at least one derived allele at each of the loci
                                        ec = *(epi_coeff_matrix + (nLOCI * locus) + locus2);
                                        //fprintf(stderr,"\nt = %li\nec = %f\n", totalGenerationsElapsed, ec);
                                        // E = 1 + ec;
                                        ecmult = epi_patch_multipliers[p];
                                        //fprintf(stderr,"ecmult = %f\n", ecmult);
                                        if ( is_reversed_locus[locus])
                                            ecmult = -ecmult;
                                        
                                        if ( ecmult < 0.0 ) { // derived allele is bad in this patch
                                            invertFitness = 1;
                                            if ( is_reversed_locus[locus] )
                                                useOppositeArray = 0; // double switch: reversed and inverted
                                            else
                                                useOppositeArray = 1; // single switch: inverted but regular
                                        }
                                        else { // derived allele is good in this patch
                                            invertFitness = 0;
                                            if ( is_reversed_locus[locus] )
                                                useOppositeArray = 1; // single switch: reversed but not inverted
                                            else
                                                useOppositeArray = 0; // no switch: not reversed, not inverted
                                        }
                                        
                                        ecmult = fabs(ecmult); // now we just want the magnitude
                                        ec = 1.0 + (ecmult * ec);
                                        
                                        //fprintf(stderr, "is_reversed_locus[%i] = %i, ecmult = %f, ec = %f\n", locus, is_reversed_locus[locus], ecmult, ec);
                                        
                                        // ecmult is used to create symmetry across the habitat: the enhancements in the right deme
                                        // are offset by reductions in the wrong deme
                                        // ec is now a multiplier that can be used below
                                        
                                        // locus 1
                                        if ( gt1 == 2 ) { // homozygous 1
                                            if ( useOppositeArray )
                                                fv1 = (*(fit0 + fitArrayIndex[locus]));
                                            else
                                                fv1 = (*(fit1 + fitArrayIndex[locus]));
                                        }
                                        else if ( gt1 == 1 ) // heterozygous
                                            fv1 = (*(fitH + fitArrayIndex[locus]));
                                        else { // homozygous 0
                                            fprintf(stderr,"\nMessed up, shouldn't have gotten here 1.\n");
                                            fprintf(stderr,"locus = %i, locus2 = %i, gt1 = %i, gt2 = %i\n", locus, locus2, gt1, gt2);
                                            exit(1);
                                        }
                                        
                                        // locus 2
                                        if ( gt2 == 2 ) { // homozygous 1
                                            if ( useOppositeArray )
                                                fv2 = (*(fit0 + fitArrayIndex[locus2]));
                                            else
                                                fv2 = (*(fit1 + fitArrayIndex[locus2]));
                                        }
                                        else if ( gt2 == 1 ) // heterozygous
                                            fv2 = (*(fitH + fitArrayIndex[locus2]));
                                        else { // homozygous 0
                                            fprintf(stderr,"\nMessed up, shouldn't have gotten here. 2\n");
                                            fprintf(stderr,"locus = %i, locus2 = %i, gt1 = %i, gt2 = %i\n", locus, locus2, gt1, gt2);
                                            exit(1);
                                        }
                                        //fprintf(stderr, "p = %i, fv1 = %f, fv2 = %f, fv1*fv2 = %f, fitVal = %f\n", p, fv1, fv2, (fv1 * fv2), fitVal);
                                        
                                        if ( fv1 < 1.0 || fv2 < 1.0 ) {
                                            fprintf(stderr,"\nError!  Bad positive interaction calc!\n\tfv1 = %f, fv2 = %f\n", fv1, fv2);
                                            exit(1);
                                        }
                                        
                                        resid = (fv1 * fv2) - 1.0; // get the increase in fitness from this pair of loci
                                        if ( invertFitness ) // negative effect of epistasis
                                            fitVal *= ((1.0 + resid) / (1.0 + (ec * resid))); // fitness diminished by percentage that epi would have increased it in the reciprocal habitat
                                        else
                                            fitVal *= (1.0 + (ec * resid));
                                        // if ec < 1, then the fitness is less than it would be
                                        //fprintf(stderr, "resid = %f, fitmult = %f, fitVal = %f\n", resid, (1.0 + (ec * resid)), fitVal);
                                        //if ( is_reversed_locus[locus] && p == 1 )
                                        //exit(0);
                                        
                                        locusTakenCareOf[locus] = 1;
                                        locusTakenCareOf[locus2] = 1;
                                    }
                                }
                            }
                        }
                    }
                }
            }
            // do the following, regardless of epistasis
            for ( l = 0; l < nSELECTED_LOCI; l++ ) {
                locus = SELECTED_LOCI[l];
                if ( !locusTakenCareOf[locus] ) { // effects of this locus on fitness have NOT yet been accounted for
                    gt = selectedLociGTs[locus];
                    if ( gt == 2 ) { // homozygous 1
                        if ( is_reversed_locus[locus] )
                            fitVal = fitVal * (*(fit0 + fitArrayIndex[locus]));
                        else
                            fitVal = fitVal * (*(fit1 + fitArrayIndex[locus]));
                    }
                    else if ( gt == 1 ) // heterozygous
                        fitVal = fitVal * (*(fitH + fitArrayIndex[locus]));
                    else { // homozygous 0
                        if ( is_reversed_locus[locus] )
                            fitVal = fitVal * (*(fit1 + fitArrayIndex[locus]));
                        else
                            fitVal = fitVal * (*(fit0 + fitArrayIndex[locus]));
                    }
                }
            }
            if ( fitVal < 0.0 ) {
                fprintf(stderr, "\nNegative fitness: fitval = %E\n", fitVal);
                fitVal = 0.0;
            }
            
            
#elif ( ADDITIVE_FITNESS )
            if ( CONSIDER_EPISTASIS ) {
                fprintf(stderr,"\nCONSIDER_EPISTASIS = %i\n", CONSIDER_EPISTASIS);
                fprintf(stderr, "\nHey dude, you didn't code all the calculations for additive fitness with epistasis yet.\n\tCode it and then compile again.\n\t... Exiting ...\n\n");
                exit(1);
            }
            // code below is old and needs to be redone for epistasis; might have also deleted some stuff from old version, so check all
            //                if ( gt == 2 ) { // homozygous 1
            //                    if ( is_reversed_locus[locus] )
            //                        fitVal = fitVal + (*(fit0 + fitArrayIndex));
            //                    else
            //                        fitVal = fitVal + (*(fit1 + fitArrayIndex));
            //                }
            //                else if ( gt == 1 ) // heterozygous
            //				fitVal = fitVal + (*(fitH + fitArrayIndex));
            //                else { // homozygous 0
            //                    if ( is_reversed_locus[locus] )
            //                        fitVal = fitVal + (*(fit1 + fitArrayIndex));
            //                    else
            //                        fitVal = fitVal + (*(fit0 + fitArrayIndex));
            //                }
            
#endif
        }
        // fitVal is now individual i's fitness over all loci
        *(fsum + p) = (*(fsum + p)) + fitVal;
        *dpt = (*(fsum + p)); // make values in fitness vector CUMULATIVE
        dpt++; // increment pointer
    }
}


void calculateLD(int gatherLDvalues)
{
    int i, j, locus, l1, l2, lociToUse[nVariableLoci];
    long int id1, id2;
    int totalVariable = 0;
    double dist, *Dvals, *DprimeVals, *DeltaVals, *dpt1, *dpt2, *dpt3;
    long int totNum = 0;
    
    //#if ( KEEP_REINTRODUCING == 1)
    //	if ( numberOfMutationTries == 2 )
    //#else
    
    for ( i = 0; i < nLOCI; i++ ) {
        locus = i;
        if ( variable_loci[locus] ) {
            lociToUse[totalVariable] = locus;
            totalVariable++;
        }
    }
    
    
    if ( gatherLDvalues >= 1 ) {
        if ( totalVariable > 1 ) {
            
            Dvals = (double *) malloc( sizeof(double) * (totalVariable * (totalVariable - 1) / 2) );
            DprimeVals = (double *) malloc( sizeof(double) * (totalVariable * (totalVariable - 1) / 2) );
            DeltaVals = (double *) malloc( sizeof(double) * (totalVariable * (totalVariable - 1) / 2) );
            dpt1 = Dvals;
            dpt2 = DprimeVals;
            dpt3 = DeltaVals;
            
            for ( i = 0; i < (totalVariable - 1); i++ ) {
                l1 = lociToUse[i];
                id1 = locusID[l1];
                for ( j = (i+1); j < totalVariable; j++ ) {
                    l2 = lociToUse[j];
                    id2 = locusID[l2];
                    
                    if ( (id1 % LD_LOCI_SUBSAMPLE == 0) && (id2 % LD_LOCI_SUBSAMPLE == 0) ) {
                        if ( chromosomeMembership[l1] == chromosomeMembership[l2] )
                            dist = fabs( MAP[l2] - MAP[l1] );
                        else
                            dist = -1.0;
                        
                        calculateLDpair(l1,l2,dist, dpt1, dpt2, dpt3, gatherLDvalues);
                        totNum++;
                        dpt1++;
                        dpt2++;
                        dpt3++;
                    }
                }
            }
            
            fprintf(effMigRates, " %E %E %E %E %E %E\n", calcMeanMagnitude(Dvals, totNum), calcVariance(Dvals, totNum), calcMeanMagnitude(DprimeVals, totNum), calcVariance(DprimeVals, totNum), calcMeanMagnitude(DeltaVals, totNum), calcVariance(DeltaVals, totNum) );
            
            free(Dvals);
            free(DprimeVals);
            free(DeltaVals);
        }
        else
            fprintf(effMigRates, " 0.0 0.0 0.0 0.0 0.0 0.0\n");
    }
    
}


void calculateLDneutralSitesOnly(int gatherLDvalues)
{
    int i, j, locus, l1, l2, lociToUse[nNEUTRAL_LOCI];
    long int id1, id2;
    int totalVariable = 0;
    double DeltaValsSame = 0.0, val, dist;
    double DeltaValsDiff = 0.0;
    double totNumSame = 0.0, totNumDiff = 0.0; // doubles since only purpose is in division operation.
    int localLDSS = LD_LOCI_SUBSAMPLE;
    
    //    if ( N <= 1000 ) { //drifty
    //        localLDSS = LD_LOCI_SUBSAMPLE * 10;
    //    }
    
    
    //#if ( KEEP_REINTRODUCING == 1)
    //	if ( numberOfMutationTries == 2 )
    //#else
    
    for ( locus = 0; locus < nLOCI; locus++ ) {
        if ( !IS_SELECTED_LOCUS[locus] && !is_reference_locus[locus] ) {
            lociToUse[totalVariable] = locus;
            totalVariable++;
        }
        if ( totalVariable > nNEUTRAL_LOCI ) {
            fprintf(stderr, "\nError in calculateLDneutralSitesOnly(): \n\t\t totalVariable %i > nNEUTRAL_LOCI = %i\n", totalVariable, nNEUTRAL_LOCI);
            exit(-1);
        }
    }
    if ( totalVariable != nNEUTRAL_LOCI ) {
        fprintf(stderr, "\nError in calculateLDneutralSitesOnly(): \n\t\t totalVariable %i != nNEUTRAL_LOCI = %i\n", totalVariable, nNEUTRAL_LOCI);
        exit(-1);
    }
    
    
    
    if ( totalVariable > 1 ) {
        
        for ( i = 0; i < (totalVariable - 1); i++ ) {
            l1 = lociToUse[i];
            id1 = locusID[l1];
            for ( j = (i+1); j < totalVariable; j++ ) {
                l2 = lociToUse[j];
                id2 = locusID[l2];
                if ( (id1 % localLDSS == 0) && (id2 % localLDSS == 0) ) {
                    val = fabs( calculateLDpairOneOff(l1, l2) );
                    if ( chromosomeMembership[l1] == chromosomeMembership[l2] ) {
                        DeltaValsSame += val;
                        totNumSame += 1.0;
                        dist = fabs( MAP[l2] - MAP[l1] );
                        if ( gatherLDvalues >= 3 )
                            fprintf(LDneutSitesSame, "%li %E %E\n", totalGenerationsElapsed, val, dist);
                    }
                    else {
                        DeltaValsDiff += val;
                        totNumDiff += 1.0;
                        if ( gatherLDvalues >= 3 )
                            fprintf(LDneutSitesDiff, "%li %E\n", totalGenerationsElapsed, val);
                    }
                }
            }
        }
        
        fprintf( LDneutSitesAvg, "%li ", totalGenerationsElapsed );
        if ( totNumSame > 0.0 )
            fprintf( LDneutSitesAvg, "%E ", (DeltaValsSame / totNumSame) );
        else
            fprintf( LDneutSitesAvg, "0.0 ");
        
        if ( totNumDiff > 0.0 )
            fprintf( LDneutSitesAvg, "%E ", (DeltaValsDiff / totNumDiff) );
        else
            fprintf( LDneutSitesAvg, "0.0 ");
        
        if ( (totNumDiff + totNumSame) > 0.0 )
            fprintf( LDneutSitesAvg, "%E\n", ((DeltaValsDiff + DeltaValsSame) / (totNumDiff + totNumSame)) );
        else
            fprintf( LDneutSitesAvg, "0.0\n");
        
    }
}



void calculateLDselectedSitesOnly(int gatherLDvalues)
{
    int i, j, locus, l1, l2, lociToUse[nSELECTED_LOCI];
    long int id1, id2;
    int totalVariable = 0;
    double DeltaValsSame = 0.0, val, dist;
    double DeltaValsDiff = 0.0;
    double totNumSame = 0.0, totNumDiff = 0.0; // doubles since only purpose is in division operation.
    
    //#if ( KEEP_REINTRODUCING == 1)
    //	if ( numberOfMutationTries == 2 )
    //#else
    
    for ( locus = 0; locus < nLOCI; locus++ ) {
        if ( IS_SELECTED_LOCUS[locus] ) {
            lociToUse[totalVariable] = locus;
            totalVariable++;
        }
    }
    if ( totalVariable != nSELECTED_LOCI ) {
        fprintf(stderr, "\nError in calculateLDselectedSitesOnly(): totalVariable %i != nSELECTED_LOCI %i\n", totalVariable, nSELECTED_LOCI);
        exit(-1);
    }
    
    
    
    if ( totalVariable > 1 ) {
        
        for ( i = 0; i < (totalVariable - 1); i++ ) {
            l1 = lociToUse[i];
            id1 = locusID[l1];
            for ( j = (i+1); j < totalVariable; j++ ) {
                l2 = lociToUse[j];
                id2 = locusID[l2];
                if ( (id1 % LD_LOCI_SUBSAMPLE == 0) && (id2 % LD_LOCI_SUBSAMPLE == 0) ) {
                    val = fabs( calculateLDpairOneOff(l1, l2) );
                    if ( chromosomeMembership[l1] == chromosomeMembership[l2] ) {
                        DeltaValsSame += val;
                        totNumSame += 1.0;
                        dist = fabs( MAP[l2] - MAP[l1] );
                        if ( gatherLDvalues >= 3 )
                            fprintf(LDselSitesSame, "%li %E %E\n", totalGenerationsElapsed, val, dist);
                    }
                    else {
                        DeltaValsDiff += val;
                        totNumDiff += 1.0;
                        if ( gatherLDvalues >= 3 )
                            fprintf(LDselSitesDiff, "%li %E\n", totalGenerationsElapsed, val);
                    }
                }
            }
        }
        
        fprintf( LDselSitesAvg, "%li ", totalGenerationsElapsed );
        if ( totNumSame > 0.0 )
            fprintf( LDselSitesAvg, "%E ", (DeltaValsSame / totNumSame) );
        else
            fprintf( LDselSitesAvg, "0.0 ");
        
        if ( totNumDiff > 0.0 )
            fprintf( LDselSitesAvg, "%E ", (DeltaValsDiff / totNumDiff) );
        else
            fprintf( LDselSitesAvg, "0.0 ");
        
        if ( (totNumDiff + totNumSame) > 0.0 )
            fprintf( LDselSitesAvg, "%E\n", ((DeltaValsDiff + DeltaValsSame) / (totNumDiff + totNumSame)) );
        else
            fprintf( LDselSitesAvg, "0.0\n");
        
    }
}


void calculateLDpair(int l1, int l2, double dist, double *Dvals, double *DprimeVals, double *DeltaVals, int gatherLDvalues)
{
    int i, j;
    double zz = 0.0, oo = 0.0, zo = 0.0, oz = 0.0;
    double Ninv = 0.5 / ((double) N);
    double exp11, DD, Dmax, Dprime, Delta;
    short int *gpt1, *gpt2;
    
    for ( j = 0; j < 2; j++ ) {
        
        gpt1 = genotypes + (2 * l1) + j; // allele of first locus on chromosome
        gpt2 = genotypes + (2 * l2) + j; // allele on second locus on chromosome
        
        for ( i = 0; i < N; i++ ) {
            if ( *gpt1 ) {
                if ( *gpt2 ) {
                    oo++;
                }
                else {
                    oz++;
                }
            }
            else {
                if ( *gpt2 ) {
                    zo++;
                }
                else {
                    zz++;
                }
            }
            gpt1 += (2 * nLOCI);
            gpt2 += (2 * nLOCI);
        }
        
    }
    
    // haplotype frequencies
    oo *= Ninv; // observed 11
    oz *= Ninv; // observed 10
    zo *= Ninv; // observed 01
    zz *= Ninv; // observed 00
    
    double A = allele_frequencies[l1], B = allele_frequencies[l2];
    
    DD = (( oo * zz ) - ( zo * oz ));
    exp11 = A * B; // expected frequency of 11 haplotype
    
    // check
    //if ( DD != (oo - exp11) )
    //	printf("bad math, DD = %G, diff = %G\n", DD, (oo - exp11));
    
    if ( DD > 0.0 )
        Dmax = fmin( (A * (1.0 - B)), ((1.0 - A) * B) );
    else
        Dmax = fmin( (A * B), ((1.0 - A) * (1.0 - B)) );
    
    Dprime = DD / Dmax;
    
    Delta = DD / sqrt( A * B * (1.0 - A) * (1.0 - B) );
    
    if ( gatherLDvalues >= 3 ) {
        if ( BeginRecordingLD ) {
            if ( fabs(DD) > LD_LowerBound ) {
                fprintf(LDfpt, "%li %li %li %E %E %E %E %E %E %i %i\n", totalGenerationsElapsed, locusID[l1], locusID[l2], DD, Dprime, Delta, dist, S_MAX1[l1], S_MAX1[l2], IS_SELECTED_LOCUS[l1], IS_SELECTED_LOCUS[l2]);
            }
        }
    }
    *DeltaVals = Delta;
    *Dvals = DD;
    *DprimeVals = Dprime;
}


double calculateLDpairOneOff(int l1, int l2)
{
    int i, j;
    double zz = 0.0, oo = 0.0, zo = 0.0, oz = 0.0;
    double Ninv = 0.5 / ((double) N);
    double exp11, DD, Dmax, Dprime, Delta;
    short int *gpt1, *gpt2;
    double A = allele_frequencies[l1], B = allele_frequencies[l2];
    
    for ( j = 0; j < 2; j++ ) {
        
        gpt1 = genotypes + (2 * l1) + j; // allele of first locus on chromosome
        gpt2 = genotypes + (2 * l2) + j; // allele on second locus on chromosome
        
        for ( i = 0; i < N; i++ ) {
            if ( *gpt1 ) {
                if ( *gpt2 ) {
                    oo++;
                }
                else {
                    oz++;
                }
            }
            else {
                if ( *gpt2 ) {
                    zo++;
                }
                else {
                    zz++;
                }
            }
            gpt1 += (2 * nLOCI);
            gpt2 += (2 * nLOCI);
        }
        
    }
    
    // haplotype frequencies
    oo *= Ninv; // observed 11
    oz *= Ninv; // observed 10
    zo *= Ninv; // observed 01
    zz *= Ninv; // observed 00
    
    DD = (( oo * zz ) - ( zo * oz ));
    Delta = DD / sqrt( A * B * (1.0 - A) * (1.0 - B) );
    
    
    
    //	exp11 = A * B; // expected frequency of 11 haplotype
    //
    //	// check
    //	//if ( DD != (oo - exp11) )
    //	//	printf("bad math, DD = %G, diff = %G\n", DD, (oo - exp11));
    //
    //	if ( DD > 0.0 )
    //		Dmax = fmin( (A * (1.0 - B)), ((1.0 - A) * B) );
    //	else
    //		Dmax = fmin( (A * B), ((1.0 - A) * (1.0 - B)) );
    //
    //	Dprime = DD / Dmax;
    
    return Delta;
}



void makeZygoteChromosomes(int parent, short int *ogtpt)
{
    int i, j, l, csome, ncos, nl, firstl, lastl, co_count, count;
    short int *ppt, *opt, *startpt;
    double ml, co, spot; // lastco;
    
    opt = ogtpt; // pointer to first allele that offspring will get for first locus of chromosome
    
    //
    // First check the simple/fast case where we're sure
    // there are no variable loci.
    //
    if (nVariableLoci == 0) {
        //
        // Following code (through while loop) is only needed
        // to ensure compatibility in deterministic mode with
        // older code.  Keeps RNG in sync, ensuring same results.
        //
        volatile double jnk;
        
        if ( GAMETE_PRODUCTION_MODE == 1 )
            i = 1;
        else
            i = nCHROMOSOMES * 2;
        
        while (i--)
            jnk  = randU();
        // end of compatibility code
        
        ppt = genotypes + ( 2 * nLOCI * parent );
        for ( i = 0; i < nLOCI; i++ ) {
            *opt = (*ppt); // put the allele in the offspring's genome
            ppt += 2;
            opt += 2;
        }
        return;
    }
    
    if ( GAMETE_PRODUCTION_MODE == 1 ) { // Genome hitchhiking only; no linkange less than 50 cM;
        ppt = genotypes + ( 2 * nLOCI * parent );
        if ( randU() < 0.5 ) {
            csome = 0;
            *opt = (*ppt);
        }
        else {
            csome = 1;
            ppt++;
            *opt = (*ppt);
        }
        opt += 2; // first locus done
        
#ifdef notdef
        int loci_count = 0;
        for ( i = 1; i < nLOCI; i++ ) {
            if ( variable_loci[i] )
                loci_count++;
        }
        if (loci_count == 0)
            loci_count_0++;
        if (loci_count == 1)
            loci_count_1++;
        if (loci_count == 2)
            loci_count_2++;
        if (loci_count > 2)
            loci_count_many++;
        
        if (loci_count != nVariableLoci)
            printf("variable loci %d != %d\n", loci_count, nVariableLoci);
        
#endif
        
        //
        // First check the simple/fast case where we're sure
        // there are no variable loci.
        //
        if (nVariableLoci == 0) {
            for ( i = 1; i < nLOCI; i++ ) {
                ppt += 2;
                *opt = (*ppt); // put the allele in the offspring's genome
                opt += 2;
            }
            return;
        }
        
        //
        // Failed the easy case, so need to check each locus
        //
        for ( i = 1; i < nLOCI; i++ ) {
            // first figure out which allele to sample from the parent
            if ( variable_loci[i] ) { // it might matter which allele is chosen
                if ( randU() < 0.5 ) { // switch to other set of alleles, like a crossover
                    if ( csome ) { // csome was 1; need to only advance 1
                        ppt++;
                        csome = 0;
                    }
                    else { // csome was 0; need to advance 3
                        ppt += 3;
                        csome = 1;
                    }
                }
                else
                    ppt += 2; // no switch, advance 2
            }
            else // both alleles same
                ppt += 2;
            
            *opt = (*ppt); // put the allele in the offspring's genome
            
            opt += 2;
        }
    }
    
    else { // linkage matters --> DH and GH (-G 0 on command line)
        // get starting chromosome
        count = 0;
        firstl = 1; // array element to start with (0 taken care of manually before relvant for loop)
        lastl = LOCI_PER_CHROMOSOME[0];
        startpt = genotypes + ( 2 * nLOCI * parent ); // pointer to first allele, first chrom, first locus of parent
        for ( i = 0; i < nCHROMOSOMES; i++ ) {
            ppt = startpt;
            nl = LOCI_PER_CHROMOSOME[i];
            // random draw for which chromosome to start with; mendel's law of independent assortment
            if ( randU() < 0.5 ) {
                csome = 0;
            }
            else {
                csome = 1;
                ppt++;
            }
            *opt = (*ppt); // first allele at first locus of this chromosome
            //			count++;
            //			fprintf(stderr, "Count = %i, off allele = %i, parent allele = %i\n",count, (*opt),(*ppt));
            opt += 2; // pointer to next locus on same offspring chromosome (which is being created and stored at the same time in the next block)
            
            //
            // First check the simple/fast case where we're sure
            // there are no variable loci.
            //
            if (nVariableLoci == 0) {
                co = randU(); // deterministic
                for (l = firstl; l < lastl; l++ ) {
                    ppt += 2;
                    *opt = (*ppt); // put the allele in the offspring's genome
                    opt += 2;
                }
            }
            else {
                co = (log(1.0 - randU())) * (-50.0); // crossovers are i.i.d. with a mean of 50 cM between them
                
                for ( l = firstl; l < lastl; l++ ) { // one locus at a time for this chromosome
                    if ( variable_loci[l] ) { // crossovers may actually matter for choosing alleles
                        spot = MAP[l];
                        if ( co < spot ) { // location of crossover is before this locus
                            co_count = 0;
                            while ( co < spot ) { // allow for MULTIPLE crossovers betweeen consecutive loci; hence the while loop
                                //lastco = co; // just for testing
                                co += (log(1.0 - randU())) * (-50.0); // spot for NEXT crossover
                                co_count++; // count number of crossovers between last locus l-1 and next locus l
                            }
                            if ( (co_count % 2) == 1 ) { // if number is odd
                                if ( csome ) {
                                    csome = 0;
                                    ppt++; // only have to advance 1 spot forward to get to 0 allele on this chromosome from 1 allele on last one
                                }
                                else {
                                    csome = 1;
                                    ppt += 3; // have to advance pointer 3 to go from 0 allele on last chromosome to 1 allele on this chromosome
                                }
                            }
                            else // even number of C-Os is same as none at all
                                ppt += 2;
                            //fprintf(stderr, "\n%i crossovers on chrom. %i. last at %G before locus %i at map loc. %G", co_count, i, lastco, l, MAP[l]);
                        }
                        else // no crossover
                            ppt += 2;
                    }
                    else // there is only one type of allele; save time by copying one
                        ppt += 2;
                    
                    *opt = (*ppt); // offspring gets selected allele
                    //				count++;
                    //				fprintf(stderr, "Count = %i, off allele = %i, parent allele = %i\n",count, (*opt),(*ppt));
                    opt += 2; // increment offspring pointer.
                }
                if ( i < (nCHROMOSOMES-1) ) { // if we are going through the loop again
                    startpt = startpt + ( 2 * nl ); // pointer to first allele on next chromosome of the PARENT
                    firstl = lastl + 1;
                    lastl = lastl + LOCI_PER_CHROMOSOME[(i+1)];
                }
            }
        }
    }
    //	exit(0);
}



void move(void)
{
    int i, j, *ipt1, *ipt2;
    double direction, distance, x, y, *xpt, *ypt;
    
    totalMigrants = 0;
    
    for ( i = 0; i < nPATCHES; i++ )
        migrationCount[i] = 0;
    
    if ( PURE_ALLOPATRY_MODEL ) {
        ipt1 = previousPatches;
        ipt2 = patch_locations;
        for ( i = 0; i < N; i++ ) {
            *ipt1 = (*ipt2);
            ipt1++;
            ipt2++;
        }
    }
    else {
        
        if ( TWO_DEME ) { // replace with preprocessor #if statment?
            ipt1 = previousPatches;
            ipt2 = patch_locations;
            xpt = x_locations;
            
            for ( i = 0; i < N; i++ ) {
                *ipt1 = (*ipt2);
                if ( totalGenerationsElapsed > END_PERIOD_ALLOPATRY ) {
                    if ( randU() < SD_MOVE ) { // in this case, SD_MOVE is treated as the probability of migration per individual per generation
                        totalMigrants++; // count the movement
                        if ( (*ipt2) ) { // patch was #1
                            *ipt2 = 0; // actuate migration event
                            *xpt = 0.25;
                            migrationCount[0]++;
                            n_in_each_patch[0]++;
                            n_in_each_patch[1]--;
                        }
                        else { // patch was #0
                            *ipt2 = 1;
                            *xpt = 0.75;
                            migrationCount[1]++;
                            n_in_each_patch[1]++;
                            n_in_each_patch[0]--;
                        }
                    }
                }
                xpt++;
                ipt1++;
                ipt2++;
            }
        }
        
        else {
            ipt1 = previousPatches;
            ipt2 = patch_locations;
            xpt = x_locations;
#if ( D == 2 )
            ypt = y_locations;
#endif
            
            for ( i = 0; i < N; i++ ) {
                *ipt1 = (*ipt2);
                
#if ( D == 1 )
                if ( totalGenerationsElapsed > END_PERIOD_ALLOPATRY ) {
                    x = boxMuller((*xpt),SD_MOVE);
                    if ( x < 0.0 ) {
                        *xpt = 0.0;
                        *ipt2 = 0;
                    }
                    else if ( x > 1.0 ) {
                        *xpt = 1.0;
                        *ipt2 = (nPATCHES-1);
                    }
                    else {
                        *xpt = x;
                        *ipt2 = (int) ( x * ((double) PATCHES) );
                        if ( (*ipt2) == PATCHES ) // would happen if *xpt == 1.0
                            *ipt2 = PATCHES - 1;
                    }
                }
#elif ( D == 2 )
                if ( totalGenerationsElapsed > END_PERIOD_ALLOPATRY ) {
                    direction = TWOPI * randU(); // direction in radians
                    distance = fabs( boxMuller(0.0, SD_MOVE) );
                    
                    x = (*xpt);
                    y = (*ypt);
                    
                    x = x + (cos(direction) * distance);
                    if ( x > 1.0 )
                        *xpt = 1.0;
                    else if ( x < 0.0 )
                        *xpt = 0.0;
                    else
                        *xpt = x;
                    
                    y = y + (sin(direction) * distance);
                    if ( y > 1.0 )
                        *ypt = 1.0;
                    else if ( y < 0.0 )
                        *ypt = 0.0;
                    else
                        *ypt = y;
                    
                    
                    *ipt2 = scalarPatchNumber((*xpt),(*ypt));
                    
                    ypt++;
                }
#endif
                
                if ( totalGenerationsElapsed > END_PERIOD_ALLOPATRY ) {
                    if ( (*ipt2) != (*ipt1) ) {
                        totalMigrants++; // count the movement
                        migrationCount[(*ipt2)]++; // counts migrants INTO each patch
                        n_in_each_patch[(*ipt2)]++; // increment n in new patch
                        n_in_each_patch[(*ipt1)]--; // decrement n in old patch
                    }
                }
                xpt++;
                ipt1++;
                ipt2++;
                
            }
        }
    }
    
    
    if ( (totalGenerationsElapsed % TS_SAMPLING_FREQUENCY == 0 && totalGenerationsElapsed > 0 ) || ( (totalGenerationsElapsed+1) == (nMUTATIONS * nGENERATIONS_MAX) ) || ( RECORDING_TIMES_IN_FILE && totalGenerationsElapsed == nextRecordingTime ) ) {
        double migrationRates[nPATCHES];
        for ( i = 0; i < nPATCHES; i++ ) {
            if ( n_in_each_patch[i] < 1 )
                fprintf(stderr,"\nERROR!  %i individuals in patch %i\n", n_in_each_patch[i],i);
            
            if ( migrationCount[i] > 0 ) {
                migrationRates[i] = ((double) migrationCount[i]) / ((double) n_in_each_patch[i]);
            }
            else {
                migrationRates[i] = 0.0;
            }
        }
        
        //fprintf(grossMigRates,"%li", totalGenerationsElapsed);
        for ( i = 0; i < nPATCHES; i++ )
            fprintf(grossMigRates,"%G ",migrationRates[i]);
        fprintf(grossMigRates,"\n");
    }
    
}



void nextMutation(void)
{
    int i, j, locus, pm, dumi, totAlleleCount = 0;
    int luckyOne, firstIndex, alleleToMutate, genotypei, locType;
    double locationOfNextMutation;
    int number_established = nVariableLoci;
    short int *gtpt;
    
    
    locationOfNextMutation = MUTATION_ORDER[m];
    if ( START_WITH_BIG_IN_PLACE && m < NUMBER_BIG_TO_ESTABLISH ) {
        // overrides if needed; may be redundant, but that is OK
        locType = 1;
        MUTATION_TYPE_SEQUENCE[m] = 1;
    }
    else
        locType = MUTATION_TYPE_SEQUENCE[m];
    
    addALocus(locationOfNextMutation, locType);
    
    locus = newestLocus;
    pm = MUTATION_LOCATIONS[m];
    
    
    variable_loci[locus] = 1;
    nVariableLoci++;
    
    
    if ( START_WITH_BIG_IN_PLACE && m < NUMBER_BIG_TO_ESTABLISH ) {
        S_MAX0[locus] = BIG_START_S_VAL;
        S_MAX1[locus] = BIG_START_S_VAL;
        for ( j = 0; j < N; j++ ) {
            gtpt = genotypes + ( 2 * nLOCI * j) + ( 2 * locus );
            if ( ((patch_locations[j] % 2) && !is_reversed_locus[locus]) || (!(patch_locations[j] % 2) && is_reversed_locus[locus])) { // is spot where 1 allele is adaptive, but this only works for 2 patches!!!!
                if ( randU() > SD_MOVE ) {
                    *gtpt = 1;
                    totAlleleCount++;
                }
                if ( randU() > SD_MOVE ) {
                    *(gtpt + 1) = 1;
                    totAlleleCount++;
                }
            }
            else {
                if ( randU() < SD_MOVE ) {
                    *gtpt = 1;
                    totAlleleCount++;
                }
                if ( randU() < SD_MOVE ) {
                    *(gtpt + 1) = 1;
                    totAlleleCount++;
                }
            }
        }
        allele_frequencies[locus] = ((double) totAlleleCount) / ( 2.0 * ((double) N));
        if ( allele_frequencies[locus] > 0.8 || allele_frequencies[locus] < 0.2 ) {  // overall should be about 0.5
            fprintf(stderr,"\nError!  Allele frequency calculated badly = %E\n\n\t\tExiting\n\n", allele_frequencies[locus]);
            exit(1);
        }
    }
    
    else { // it is NOT a "big one" to start near migration-selection balance
        // introduce it in a random patch
        firstIndex = 0;
        if ( pm > 0 ) {
            for ( i = 0; i < pm; i++ ) // assumes binned sort of population by patches
                firstIndex += n_in_each_patch[i]; // first index of first individual of each patch
            // /* // test check
            if ( firstIndex > ( N - n_in_each_patch[(nPATCHES-1)] )  ) {
                fprintf(stderr, "\n\nERROR!  Bad index for mutation patch location!!\n\n");
                exit(1);
            }
            // */ // end test check
        }
        
        // test check
        if ( n_in_each_patch[pm] < 1 )
            fprintf(stderr,"\nERROR!  %i individuals in patch %i\n", n_in_each_patch[pm], pm);
        
        luckyOne = firstIndex + ( randI() % n_in_each_patch[pm] ); // randomly chosen individual
        
        if ( randU() < 0.5 )
            alleleToMutate = 2 * locus;
        else
            alleleToMutate = ( 2 * locus ) + 1;
        
        genotypei = (luckyOne * 2 * nLOCI) + alleleToMutate; // the index of the entry to change in the genotypes array
        
        if ( genotypes[genotypei] ) { // if the 1 allele fixed
            genotypes[genotypei] = 0;
            allele_frequencies[locus] = 1.0 - ( 0.5 / ((double) N) );
        }
        else {
            genotypes[genotypei] = 1; // the mutation is accomplished
            allele_frequencies[locus] = 0.5 / ((double) N);
        }
    }
    
    /* // test
     short int *newGtpt;
     int j;
     fprintf(stderr,"\ngenotypes:\n");
     newGtpt = genotypes;
     for ( i = 0; i < N; i++ ) {
     for ( j = 0; j < (nLOCI*2); j++ ) {
     fprintf(stderr,"%i ", (*newGtpt));
     newGtpt++;
     }
     fprintf(stderr,"\n");
     }
     
     fprintf(stderr, "\n\nMutation at locus %i\n", locus);
     fprintf(stderr, "Lucky one = index %i in patch %i:\n", luckyOne, patch_locations[luckyOne]);
     for ( i = 0; i < (nLOCI * 2); i++ ) {
     fprintf(stderr, "%i  ", genotypes[((luckyOne * 2 * nLOCI)+i)]);
     if ( (i-19) % 20 == 0 )
     fprintf(stderr, "\n");
     }
     fprintf(stderr, "\n\n");
     //exit(-1);
     // */ // end test
}


void removeALocus(int locusToRemove)
{
    int i, j, focalChrom = chromosomeMembership[locusToRemove];
    int firstIndexToMove = locusToRemove + 1;
    int locType = IS_SELECTED_LOCUS[locusToRemove];
    int nMoving = (nLOCI - 1) - locusToRemove;
    short int *srcPt, *destPt;
    double *ept;
    size_t sizeNeeded;
    long int dum = locusID[locusToRemove];
    
    if ( locType )
        nSELECTED_LOCI--;
    else
        nNEUTRAL_LOCI--;
    
    if ( is_reversed_locus[locusToRemove] )
        nReversedLoci--;
    else
        nRegularLoci--;
    
    fprintf(logOfRemovedLoci, "%li %li %E %i %E %E %E %i %i\n", totalGenerationsElapsed, dum, MAP[locusToRemove], focalChrom, S_MAX0[locusToRemove], S_MAX1[locusToRemove], H[locusToRemove], is_reversed_locus[locusToRemove], locType);
    
    if ( nMoving > 0 ) {
        memmove( (allele_frequencies + locusToRemove), (allele_frequencies + firstIndexToMove), (sizeof(double) * nMoving) );
        memmove( (MAP + locusToRemove), (MAP + firstIndexToMove), (sizeof(double) * nMoving) );
        memmove( (S_MAX0 + locusToRemove), (S_MAX0 + firstIndexToMove), (sizeof(double) * nMoving) );
        memmove( (S_MAX1 + locusToRemove), (S_MAX1 + firstIndexToMove), (sizeof(double) * nMoving) );
        memmove( (H + locusToRemove), (H + firstIndexToMove), (sizeof(double) * nMoving) );
        memmove( (IS_SELECTED_LOCUS + locusToRemove), (IS_SELECTED_LOCUS + firstIndexToMove), (sizeof(int) * nMoving) );
        memmove( (chromosomeMembership + locusToRemove), (chromosomeMembership + firstIndexToMove), (sizeof(int) * nMoving) );
        memmove( (variable_loci + locusToRemove), (variable_loci + firstIndexToMove), (sizeof(int) * nMoving) );
        memmove( (is_reference_locus + locusToRemove), (is_reference_locus + firstIndexToMove), (sizeof(int) * nMoving) );
        memmove( (fixed_allele + locusToRemove), (fixed_allele + firstIndexToMove), (sizeof(short int) * nMoving) );
        memmove( (locusID + locusToRemove), (locusID + firstIndexToMove), (sizeof(long int) * nMoving) );
        
        
        sizeNeeded = nPATCHES * nMoving * sizeof(double);
        
        memmove( (fit0 + (nPATCHES * locusToRemove)), (fit0 + ( nPATCHES * (locusToRemove + 1) )), sizeNeeded );
        memmove( (fit1 + (nPATCHES * locusToRemove)), (fit1 + ( nPATCHES * (locusToRemove + 1) )), sizeNeeded );
        memmove( (fitH + (nPATCHES * locusToRemove)), (fitH + ( nPATCHES * (locusToRemove + 1) )), sizeNeeded );
    }
    
    LOCI_PER_CHROMOSOME[focalChrom] = LOCI_PER_CHROMOSOME[focalChrom] - 1;
    
    // adjust genotypes!
    srcPt = genotypes + (locusToRemove * 2) + 2;
    destPt = genotypes + (locusToRemove * 2);
    
    sizeNeeded = 2 * (nLOCI - 1) * sizeof(short int);
    
    for ( i = 0; i < N; i++ ) {
        
        if ( i < (N-1) )
            memmove( destPt, srcPt, sizeNeeded );
        else
            memmove( destPt, srcPt, (2 * nMoving * sizeof(short int)) );
        
        destPt += (2 * ( nLOCI - 1 ));
        srcPt += (2 * nLOCI);
    }
    
    
    if ( CONSIDER_EPISTASIS )
        shrinkEpistasisMatrix(locusToRemove);
    
    if ( nMoving > 0 )
        memmove( (is_reversed_locus + locusToRemove), (is_reversed_locus + firstIndexToMove), (sizeof(_Bool) * nMoving) );
    
    if ( CONSIDER_EPISTASIS ) {
        srcPt = epistasisMatrix;
        for ( i = 0; i < (nLOCI-1); i++ )  // use nLOCI-1 since epistasisMatrix is already "shrunk" but nLOCI not yet decremented
            any_epis_for_this_locus[i] = 0; // going to do a complete redo of this array since there could be losses anywhere
        for ( i = 0; i < (nLOCI-2); i++ ) {
            srcPt = epistasisMatrix + (i * (nLOCI-1)) + (i + 1);
            for ( j = i+1; j < (nLOCI-1); j++ ) {
                if ( (*srcPt) ) {
                    any_epis_for_this_locus[i] = 1;
                    any_epis_for_this_locus[j] = 1;
                }
                srcPt++;
            }
        }
    }
    
    
    nLOCI--;
    
    
    /*
     // test prints
     // the map
     fprintf(stderr,"\nThe NEW map:\n");
     for ( i = 0; i < nLOCI; i++ )
     fprintf(stderr,"%.2f  ", MAP[i]);
     fprintf(stderr, "\n\nThe NEW number of LOCI_PER_CHROMOSOME:\n");
     for (i=0; i<nCHROMOSOMES; i++) {
     fprintf(stderr, "%i  ", LOCI_PER_CHROMOSOME[i]);
     }
     
     fprintf(stderr, "\n\nThe NEW is_reference_locus:\n");
     for (i=0; i<nLOCI; i++) {
     fprintf(stderr, "%i ", is_reference_locus[i]);
     }
     
     fprintf(stderr,"\nLocus IDs:\n");
     for ( i = 0; i < nLOCI; i++ )
     fprintf(stderr,"%li ", locusID[i]);
     fprintf(stderr,"\n");
     
     fprintf(stderr,"\nChromosome Membership\n");
     for ( i = 0; i < nLOCI; i++ )
     fprintf(stderr,"%i ", chromosomeMembership[i]);
     fprintf(stderr,"\n");
     
     fprintf(stderr,"\nVariable Loci:\n");
     for ( i = 0; i < nLOCI; i++ )
     fprintf(stderr,"%i ", variable_loci[i]);
     fprintf(stderr,"\n");
     for ( i = 0; i < nLOCI; i++ )
     if ( variable_loci[i] )
     fprintf(stderr,"%i ", i);
     fprintf(stderr,"\n");
     fprintf(stderr,"nVariableLoci = %i\n",nVariableLoci);
     
     fprintf(stderr,"\nlocus\tfit0\t\t\tfit1\t\t\tfitH\n");
     double *newFit0 = fit0;
     double *newFit1 = fit1;
     double *newFitH = fitH;
     for ( i = 0; i < nLOCI; i++ ) {
     fprintf(stderr,"%i\t",i);
     for ( j = 0; j < nPATCHES; j++ ) {
     fprintf(stderr,"%.3f ", (*newFit0));
     newFit0++;
     }
     fprintf(stderr,"\t\t");
     for ( j = 0; j < nPATCHES; j++ ) {
     fprintf(stderr,"%.3f ", (*newFit1));
     newFit1++;
     }
     fprintf(stderr,"\t\t");
     for ( j = 0; j < nPATCHES; j++ ) {
     fprintf(stderr,"%.3f ", (*newFitH));
     newFitH++;
     }
     fprintf(stderr,"\n");
     }
     
     
     fprintf(stderr,"\n");
     
     //exit(-1);
     // end of test prints
     //  */
}



void reproduce(int gatherLDvalues)
{
    int i, j, n_off_born[nPATCHES], newN = 0, patch, no, nj, ii;
    int firsti, lasti, mom, dad, offCount, count;
    double dumd, fitVal, fs;
    // double nEffectiveMigrants[nPATCHES],
    double dN, fitnesses[N], fitSum[nPATCHES], *dpt1, *dpt2, *dpt3, *dpt4;
    short int *off_genotypes, *spt1, *spt2;
    double off_x[(3*N)]; // more than needed for sure
#if ( D == 2 )
    double off_y[(3*N)];
#endif
    int	*ipt1, *ipt2;
    double offPerParent[N];
    long int nReproducers;
    
    //	for ( i = 0; i < nPATCHES; i++ )
    //		nEffectiveMigrants[i] = 0.0;
    
    
    // first, determine how many new offspring will be born
    if ( FIXED_N )
        newN = N;
    
    for ( i = 0; i < nPATCHES; i++ ) {
        if ( FIXED_N )
            n_off_born[i] = n_in_each_patch[i];
        else { // density dependent (logistic) growth
            // determine growth rate and number of offspring
            // logistic growth in discrete time:  N(t+1) = N(t) + N(t) * 1.0 * ( ( K - N ) / K );
            if ( n_in_each_patch[i] < 2 )
                n_off_born[i] = 0; // selfing not allowed in this case
            else {
                dN = (double) n_in_each_patch[i];
                n_off_born[i] = (int) ( dN + ( dN * ( ( K - dN )/K ) ) ); /* leaving out demographic stochasticity here for now.
                                                                           Could be added by making K a function of time or throwing
                                                                           in some white noise */
                if ( n_off_born[i] < 0 )
                    n_off_born[i] = 0;
                
                newN += n_off_born[i];
            }
        }
    }
    
    off_genotypes = (short int*) palloc(POOL_GENOTYPES, ( sizeof(short int) * newN * 2 * nLOCI ) );
    
    ///*	off_x = (double*) palloc( ( sizeof(double) * newN ) );
    //#if (D == 2)
    //	off_y = (double*) palloc( ( sizeof(double) * newN ) );
    //#endif
    //*/
    // now we know how many offspring to produce in each patch
    // figure out fitnesses
    
    //fprintf(stderr, "\nstart = %G",randU());
    
    if ( GAMETE_PRODUCTION_MODE == 2 ) { // bag of genes
        if ( ((totalGenerationsElapsed % TS_SAMPLING_FREQUENCY) == 0 && totalGenerationsElapsed > 0 ) || ( (totalGenerationsElapsed+1) == (nMUTATIONS * nGENERATIONS_MAX) ) || ( RECORDING_TIMES_IN_FILE && totalGenerationsElapsed == nextRecordingTime ) ) {
            for ( j = 0; j < nPATCHES; j++ )
                fitSum[j] = 0.0;
            calculateFitnesses(fitnesses,fitSum); // think of fitnesses as the vector of individuals, but with each individual weighted by its fitness
            calcExpectedME(fitnesses, fitSum, gatherLDvalues);
        }
        bagOfGenes(off_genotypes, n_off_born, newN); // make the offspring genotypes
    }
    else {
        for ( j = 0; j < nPATCHES; j++ )
            fitSum[j] = 0.0;
        
        calculateFitnesses(fitnesses,fitSum); // think of fitnesses as the vector of individuals, but with each individual weighted by its fitness
        
        if ( ((totalGenerationsElapsed % TS_SAMPLING_FREQUENCY) == 0 && totalGenerationsElapsed > 0 ) || ( (totalGenerationsElapsed+1) == (nMUTATIONS * nGENERATIONS_MAX) ) || ( RECORDING_TIMES_IN_FILE && totalGenerationsElapsed == nextRecordingTime ) )
            calcExpectedME(fitnesses, fitSum, gatherLDvalues);
        
        /* // test
         fprintf(stderr, "\n\n");
         for ( i = 0; i < nPATCHES; i++ )
         fprintf(stderr, "%i, %G\n", n_in_each_patch[i], fitSum[i]);
         fprintf(stderr, "\n\n");
         for ( i = 0; i < N; i++ )
         fprintf(stderr, "%G, ", fitnesses[i]);
         fprintf(stderr, "\n\n");
         exit(0);
         // */ // end test
        if ( ((totalGenerationsElapsed % TS_SAMPLING_FREQUENCY) == 0  && totalGenerationsElapsed > 0) || ( RECORDING_TIMES_IN_FILE && totalGenerationsElapsed == nextRecordingTime )  ) {
            fprintf(effPopSizeData, "%li", totalGenerationsElapsed);
        }
        
        // now use fitnesses to make offspring
        firsti = 0;
        lasti = n_in_each_patch[0] - 1;
        if ( lasti < 0 )
            fprintf(stderr,"\nError in reproduction(); %i individuals in patch 0\n",n_in_each_patch[0]);
        offCount = 0;
        for ( j = 0; j < nPATCHES; j++ ) {
            no = n_off_born[j];
            nj = n_in_each_patch[j];
            
            if ( ((totalGenerationsElapsed % TS_SAMPLING_FREQUENCY) == 0  && totalGenerationsElapsed > 0) || ( RECORDING_TIMES_IN_FILE && totalGenerationsElapsed == nextRecordingTime ) ) {
                nReproducers = 0;
                for ( i = firsti; i <= lasti; i++ )
                    offPerParent[i] = 0.0;
            }
            
            if ( no > 0 ) {
                dN = ((double) nj) - 1.0;
                fs = fitSum[j];
                for ( i = 0; i < no; i++ ) {
                    if ( nj > 2 ) {
                        // we need two parents whose indexes are between firsti and lasti, inclusive
                        // choose first parent
                        if ( PURE_NEUTRAL_MODEL ) {
                            mom = firsti + (randI() % nj);
                            do {
                                dad = firsti + (randI() % nj);
                            } while (dad == mom);
                            // /* // test check
                            if ( dad > lasti || mom > lasti || dad < firsti || mom < firsti ) {
                                fprintf(stderr, "\nError in reproduce() in PURE_NEUTRAL_MODEL parent choice:\n\tIndexes out of bounds.\n\t");
                                fprintf(stderr, "firsti = %i, lasti = %i, dad = %i, mom = %i, j (patch) = %i, nj = %i\n", firsti, lasti, dad, mom, j, nj);
                                exit(1);
                            }
                        }
                        else {
                            dumd = randU();
                            fitVal = dumd * fs;
                            mom = ((int) (dumd * dN)) + firsti; // an approximate point to start looking in the fitnesses array
                            /* // test check,
                             if ( mom < firsti || mom > lasti ) {
                             fprintf(stderr, "\n\nERROR!  Mom out of bounds!\n\n");
                             exit(1);
                             }
                             // */ // end test check
                            while ( fitnesses[mom] < fitVal ) {
                                mom++;
                                // /* // test check
                                if ( mom > lasti ) {
                                    fprintf(stderr, "\n\nERROR!  Mom too big!\n\n");
                                    exit(1);
                                }
                                // */ // end test check
                            }
                            while ( (mom > firsti) && (fitnesses[(mom-1)] > fitVal) ) {
                                mom--;
                            }
                            
                            // choose second parent
                            do {
                                dumd = randU();
                                fitVal = dumd * fs; // relative position of random number in the cumulative distribution of this patch's fitness values.
                                dad = ((int) (dumd * dN)) + firsti;  // approximate place to start looking in the array
                                /* // test check
                                 if ( dad < firsti || dad > lasti ) {
                                 fprintf(stderr, "\n\nERROR!  Dad out of bounds!");
                                 exit(1);
                                 }
                                 // */ // end test check
                                while ( fitnesses[dad] < fitVal ) {
                                    dad++;
                                    // /* // test check
                                    if ( dad > lasti ) {
                                        fprintf(stderr, "\n\nERROR!  Dad too big!\n\n");
                                        exit(1);
                                    }
                                    // */ // end test check
                                }
                                while ( (dad > firsti) && (fitnesses[(dad-1)] > fitVal) ) {
                                    dad--;
                                }
                            } while ( dad == mom ); // end of do loop; condition makes sure selfing is NOT allowed
                            //fprintf(stderr, "dad and mom = \t%i\t%i\n", dad, mom);
                        }
                    }
                    
                    else if ( nj == 2 ) {  // when there are exactly two present in the patch, it's easy to choose.
                        if ( randU() < 0.5 ) {
                            mom = firsti;
                            dad = lasti;
                        }
                        else {
                            mom = lasti;
                            dad = firsti;
                        }
                    }
                    
                    else if ( nj == 1 ) { // selfing only allowed when FIXED_N == 1, otherwise we would't have entered this loop.
                        mom = firsti;
                        dad = firsti;
                    }
                    
                    // now mom and dad are established for this offspring:
                    //					if ( previousPatches[mom] != j )
                    //						nEffectiveMigrants[j] += 1.0;
                    //
                    //					if ( previousPatches[dad] != j )
                    //						nEffectiveMigrants[j] += 1.0;
                    
                    // now use the dad and mom to make a zygote
                    
                    /* // db
                     for ( ii = 0; ii < nLOCI; ii++ ) {
                     if ( variable_loci[ii] != 0 && variable_loci[ii] != 1 ) {
                     fprintf(stderr,"\n bad value of variable_loci[%i] = %i\n\tnLOCI = %li\n\n\tExiting\n\n", ii, variable_loci[ii], nLOCI);
                     exit(1);
                     }
                     }
                     // */ // end db
                    
                    if ( ((totalGenerationsElapsed % TS_SAMPLING_FREQUENCY) == 0 && totalGenerationsElapsed > 0) || ( RECORDING_TIMES_IN_FILE && totalGenerationsElapsed == nextRecordingTime ) ) {
                        if ( offPerParent[mom] < 1.0 )
                            nReproducers++;
                        if ( offPerParent[dad] < 1.0 )
                            nReproducers++;
                        offPerParent[mom] = offPerParent[mom] + 1.0;
                        offPerParent[dad] = offPerParent[dad] + 1.0;
                    }
                    
                    makeZygoteChromosomes(mom,(off_genotypes + (2 * nLOCI * offCount)));
                    makeZygoteChromosomes(dad,(off_genotypes + (2 * nLOCI * offCount) + 1));
                    
                    /* // db
                     for ( ii = 0; ii < nLOCI; ii++ ) {
                     if ( variable_loci[ii] != 0 && variable_loci[ii] != 1 ) {
                     fprintf(stderr,"\n bad value of variable_loci[%i] = %i\n\tnLOCI = %li\n\n\tExiting\n\n", ii, variable_loci[ii], nLOCI);
                     exit(1);
                     }
                     }
                     // */ // end db
                    
                    if ( !OFFSPRING_IN_RANDOM_LOCATIONS ) {
                        off_x[offCount] = x_locations[mom];
#if ( D == 2 )
                        off_y[offCount] = y_locations[mom];
#endif
                    }
                    offCount++;
                    // /* // test check
                    if ( offCount > (3 * N) ) {
                        fprintf(stderr, "\n\nERROR!  Too many offspring!!!");
                        exit(1);
                    }
                    // */ // end test check
                }
            }
            
            if ( ((totalGenerationsElapsed % TS_SAMPLING_FREQUENCY) == 0 && totalGenerationsElapsed > 0) || ( RECORDING_TIMES_IN_FILE && totalGenerationsElapsed == nextRecordingTime ) ) {
                if ( no > 0 )
                    fprintf(effPopSizeData, " %i %i %li %E %E", j, no, nReproducers, calcMeanMagnitude(&offPerParent[firsti], nj), calcVariance(&offPerParent[firsti], nj));
                else
                    fprintf(effPopSizeData, " 0 0 0 0.0 0.0");
                // mean should always be 2 offspring/parent unless population size is not fixed, so that is really just a check on the system
            }
            
            firsti = lasti + 1;
            if ( j < (nPATCHES-1) )
                lasti = lasti + n_in_each_patch[(j+1)];
        }
        if ( ((totalGenerationsElapsed % TS_SAMPLING_FREQUENCY) == 0 && totalGenerationsElapsed > 0) || ( RECORDING_TIMES_IN_FILE && totalGenerationsElapsed == nextRecordingTime ) )
            fprintf(effPopSizeData, "\n");
    }
    
    
    // copy from offspring array to parent array
    
    // re-size arrays if necessary
    if ( !FIXED_N ) {
        pfree(patch_locations);
        patch_locations = (int*) palloc(POOL_PATCH_LOCATIONS,  sizeof(int) * newN );
        pfree(x_locations);
        x_locations = (double*) palloc(POOL_X_LOCATIONS, sizeof(double) * newN );
#if ( D == 2 )
        pfree(y_locations);
        y_locations = (double*) palloc(POOL_Y_LOCATIONS, sizeof(double) * newN );
#endif
    }
    
    // deal with locations
    if ( GAMETE_PRODUCTION_MODE == 2 || OFFSPRING_IN_RANDOM_LOCATIONS ) {
        if ( !FIXED_N ) { /* here we can't just re-use the old locations since the numbers may differ; hence, we choose random
                           also, with the "bag of gametes", there really aren't parental locations */
            int row, col;
            double mult = 1.0 / ((double) PATCHES); // "mult" is a multiplier that gives the width of a patch on the unit line/square
            dpt1 = x_locations;
            ipt1 = patch_locations;
#if ( D == 2 )
            dpt3 = y_locations;
#endif
            for ( i = 0; i < nPATCHES; i++ ) {
                count = n_off_born[i];
                for ( j = 0; j < count; j++ ) {
#if ( D == 1 )
                    *dpt1 = (randU() * mult) + i; // random x coordinate in patch i
#elif ( D == 2 )
                    row = i % PATCHES;
                    col = i / PATCHES;
                    *dpt1 = (randU() * mult) + col;
                    *dpt3 = (randU() * mult) + row;
                    dpt3++;
#endif
                    *ipt1 = i;
                    dpt1++;
                    ipt1++;
                }
            }
            // the new random locations are now assigned
        }
        // no "else" is needed here because if FIXED_N == 1, then just use the same x and y locations
    }
    
    else { //  offspring are born where mom was
        dpt1 = x_locations;
        dpt2 = off_x;
        ipt1 = patch_locations;
#if ( D == 2 )
        dpt3 = y_locations;
        dpt4 = off_y;
#endif
        for ( i = 0; i < nPATCHES; i++ ) {
            no = n_off_born[i];
            for ( j = 0; j < no; j++ ) {
                *dpt1 = (*dpt2);
                *ipt1 = i;
                dpt1++;
                dpt2++;
                ipt1++;
#if ( D == 2 )
                *dpt3 = (*dpt4);
                dpt3++;
                dpt4++;
#endif
            }
        }
    }
    
    // copy the new genotypes
    pfree(genotypes);  // free the old block that is no longer needed
    genotypes = off_genotypes; // assign the global to the new array
    
    if ( !FIXED_N ) {
        N = newN;
        for ( i = 0; i < nPATCHES; i++ )
            n_in_each_patch[i] = n_off_born[i];
    }
    
    //	if ( GAMETE_PRODUCTION_MODE < 2 ) {
    //		for ( j = 0; j < nPATCHES; j++ )
    //			fprintf(effMigRates,"%G, ", (nEffectiveMigrants[j] / ( 2.0 * ((double) (n_off_born[j])))) );
    //		// former now converted to a rate; denom multiplied by 2 because each has 2 parents
    //		fprintf(effMigRates,"\n");
    //	}
}


void setUpFitnesses(void)
{
    int i, j, dum;
    int fitArrayIndex;
    size_t sizeNeeded;
    double incr;
    
    sizeNeeded = (size_t) (nLOCI * nPATCHES * sizeof(double));
    
    fit0 = (double *) malloc( sizeNeeded );
    fit1 = (double *) malloc( sizeNeeded );
    fitH = (double *) malloc( sizeNeeded );
    
    
    for ( i = 0; i < nLOCI; i++ )
        H[i] = HVAL; // could be changed later if there is different H for each locus
    
    memset(fit0, 0, sizeNeeded);
    memset(fit1, 0, sizeNeeded);
    memset(fitH, 0, sizeNeeded);
    
    if ( nPATCHES == 2 ) {
        for ( i = 0; i < nLOCI; i++ ) {
            if ( IS_SELECTED_LOCUS[i] ) {
                // patch 0, where 0 is the best allele
                fit0[(i * nPATCHES)] = S_MAX0[i];
                fit1[(i * nPATCHES)] = 0.0;
                fitH[(i * nPATCHES)] = ( 1.0 - H[i] ) * S_MAX0[i]; // H[i] is the dominance coefficient for 1 alleles
                
                //fprintf(stderr, "I'm here. locus = %i, fit0 = %G, fit1 = %G, fitH = %G\n",i,fit0[(i * nPATCHES)],fit1[(i * nPATCHES)],fitH[(i * nPATCHES)]);
                
                // patch 1, where 1 is the best allele
                fit0[((i * nPATCHES) + 1)] = 0.0;
                fit1[((i * nPATCHES) + 1)] = S_MAX1[i];
                fitH[((i * nPATCHES) + 1)] = H[i] * S_MAX1[i];
            }
        }
        epi_patch_multipliers[0] = -1.0;
        epi_patch_multipliers[1] = 1.0;
    }
    
    else if ( D == 2 && nPATCHES == 4 ) {
        for ( i = 0; i < nLOCI; i++ ) {
            if ( IS_SELECTED_LOCUS[i] ) {
                // patch 0, where 0 is the best allele
                fit0[(i * nPATCHES)] = S_MAX0[i];
                fit0[((i * nPATCHES) + 2)] = fit0[(i * nPATCHES)];
                fit1[(i * nPATCHES)] = 0.0;
                fit1[((i * nPATCHES) + 2)] = fit1[(i * nPATCHES)];
                fitH[(i * nPATCHES)] = ( 1.0 - H[i] ) * S_MAX0[i]; // H[i] is the dominance coefficient for 1 alleles
                fitH[((i * nPATCHES) + 2)] = fitH[(i * nPATCHES)];
                //fprintf(stderr, "I'm here. locus = %i, fit0 = %G, fit1 = %G, fitH = %G\n",i,fit0[(i * nPATCHES)],fit1[(i * nPATCHES)],fitH[(i * nPATCHES)]);
                
                // patch 1, where 1 is the best allele
                fit0[((i * nPATCHES) + 1)] = 0.0;
                fit0[((i * nPATCHES) + 3)] = fit0[((i * nPATCHES) + 1)];
                fit1[((i * nPATCHES) + 1)] = S_MAX1[i];
                fit1[((i * nPATCHES) + 3)] = fit1[((i * nPATCHES) + 1)];
                fitH[((i * nPATCHES) + 1)] = H[i] * S_MAX1[i];
                fitH[((i * nPATCHES) + 3)] = fitH[((i * nPATCHES) + 1)];
            }
        }
        epi_patch_multipliers[0] = -1.0;
        epi_patch_multipliers[1] = 1.0;
        epi_patch_multipliers[2] = -1.0;
        epi_patch_multipliers[3] = 1.0;
    }
    
    else if ( MOSAIC ) {
        // set up the random mosaic
        int count1 = 0, count0 = 0, dum;
        BEST_ALLELE_IN_PATCH[0] = 0;
        BEST_ALLELE_IN_PATCH[(nPATCHES-1)] = 1;
        for ( i = 1; i < (nPATCHES-1); i++ ) {
            if ( randU() < 0.5 ) {
                BEST_ALLELE_IN_PATCH[i] = 0;
                count0++;
            }
            else {
                BEST_ALLELE_IN_PATCH[i] = 1;
                count1++;
            }
        }
        while ( (count0 - count1) > 1 ) { // imbalance with too many allotted to 0
            dum = 1 + ( randI() % (nPATCHES-2) ); // choose random patch whose index is between 1 and nPATCHES-2 (leave 0 and nPATCHES-1 alone)
            if ( BEST_ALLELE_IN_PATCH[dum] == 0 ) { // check if it is allotted to zero
                BEST_ALLELE_IN_PATCH[dum] = 1;
                count0--;
                count1++;
            }
        }
        while ( (count1 - count0) > 1 ) { // imbalance with too many allotted to 1
            dum = randI() % nPATCHES; // choose random patch
            if ( BEST_ALLELE_IN_PATCH[dum] == 1 ) { // check if it is allotted to 1
                BEST_ALLELE_IN_PATCH[dum] = 0;
                count0++;
                count1--;
            }
        }
        
        for ( j = 0; j < nPATCHES; j++ ) {
            if ( BEST_ALLELE_IN_PATCH[j] ) { // 1 is best allele
                for ( i = 0; i < nLOCI; i++ ) {
                    if ( IS_SELECTED_LOCUS[i] ) {
                        fit1[((i * nPATCHES) + j)] = S_MAX1[i];
                        fit0[((i * nPATCHES) + j)] = 0.0;
                        fitH[((i * nPATCHES) + j)] = H[i] * S_MAX1[i];
                    }
                }
            }
            else { // 0 is best allele
                for ( i = 0; i < nLOCI; i++ ) {
                    if ( IS_SELECTED_LOCUS[i] ) {
                        fit1[((i * nPATCHES) + j)] = 0.0;
                        fit0[((i * nPATCHES) + j)] = S_MAX0[i];
                        fitH[((i * nPATCHES) + j)] = (1.0 - H[i]) * S_MAX0[i];
                    }
                }
            }
        }
        fprintf(stderr, "\nBest allele in each patch:\n");
        for ( i = 0; i < nPATCHES; i++ )
            fprintf(stderr, "%i, ", BEST_ALLELE_IN_PATCH[i]);
        fprintf(stderr, "\n");
        for ( i = 0; i < nPATCHES; i++ ) {
            if ( BEST_ALLELE_IN_PATCH[i] )
                epi_patch_multipliers[i] = 1.0;
            else
                epi_patch_multipliers[i] = -1.0;
        }
    }
    
    else {
        // gradient
        double step0, step1, maxH, HstepUp, HstepDown, avg;
        int maxHpatch, rowNum, colNum, nSteps, midPatch1, midPatch2;
        
        for ( i = 0; i < nLOCI; i++ ) {
            if ( IS_SELECTED_LOCUS[i] ) {
                
                // increment of change
                if ( D == 1) {
                    step0 = S_MAX0[i] / ( (double) (PATCHES - 1) );
                    step1 = S_MAX1[i] / ( (double) (PATCHES - 1) );
                    maxHpatch = (int) ( 0.5 + (((double) (PATCHES - 1)) * H[i]) );
                    if  ( H[i] < 1.0 && maxHpatch == (PATCHES - 1) )
                        maxHpatch--;
                    else if ( H[i] > 0.0 && maxHpatch == 0 )
                        maxHpatch++;
                }
                else if ( D == 2 ) {
                    step0 = S_MAX0[i] / ( (double) ( (PATCHES - 1) ) );
                    step1 = S_MAX1[i] / ( (double) ( (PATCHES - 1) ) );
                    maxHpatch = (int) ( 0.5 + (((double) (PATCHES - 1)) * H[i]) ); /* this is not a patch number
                                                                                    per se but rather the number of steps away from
                                                                                    the upper row */
                    if  ( H[i] < 1.0 && maxHpatch == ((PATCHES - 1)) )
                        maxHpatch--;
                    else if ( H[i] > 0.0 && maxHpatch == 0 )
                        maxHpatch++;
                }
                
                // upper row; best for 0
                fit0[((i * nPATCHES))] = S_MAX0[i];
                fit1[((i * nPATCHES))] = 0.0;
                fitH[((i * nPATCHES))] = ( 1.0 - H[i] ) * S_MAX0[i];
                
                // lower row; best for 1
                fit0[((i * nPATCHES) + nPATCHES - 1)] = 0.0;
                fit1[((i * nPATCHES) + nPATCHES - 1)] = S_MAX1[i];
                fitH[((i * nPATCHES) + nPATCHES - 1)] = H[i] * S_MAX1[i];
                
                maxH = ( (1.0 - H[i]) * S_MAX0[i] ) + ( H[i] * S_MAX1[i] );
                HstepUp = (maxH - fitH[(i * nPATCHES)]) / ( (double) maxHpatch );
                HstepDown = (maxH - fitH[((i * nPATCHES) + nPATCHES - 1)]) / ( (double) ((PATCHES - 1) - maxHpatch) );
                
                for ( j = 1; j < nPATCHES; j++ ) {
                    if ( D == 1 ) {
                        fit0[(i * nPATCHES)] = fit0[(i * nPATCHES)] - (step0 * ((double) j));
                        fit1[((i * nPATCHES) + j)] = fit1[(i * nPATCHES)] + (step1 * ((double) j));
                        if ( j <= maxHpatch )
                            fitH[((i * nPATCHES) + j)] = fitH[(i * nPATCHES)] + (HstepUp * ((double) j));
                        else
                            fitH[((i * nPATCHES) + j)] = fitH[((i * nPATCHES) + maxHpatch)] - (HstepDown * ((double) (j-maxHpatch)));
                    }
                    else if ( D == 2 ) {
                        rowNum = j % PATCHES;
                        nSteps = rowNum;
                        fit0[((i * nPATCHES) + j)] = fit0[(i * nPATCHES)] - (step0 * ((double) nSteps));
                        fit1[((i * nPATCHES) + j)] = fit1[(i * nPATCHES)] + (step1 * ((double) nSteps));
                        if ( nSteps <= maxHpatch )
                            fitH[((i * nPATCHES) + j)] = fitH[(i * nPATCHES)] + ( HstepUp * ((double) nSteps) );
                        else
                            fitH[((i * nPATCHES) + j)] = fitH[((i * nPATCHES) + maxHpatch)] - ( HstepUp * ((double) (nSteps - maxHpatch)) );
                    }
                }
                // handle the symmetric case for even number of patches in 1-D (i.e. when there technically is no "middle" patch)
                if ( (PATCHES > 2) && (D == 1) && ((PATCHES % 2) == 0) && (H[i] == 0.5) ) {
                    midPatch1 = PATCHES / 2;
                    midPatch2 = midPatch1 - 1;
                    avg = ( fitH[((i * nPATCHES) + midPatch1)] + fitH[((i * nPATCHES) + midPatch2)] ) / 2.0;
                    fitH[((i * nPATCHES) + midPatch1)] = avg;
                    fitH[((i * nPATCHES) + midPatch2)] = avg;
                }
            }
        }
        
        // set up epi_patch_multipliers
        if ( D == 1) {
            if ( (PATCHES % 2) == 0 ) { // even number of patches
                incr = 1.0 / ((double) (PATCHES / 2));
                midPatch1 = PATCHES / 2;
                epi_patch_multipliers[0] = -1.0;
                for ( i = 1; i < PATCHES; i++ ) {
                    epi_patch_multipliers[i] = epi_patch_multipliers[(i-1)] + incr;
                    if ( i == midPatch1 )
                        epi_patch_multipliers[i] = epi_patch_multipliers[i] + incr; // necessary since there is no middle patch
                }
            }
            else {
                incr = 1.0 / ((double) ((PATCHES - 1) / 2));
                epi_patch_multipliers[0] = -1.0;
                for ( i = 1; i < PATCHES; i++ ) {
                    epi_patch_multipliers[i] = epi_patch_multipliers[(i-1)] + incr;
                }
            }
        }
        else { // 2-D habitat
            if ( (PATCHES % 2) == 0 ) {
                incr = 1.0 / ((double) (PATCHES / 2));
                midPatch1 = PATCHES / 2;
            }
            else
                incr = 1.0 / ((double) ((PATCHES - 1) / 2));
            for ( i = 0; i < nPATCHES; i++ ) {
                rowNum = i % PATCHES;
                if ( (PATCHES % 2) == 0 ) {
                    if ( rowNum < midPatch1 )
                        epi_patch_multipliers[i] = -1.0 + ( ((double) rowNum) * incr );
                    else
                        epi_patch_multipliers[i] = -1.0 + ( ((double) (rowNum + 1)) * incr);
                }
                else
                    epi_patch_multipliers[i] = -1.0 + (rowNum * incr);
            }
            
        }
        
    }
    
    
    if ( MULTIPLICATIVE_FITNESS ) {  // adjust values for multiplying by a factor rather than adding selection coefficients in the additive case
        for ( i = 0; i < nLOCI; i++ ) {
            for ( j = 0; j < nPATCHES; j++ ) {
                fit0[((i * nPATCHES) + j)] = fit0[((i * nPATCHES) + j)] + 1.0;
                fit1[((i * nPATCHES) + j)] = fit1[((i * nPATCHES) + j)] + 1.0;
                fitH[((i * nPATCHES) + j)] = fitH[((i * nPATCHES) + j)] + 1.0;
            }
        }
    }
    
}


void setUpMap(void)
{
    int i, j, locus, count;
    
    if ( MAP_TYPE == 1 ) { // "infinite sites" random map
        int foo, diff, putHere, loci, startLocus;
        double spot, ml, meanDist;
        
        nLOCI = 2 * nCHROMOSOMES;
        
        allocateGlobals();
        
        for ( i = 0; i < nCHROMOSOMES; i++ ) {
            MAP_LENGTHS[i] = TOTAL_MAP_LENGTH / ((double) nCHROMOSOMES);
            LOCI_PER_CHROMOSOME[i] = 2;
            MAP[(2*i)] = 0.0;
            MAP[((2*i) + 1)] = MAP_LENGTHS[i];
        }
        
        for ( i = 0; i < nLOCI; i++ ) {
            is_reference_locus[i] = 1;
            IS_SELECTED_LOCUS[i] = 0;
            locusID[i] = -(i+1); // all reference loci assigned negative numbers
            is_reversed_locus[i] = 0;
        }
    }
    
    else if ( MAP_TYPE == 2 ) {
        loadMapFromFile();
        fprintf(stderr, "\n\nError!  Function not yet defined!! Exiting\n\n");
        exit(1);
    }
    
    else {
        // default map, in cM, based on previous simulations
        fprintf(stderr, "\n\nError!  MAP_TYPE Not defined!! Exiting\n\n");
        exit(1);
    }
    
    
    j = LOCI_PER_CHROMOSOME[0];
    int chrom = 0;
    for ( i = 0; i < nLOCI; i++ ) {
        if ( i == j ) {
            chrom++;
            j += LOCI_PER_CHROMOSOME[chrom];
        }
        chromosomeMembership[i] = chrom;
        //printf("chrom member %i\n", chromosomeMembership[i]);
    }
    
}


void setUpMutationSequence(void)
{
    int j, count_n, candidate;
    long int i;
    
    MUTATION_ORDER = (double *) malloc( sizeof(double) * nMUTATIONS );
    MUTATION_TYPE_SEQUENCE = (int *) malloc( sizeof(int) * nMUTATIONS );
    MUTATION_LOCATIONS = (int *) malloc( sizeof(int) * nMUTATIONS );
    H_SEQUENCE = (double *) malloc( sizeof(double) * nMUTATIONS );
    REVERSAL_SEQUENCE = (_Bool *) malloc( sizeof(_Bool) * nMUTATIONS );
    
    if (USE_MUTATIONS_FROM_FILES) {
        // read in order of loci affected and values of S_MAX
        fprintf(stderr, "\n\nERROR!  USE_MUTATIONS_FROM_FILES algorithm not defined! \n\n\tExiting from function setUpMutationSequence()...\n\n");
        exit(1);
    }
    else {
        // make sequence of loci affected and associated values of S_MAX
        _Bool keepLooking = 1;
        
        for ( i = 0; i < nLOCI; i++ ) {
            S_MAX1[i] = 0.0;
            S_MAX0[i] = 0.0;
        }
        
        for ( i = 0; i < nMUTATIONS; i++ ) {
            MUTATION_ORDER[i] = randU() * TOTAL_MAP_LENGTH;
            MUTATION_LOCATIONS[i] = randI() % nPATCHES;
            if ( ( i < NUMBER_BIG_TO_ESTABLISH ) && ( START_WITH_BIG_IN_PLACE ) )
                MUTATION_TYPE_SEQUENCE[i] = 1;
            else {
                if ( randU() < FRACTION_SELECTED_LOCI )
                    MUTATION_TYPE_SEQUENCE[i] = 1;
                else
                    MUTATION_TYPE_SEQUENCE[i] = 0;
            }
            
            if ( randU() < PROBABILITY_DERIVED_REVERSED )
                REVERSAL_SEQUENCE[i] = 1;
            else
                REVERSAL_SEQUENCE[i] = 0;
            
            if ( MUTATION_TYPE_SEQUENCE[i] ) {
                if ( MUTATION_DISTRIBUTION == 3 ) {
                    S_MAX1_SEQUENCE[i] = DEME1_CONSTANT_S;
                    H_SEQUENCE[i] = HVAL;
                }
                else {
                    if ( ( i < NUMBER_BIG_TO_ESTABLISH ) && ( START_WITH_BIG_IN_PLACE ) )
                        S_MAX1_SEQUENCE[i] = BIG_START_S_VAL;
                    else
                        S_MAX1_SEQUENCE[i] = (log(1.0 - randU())) * (-MEAN_S); // negative exponential distribution of S sizes
                    H_SEQUENCE[i] = HVAL;
                }
            }
            else {
                S_MAX1_SEQUENCE[i] = 0.0;
                H_SEQUENCE[i] = HVAL;
            }
        }
        
        double *sortedS, valAtPercent, width, maxVal;
        int numInTail, numToFlatten;
        
        if ( MUTATION_DISTRIBUTION == 1 ) { // fatten tail of mutations of large effect if desired
            sortedS = (double *) malloc(sizeof(double) * nMUTATIONS);
            // first, determine what the actual upper FATTEN_TAIL_PROPORTION of the distribution is
            numInTail = (int) ( (FATTEN_TAIL_PROPORTION * ((double) nMUTATIONS) * FRACTION_SELECTED_LOCI ) + 0.5 );
            for ( i = 0; i < nMUTATIONS; i++ )
                sortedS[i] = S_MAX1_SEQUENCE[i];
            qsort(sortedS, nMUTATIONS, sizeof(double), compare_doubles);
            valAtPercent = sortedS[(nMUTATIONS - numInTail)];
            width = FATTEN_TAIL_MAX - valAtPercent;
            maxVal = sortedS[(nMUTATIONS - 1)];
            fprintf(stderr, "\nvalAtPercent = %G\n",valAtPercent);
            // now fatten the tail
            if ( width > 0.0 && maxVal < FATTEN_TAIL_MAX ) {
                for ( i = 0; i < nMUTATIONS; i++ ) {
                    if ( S_MAX1_SEQUENCE[i] >= valAtPercent ) {
                        S_MAX1_SEQUENCE[i] = ( randU() * width ) + valAtPercent; // uniform on interval [valAtPercent,1.5]
                        fprintf(stderr, "\nChanged S_MAX1_SEQUENCE[%li] to %E",i,S_MAX1_SEQUENCE[i]);
                    }
                }
            }
        }
        //		else if ( MUTATION_DISTRIBUTION == 2 ) { // flatten distribution to have more large effect mutations early on
        //			// first, determine what the actual upper FATTEN_TAIL_PROPORTION of the distribution is
        //			numToFlatten = (int) ( (FATTEN_TAIL_PROPORTION * ((double) nMUTATIONS)) + 0.5 );
        //			fprintf(stderr, "\nnum to flatten = %i\n",numToFlatten);
        //			for ( i = 0; i < numToFlatten; i++ )
        //				S_MAX1[(MUTATION_ORDER[i])] = randU() * FATTEN_TAIL_MAX; // uniform random on [0,FATTEN_TAIL_MAX]
        //		}
        
        /* // test check
         if ( MUTATION_DISTRIBUTION ) {
         fprintf(stderr, "\n\nThe altered distribution of S_MAX1:\n");
         for ( i = 0; i < nLOCI; i++ )
         fprintf(stderr, "%G ",S_MAX1[i]);
         fprintf(stderr, "\n\n");
         }
         //exit(0);
         // */ //end test check
        if ( SYMMETRIC_MUTATIONS ) {
            for ( i = 0; i < nMUTATIONS; i++ )
                S_MAX0_SEQUENCE[i] = S_MAX1_SEQUENCE[i];
            if ( MUTATION_DISTRIBUTION == 1 )
                free(sortedS);
        }
        else {
            if ( MUTATION_DISTRIBUTION == 3 ) {
                for ( i = 0; i < nMUTATIONS; i++ ) {
                    if ( MUTATION_TYPE_SEQUENCE[i] )
                        S_MAX0_SEQUENCE[i] = DEME0_CONSTANT_S;
                    else
                        S_MAX0_SEQUENCE[i] = 0.0;
                }
            }
            else {
                for ( i = 0; i < nMUTATIONS; i++ ) {
                    if ( MUTATION_TYPE_SEQUENCE[i] )
                        S_MAX0_SEQUENCE[i] = (log(1.0 - randU())) * (-MEAN_S);
                    else
                        S_MAX0_SEQUENCE[i] = 0.0;
                }
                if ( MUTATION_DISTRIBUTION == 1 ) { // fatten tail of mutations of large effect if desired; repeat as above to keep things fair
                    for ( i = 0; i < nMUTATIONS; i++ )
                        sortedS[i] = S_MAX0_SEQUENCE[i];
                    qsort(sortedS, nMUTATIONS, sizeof(double), compare_doubles);
                    valAtPercent = sortedS[(nMUTATIONS - numInTail)];
                    width = FATTEN_TAIL_MAX - valAtPercent;
                    maxVal = sortedS[(nMUTATIONS - 1)];
                    //fprintf(stderr, "\nvalAtPercent = %G\n",valAtPercent);
                    // now fatten the tail
                    if ( width > 0.0 && maxVal < FATTEN_TAIL_MAX ) {
                        for ( i = 0; i < nMUTATIONS; i++ ) {
                            if ( S_MAX0_SEQUENCE[i] >= valAtPercent ) {
                                S_MAX0_SEQUENCE[i] = ( randU() * width ) + valAtPercent; // uniform on interval [valAtPercent,FATTEN_TAIL_MAX]
                                //fprintf(stderr, "\nChanged S_MAX0[%i] to %G",i,S_MAX0[i]);
                            }
                        }
                    }
                    free(sortedS);
                }
            }
            //			else if ( MUTATION_DISTRIBUTION == 2 ) { // flatten distribution to have more large effect mutations early on
            //				// first, determine what the actual upper FATTEN_TAIL_PROPORTION of the distribution is
            //				for ( i = 0; i < numToFlatten; i++ )
            //					S_MAX0[(MUTATION_ORDER[i])] = randU() * FATTEN_TAIL_MAX; // uniform random on [0,FATTEN_TAIL_MAX]
            //			}
        }
    }
    
    // MUTATION_ORDER now specifies which loci are under selection, and in which order those loci will mutate from 0 to 1
    // make some matlab scripts for easy loading of records
    
    FILE *fpt;
    fpt = fopen("Mutations.m","w");
    fprintf(fpt,"MUTATION_ORDER = [");
    for ( i = 0; i < nMUTATIONS; i++ )
        fprintf(fpt, "%E ", MUTATION_ORDER[i]);
    fprintf(fpt,"];\n");
    fprintf(fpt,"MUTATION_TYPE_SEQUENCE = [");
    for ( i = 0; i < nMUTATIONS; i++ )
        fprintf(fpt, "%i ", MUTATION_TYPE_SEQUENCE[i]);
    fprintf(fpt,"];\n");
    fprintf(fpt,"MUTATION_LOCATIONS = [");
    for ( i = 0; i < nMUTATIONS; i++ )
        fprintf(fpt, "%i ", MUTATION_LOCATIONS[i]);
    fprintf(fpt,"];\n");
    fprintf(fpt,"MUTATION_S_VALUES_0 = [");
    for ( i = 0; i < nMUTATIONS; i++ )
        fprintf(fpt, "%G ", S_MAX0_SEQUENCE[i]);
    fprintf(fpt,"];\n");
    fprintf(fpt,"MUTATION_S_VALUES_1 = [");
    for ( i = 0; i < nMUTATIONS; i++ )
        fprintf(fpt, "%G ", S_MAX1_SEQUENCE[i]);
    fprintf(fpt,"];\n");
    fprintf(fpt,"REVERSAL_SEQUENCE = [");
    for ( i = 0; i < nMUTATIONS; i++ )
        fprintf(fpt, "%i ", REVERSAL_SEQUENCE[i]);
    fprintf(fpt,"];\n");
    fclose(fpt);
    
    
}






/*  ***************************************************** */
/*  **************   UTILITY FUNCTIONS   **************** */
/*  ***************************************************** */


void adjustFitnesses(void)
{
    int j, dum, i;
    long int locus = newestLocus;
    
    if ( nPATCHES == 2 ) {
        // patch 0, where 0 is the best allele
        fit0[(locus * nPATCHES)] = S_MAX0[locus];
        fit1[(locus * nPATCHES)] = 0.0;
        fitH[(locus * nPATCHES)] = ( 1.0 - H[locus] ) * S_MAX0[locus]; // H[locus] is the dominance coefficient for 1 alleles
        
        //fprintf(stderr, "I'm here. locus = %i, fit0 = %G, fit1 = %G, fitH = %G\n",locus,fit0[(locus * nPATCHES)],fit1[(locus * nPATCHES)],fitH[(locus * nPATCHES)]);
        
        // patch 1, where 1 is the best allele
        fit0[((locus * nPATCHES) + 1)] = 0.0;
        fit1[((locus * nPATCHES) + 1)] = S_MAX1[locus];
        fitH[((locus * nPATCHES) + 1)] = H[locus] * S_MAX1[locus];
    }
    
    else if ( D == 2 && nPATCHES == 4 ) {
        i = locus;
        fit0[(i * nPATCHES)] = S_MAX0[i];
        fit0[((i * nPATCHES) + 2)] = fit0[(i * nPATCHES)];
        fit1[(i * nPATCHES)] = 0.0;
        fit1[((i * nPATCHES) + 2)] = fit1[(i * nPATCHES)];
        fitH[(i * nPATCHES)] = ( 1.0 - H[i] ) * S_MAX0[i]; // H[i] is the dominance coefficient for 1 alleles
        fitH[((i * nPATCHES) + 2)] = fitH[(i * nPATCHES)];
        fit0[((i * nPATCHES) + 1)] = 0.0;
        fit0[((i * nPATCHES) + 3)] = fit0[((i * nPATCHES) + 1)];
        fit1[((i * nPATCHES) + 1)] = S_MAX1[i];
        fit1[((i * nPATCHES) + 3)] = fit1[((i * nPATCHES) + 1)];
        fitH[((i * nPATCHES) + 1)] = H[i] * S_MAX1[i];
        fitH[((i * nPATCHES) + 3)] = fitH[((i * nPATCHES) + 1)];
    }
    
    else if ( MOSAIC ) {
        // mosaic
        for ( j = 0; j < nPATCHES; j++ ) {
            if ( BEST_ALLELE_IN_PATCH[locus] ) { // 1 is best allele
                fit1[((locus * nPATCHES) + j)] = S_MAX1[locus];
                fit0[((locus * nPATCHES) + j)] = 0.0;
                fitH[((locus * nPATCHES) + j)] = H[locus] * S_MAX1[locus];
            }
            else { // 0 is best allele
                fit1[((locus * nPATCHES) + j)] = 0.0;
                fit0[((locus * nPATCHES) + j)] = S_MAX0[locus];
                fitH[((locus * nPATCHES) + j)] = (1.0 - H[locus]) * S_MAX0[locus];
            }
        }
    }
    
    else {
        // gradient
        double step0, step1, maxH, HstepUp, HstepDown, avg;
        int maxHpatch, rowNum, colNum, nSteps, midPatch1, midPatch2;
        
        // increment of change
        if ( D == 1) {
            step0 = S_MAX0[locus] / ( (double) (PATCHES - 1) );
            step1 = S_MAX1[locus] / ( (double) (PATCHES - 1) );
            maxHpatch = (int) ( 0.5 + (((double) (PATCHES - 1)) * H[locus]) );
            if  ( H[locus] < 1.0 && maxHpatch == (PATCHES - 1) )
                maxHpatch--;
            else if ( H[locus] > 0.0 && maxHpatch == 0 )
                maxHpatch++;
        }
        else if ( D == 2 ) {
            step0 = S_MAX0[locus] / ( (double) ( (PATCHES - 1) ) );
            step1 = S_MAX1[locus] / ( (double) ( (PATCHES - 1) ) );
            maxHpatch = (int) ( 0.5 + (((double) (PATCHES - 1)) * H[locus]) ); /* this is not a patch number
                                                                                per se but rather the number of steps away from
                                                                                the upper row */
            if  ( H[locus] < 1.0 && maxHpatch == ( (PATCHES - 1)) )
                maxHpatch--;
            else if ( H[locus] > 0.0 && maxHpatch == 0 )
                maxHpatch++;
        }
        
        // upper row; best for 0
        fit0[(locus * nPATCHES)] = S_MAX0[locus];
        fit1[(locus * nPATCHES)] = 0.0;
        fitH[(locus * nPATCHES)] = ( 1.0 - H[locus] ) * S_MAX0[locus];
        
        // lower row; best for 1
        fit0[((locus * nPATCHES) + nPATCHES - 1)] = 0.0;
        fit1[((locus * nPATCHES) + nPATCHES - 1)] = S_MAX1[locus];
        fitH[((locus * nPATCHES) + nPATCHES - 1)] = H[locus] * S_MAX1[locus];
        
        maxH = ( (1.0 - H[locus]) * S_MAX0[locus] ) + ( H[locus] * S_MAX1[locus] );
        HstepUp = (maxH - fitH[(locus * nPATCHES)]) / ( (double) maxHpatch );
        HstepDown = (maxH - fitH[((locus * nPATCHES) + nPATCHES - 1)]) / ( (double) ((PATCHES - 1) - maxHpatch) );
        
        for ( j = 1; j < nPATCHES; j++ ) {
            if ( D == 1 ) {
                fit0[((locus * nPATCHES) + j)] = fit0[(locus * nPATCHES)] - (step0 * ((double) j));
                fit1[((locus * nPATCHES) + j)] = fit1[(locus * nPATCHES)] + (step1 * ((double) j));
                if ( j <= maxHpatch )
                    fitH[((locus * nPATCHES) + j)] = fitH[(locus * nPATCHES)] + (HstepUp * ((double) j));
                else
                    fitH[((locus * nPATCHES) + j)] = fitH[((locus * nPATCHES) + maxHpatch)] - (HstepDown * ((double) (j-maxHpatch)));
            }
            else if ( D == 2 ) {
                rowNum = j % PATCHES;
                nSteps = rowNum;
                fit0[((locus * nPATCHES) + j)] = fit0[(locus * nPATCHES)] - (step0 * ((double) nSteps));
                fit1[((locus * nPATCHES) + j)] = fit1[(locus * nPATCHES)] + (step1 * ((double) nSteps));
                if ( nSteps <= maxHpatch )
                    fitH[((locus * nPATCHES) + j)] = fitH[(locus * nPATCHES)] + ( HstepUp * ((double) nSteps) );
                else
                    fitH[((locus * nPATCHES) + j)] = fitH[((locus * nPATCHES) + maxHpatch)] - ( HstepUp * ((double) (nSteps - maxHpatch)) );
            }
        }
        // handle the symmetric case for even number of patches in 1-D (i.e. when there technically is no "middle" patch)
        if ( (PATCHES > 2) && (D == 1) && ((PATCHES % 2) == 0) && (H[locus] == 0.5) ) {
            midPatch1 = PATCHES / 2;
            midPatch2 = midPatch1 - 1;
            avg = ( fitH[((locus * nPATCHES) + midPatch1)] + fitH[((locus * nPATCHES) + midPatch2)] ) / 2.0;
            fitH[((locus * nPATCHES) + midPatch1)] = avg;
            fitH[((locus * nPATCHES) + midPatch2)] = avg;
        }
    }
    
    
    if ( MULTIPLICATIVE_FITNESS ) {  // adjust values for multiplying by a factor rather than adding selection coefficients in the additive case
        for ( j = 0; j < nPATCHES; j++ ) {
            fit0[((locus * nPATCHES) + j)] = fit0[((locus * nPATCHES) + j)] + 1.0;
            fit1[((locus * nPATCHES) + j)] = fit1[((locus * nPATCHES) + j)] + 1.0;
            fitH[((locus * nPATCHES) + j)] = fitH[((locus * nPATCHES) + j)] + 1.0;
        }
    }
}




void allocateGlobals(void)
{
    // other globals
    
    allele_frequencies = (double *) malloc( sizeof(double) * nLOCI ); // frequency of "1" allele at each locus across the whole population
    MAP = (double *) malloc( sizeof(double) * nLOCI ); // genetic map of all chromosomes; each zero is the start of a new chromosome
    S_MAX1 = (double *) malloc( sizeof(double) * nLOCI ); // values of selection coefficients for 1 alleles
    S_MAX0 = (double *) malloc( sizeof(double) * nLOCI ); // values of selection coefficients for 0 alleles; may be symmetric with S_MAX1
    S_MAX1_SEQUENCE = (double *) malloc( sizeof(double) * nMUTATIONS );
    S_MAX0_SEQUENCE = (double *) malloc( sizeof(double) * nMUTATIONS );
    IS_SELECTED_LOCUS = (int *) malloc( sizeof(int) * nLOCI );
    variable_loci = (int *) malloc( sizeof(int) * nLOCI ); // record which loci have more than one allele at a given time
    fixed_allele = (short int *) malloc( sizeof(short int) * nLOCI ); // will only be used for loci that are NOT variable
    chromosomeMembership = (int *) malloc( sizeof(int) * nLOCI );
    H = (double *) malloc( sizeof(double) * nLOCI );
    LOCI_PER_CHROMOSOME = (int*) malloc( sizeof(int) * nCHROMOSOMES );
    MAP_LENGTHS = (double*) malloc( sizeof(double) * nCHROMOSOMES );
    is_reference_locus = (int *) malloc( sizeof(int) * nLOCI );
    locusID = (long int *) malloc( sizeof(long int) * nLOCI);
    is_reversed_locus = (_Bool *) malloc( (sizeof(_Bool) * nLOCI) );
    if ( CONSIDER_EPISTASIS ) {
        any_epis_for_this_locus = (_Bool *) malloc( (sizeof(_Bool) * nLOCI) );
        epi_coeff_matrix = (double *) palloc( POOL_EPI_COEFFS, (sizeof(double) * (nLOCI * nLOCI)) );
        epistasisMatrix = (short int *) palloc( POOL_EPISTASIS, (sizeof(short int) * ( nLOCI * nLOCI )) );
    }
    
}


double boxMuller(double mu, double sd)
{
    /* boxmuller.c           Implements the Polar form of the Box-Muller
     Transformation
     
     (c) Copyright 1994, Everett F. Carter Jr.
     Permission is granted by the author to use
     this software for any application provided this
     copyright notice is preserved.
     [ smf's note: accessed online 12.03.06 at
     http://www.taygeta.com/random/boxmuller.html ]
     
     */
    
    double x1, x2, w, y1;
    static double y2;
    static int use_last = 0;
    
    if (use_last)		        /* use value from previous call */
    {
        y1 = y2;
        use_last = 0;
    }
    else
    {
        do {
            x1 = 2.0 * ( randU() ) - 1.0;
            x2 = 2.0 * ( randU() ) - 1.0;
            w = x1 * x1 + x2 * x2;
        } while ( w >= 1.0 );
        
        w = sqrt( (-2.0 * log( w ) ) / w );
        y1 = x1 * w;
        y2 = x2 * w;
        use_last = 1;
    }
    
    return ( mu + y1 * sd );
}


void calcDXY(int nSamples)
{
    // this function calculates Dxy using actual samples of sequences and counting the differences among them.
    int i;
    
    for ( i = 0; i < nPATCHES; i++ ) {
        if ( n_in_each_patch[i] <= nSamples ) {
            nSamples = n_in_each_patch[i] / 2;
        }
    }
    
    int j, xi, yj, locus;
    int deme0sample[nSamples], deme1sample[nSamples];
    short int *gtpt0, *g0, *g1;
    int haplo0[nSamples], haplo1[nSamples];
    double dum, mult, DXY[3];
    long int sampleDiffs, neutralDiffs, selectedDiffs;
    
    // choose haplotypes of diploid individuals at random
    for ( i = 0; i < nSamples; i++ ) {
        dum = randU();
        if ( dum < 0.25 ) {
            haplo0[i] = 0;
            haplo1[i] = 0;
        }
        else if ( dum < 0.5 ) {
            haplo0[i] = 0;
            haplo1[i] = 1;
        }
        else if ( dum < 0.75 ) {
            haplo0[i] = 1;
            haplo1[i] = 0;
        }
        else {
            haplo0[i] = 1;
            haplo1[i] = 1;
        }
    }
    
    chooseDemeSample(deme0sample, 0, nSamples);
    chooseDemeSample(deme1sample, 1, nSamples);
    
    DXY[0] = 0.0; // all sites
    DXY[1] = 0.0; // selected sites only
    DXY[2] = 0.0; // neutral sites only
    
    for ( i = 0; i < nSamples; i++ ) {
        xi = deme0sample[i];
        gtpt0 = genotypes + (2 * nLOCI * xi) + haplo0[i];
        for ( j = 0; j < nSamples; j++ ) {
            sampleDiffs = 0;
            selectedDiffs = 0;
            neutralDiffs = 0;
            yj = deme1sample[j];
            g0 = gtpt0;
            g1 = genotypes + (2 * nLOCI * yj) + haplo1[j];
            for ( locus = 0; locus < nLOCI; locus++ ) {
                if ( (*g0) != (*g1) ) {
                    sampleDiffs++;
                    if ( IS_SELECTED_LOCUS[locus] )
                        selectedDiffs++;
                    else
                        neutralDiffs++;
                }
                g0 += 2;
                g1 += 2;
            }
            DXY[0] = DXY[0] + sampleDiffs;
            DXY[1] = DXY[1] + selectedDiffs;
            DXY[2] = DXY[2] + neutralDiffs;
        }
    }
    
    mult = 1.0 / (nSamples * nSamples); // frequency of each haplotype pair assuming each one is unique
    DXY[0] = (mult * DXY[0]) / ((double)( nSELECTED_LOCI + nNEUTRAL_LOCI )); // DXY for all sites
    DXY[1] = (mult * DXY[1]) / ((double) nSELECTED_LOCI); // DXY if only selected sites were considered
    DXY[2] = (mult * DXY[2]) / ((double) nNEUTRAL_LOCI); // DXY if only neutral sites were considered
    
    fprintf(DXYTS, "%li %E %E %E\n", totalGenerationsElapsed, DXY[0], DXY[1], DXY[2]);
    
}


void calcDXY2(double *alleleFrequenciesByPatch)
{
    int i, j, counter, locus, lcount = 0, selcount = 0, neutcount = 0;
    double p1, p2, q1, q2, eDiff;
    int nPairs;
    _Bool yesSel;
    nPairs = (nPATCHES * (nPATCHES - 1) / 2);
    
    double expectedDiffs[nPairs];
    double expectedSelDiffs[nPairs];
    double expectedNeutDiffs[nPairs];
    for ( i = 0; i < nPairs; i++ ) {
        expectedDiffs[i] = 0.0;
        expectedSelDiffs[i] = 0.0;
        expectedNeutDiffs[i] = 0.0;
    }
    
    for ( locus = 0; locus < nLOCI; locus++ ) {
        if ( variable_loci[locus] ) {
            lcount++;
            if ( IS_SELECTED_LOCUS[locus] ) {
                selcount++;
                yesSel = 1;
            }
            else {
                neutcount++;
                yesSel = 0;
            }
            
            counter = 0;
            for ( i = 0; i < (nPATCHES-1); i++ ) {
                for ( j = (i + 1); j < nPATCHES; j++ ) {
                    p1 = *(alleleFrequenciesByPatch + (locus * nPATCHES) + i);
                    q1 = 1.0 - p1;
                    p2 = *(alleleFrequenciesByPatch + (locus * nPATCHES) + j);
                    q2 = 1.0 - p2;
                    eDiff = (p1 * q2) + (p2 * q1);
                    
                    expectedDiffs[counter] += eDiff;
                    if ( yesSel )
                        expectedSelDiffs[counter] += eDiff;
                    else
                        expectedNeutDiffs[counter] += eDiff;
                    
                    counter++;
                }
            }
            
        }
        
        
    }
    
    fprintf(DXYTS, "%li", totalGenerationsElapsed);
    for ( i = 0; i < nPairs; i++ ) {
        expectedDiffs[i] = expectedDiffs[i] / ((double) lcount);
        expectedSelDiffs[i] = expectedSelDiffs[i] / ((double) selcount);
        expectedNeutDiffs[i] = expectedNeutDiffs[i] / ((double) neutcount);
        fprintf(DXYTS, ",%E,%E,%E", expectedDiffs[i], expectedSelDiffs[i], expectedNeutDiffs[i]);
    }
    fprintf(DXYTS, "\n");
    
}



void calcExpectedME(double *fitpt, double *fitsumpt, int gatherLDvalues)
{
    int i, j, count = 0, nj, k, nImmigrants, nResidents;
    int *pppt = previousPatches; // previousPatches points to array of length N storing previous locations
    double immFitSum, *fvpt = fitpt, eme[nPATCHES], *fspt = fitsumpt, avgImmFit[nPATCHES], avgResFit[nPATCHES];
    double focalFitVal, maxFit[nPATCHES], randomFit[nPATCHES], minFit[nPATCHES], avgPatchFitness[nPATCHES];
    _Bool firstOne, allPatchesDiverged, stillLooking;
    double *fitTSvals;
    
    if ( PURE_NEUTRAL_MODEL ) {
        for ( i = 0; i < nPATCHES; i++ ) {
            
            nImmigrants = migrationCount[i];
            nResidents = n_in_each_patch[i] - migrationCount[i];
            eme[i] = ((double) nImmigrants) / ((double) n_in_each_patch[i]);
            
            if ( nImmigrants > 0 )
                avgImmFit[i] = 1.0;
            else
                avgImmFit[i] = -1.0;
            
            
            if ( nResidents > 0 )
                avgResFit[i] = 1.0;
            else
                avgResFit[i] = -1.0;
            
            if ( n_in_each_patch[i] > 0 )
                avgPatchFitness[i] = 1.0;
            else
                avgPatchFitness[i] = -1.0;
        }
    }
    else {
        for ( i = 0; i < nPATCHES; i++ ) {
            nj = n_in_each_patch[i];
            immFitSum = 0.0;
            fitTSvals = fvpt;
            for ( j = 0; j < nj; j++ ) {
                if ( (*pppt) != i ) { // this is a migrant
                    if ( j > 0 )
                        focalFitVal = ((*fvpt) - (*(fvpt-1))); // because the fitnesses vector hold cumulative sums
                    else
                        focalFitVal = (*fvpt);
                    immFitSum += focalFitVal;
                    if ( focalFitVal < 0.0 ) {
                        fprintf(stderr,"\nFit patch %i, individual %i = %E\nAdjacent fitness values:  ", i, j, focalFitVal);
                        for ( k = -1; k < 2; k++ ) {
                            fprintf(stderr,"%E  ", (*(fvpt+k)) );
                        }
                        fprintf(stderr, "\n");
                    }
                }
                fvpt++;
                pppt++;
            }
            // some data recording
            if ( RECORD_FIT_TS && nj >= nRECORD_FIT) { // want nRECORD_FIT recorded points
                int stepSize, lastOne, firstOne;
                if ( i == 0 ) {
                    stepSize = nj / nRECORD_FIT; // want nRECORD_FIT recorded points
                    lastOne = nj - 1; // want to sample as many immigrants as possible, and they should be near the end
                    firstOne = lastOne - ( (nRECORD_FIT - 1) * stepSize );
                    if ( firstOne < 0 ) {
                        fprintf(stderr, "\n Your algorithm for deme0 fitness sampling sucks.\n");
                        exit(1);
                    }
                    fitTSvals += firstOne;
                    fprintf(fitTSdeme0, "%li", totalGenerationsElapsed);
                    for ( j = firstOne; j <= lastOne; j += stepSize ) {
                        if ( j > 0 )
                            fprintf( fitTSdeme0, " %E", ((*fitTSvals) - (*(fitTSvals-1))) );
                        else
                            fprintf( fitTSdeme0, " %E", (*fitTSvals) );
                        fitTSvals += stepSize;
                    }
                }
                else if ( i == (nPATCHES-1) ) {
                    stepSize = nj / nRECORD_FIT; // want nRECORD_FIT recorded points
                    firstOne = 0; // want to sample as many immigrants as possible, and they should be near the beginning, I think...
                    lastOne = ( (nRECORD_FIT - 1) * stepSize );
                    if ( lastOne >= nj ) {
                        fprintf(stderr, "\n Your algorithm for deme1 fitness sampling sucks.\n");
                        exit(1);
                    }
                    fprintf(fitTSdeme1, "%li", totalGenerationsElapsed);
                    for ( j = firstOne; j <= lastOne; j += stepSize ) {
                        if ( j > 0 )
                            fprintf( fitTSdeme1, " %E", ((*fitTSvals) - (*(fitTSvals-1))) );
                        else
                            fprintf( fitTSdeme1, " %E", (*fitTSvals) );
                        fitTSvals += stepSize;
                    }
                }
            }
            if ( immFitSum > 0.0 && (*fspt) > 0.0 )
                eme[i] = immFitSum / (*fspt); /*  the expected effective migration rate is the proportion of all fitness
                                               in the patch that is owned by migrants */
            else
                eme[i] = 0.0;
            
            nImmigrants = migrationCount[i];
            if ( nImmigrants > 0 )
                avgImmFit[i] = immFitSum / ((double) nImmigrants);
            else
                avgImmFit[i] = -1.0;
            
            nResidents = n_in_each_patch[i] - migrationCount[i];
            if ( nResidents > 0 )
                avgResFit[i] = ((*fspt) - immFitSum) / ((double) nResidents);
            else
                avgResFit[i] = -1.0;
            
            if ( n_in_each_patch[i] > 0 )
                avgPatchFitness[i] = (*fspt) / ((double) n_in_each_patch[i]);
            else
                avgPatchFitness[i] = -1.0;
            
            
            fspt++;
        }
    }
    
    // how about calculating the maximum possible fitness in a patch?
    // then average as a proportion of the maximum?
    // how about average LD?
    
    fprintf(effMigRates,"%li", totalGenerationsElapsed);
    for ( j = 0; j < nPATCHES; j++ )
        fprintf(effMigRates," %G", eme[j]);
    fprintf(effMigRates," %i", nVariableLoci);
    for ( j = 0; j < nPATCHES; j++ ) {
        if ( PURE_NEUTRAL_MODEL ) {
            maxFit[j] = 1.0;
            minFit[j] = 1.0;
            randomFit[j] = 1.0;
        }
        else {
            // print average fitness of residents and recent immigrants
            maxFit[j] = calculateMaxPossibleFitness(j);
            minFit[j] = calculateMinPossibleFitness(j);
            if ( CONSIDER_EPISTASIS ) {  // these calculations are back of the envelope approximations
                if ( MULTIPLICATIVE_FITNESS )
                    randomFit[j] = sqrt( (maxFit[j] * minFit[j]) ); // geometric mean
                else
                    randomFit[j] = (maxFit[j] + minFit[j]) / 2.0; // arithmetic mean
            }
            else
                randomFit[j] = calcRandomHWFitness(j);
        }
        nImmigrants = migrationCount[j];
        nResidents = n_in_each_patch[j] - migrationCount[j];
        fprintf(effMigRates, " %i %i %E %E %E %E %E %E", nResidents, nImmigrants, avgResFit[j], avgImmFit[j], maxFit[j], minFit[j], randomFit[j], avgPatchFitness[j] );
    }
    j = nPATCHES - 1;
    if ( RECORD_FIT_TS ) {
        if ( n_in_each_patch[0] >= nRECORD_FIT )
            fprintf(fitTSdeme0, " %E %E %E %E %E %E\n", avgResFit[0], avgImmFit[0], maxFit[0], minFit[0], randomFit[0], avgPatchFitness[0]);
        if ( n_in_each_patch[j] >= nRECORD_FIT )
            fprintf(fitTSdeme1, " %E %E %E %E %E %E\n", avgResFit[j], avgImmFit[j], maxFit[j], minFit[j], randomFit[j], avgPatchFitness[j]);
    }
    
    if ( RECORD_LD_VALUES && (totalGenerationsElapsed > END_PERIOD_ALLOPATRY) ) {
        if ( BeginRecordingLD ) { // we already started recorded before now
            if ( (avgImmFit[0]/avgResFit[0] <= END_THRESH_FOR_LD) && migrationCount[0] && (avgImmFit[j]/avgResFit[j] <= END_THRESH_FOR_LD) && migrationCount[j] && !RECORDING_TIMES_IN_FILE ) {
                RECORD_LD_VALUES = 0;
                fprintf(stderr, "\nStopped recording LD at generation %li\n", totalGenerationsElapsed);
            }
        }
        else { // check to see if we should start recording
            if ( migrationCount[0] && migrationCount[j] ) {
                if ( (avgImmFit[0]/avgResFit[0] <= START_THRESH_FOR_LD) || (avgImmFit[j]/avgResFit[j] <= START_THRESH_FOR_LD) ) {
                    BeginRecordingLD = 1;
                    fprintf(stderr, "\nStarted recording LD at generation %li\n", totalGenerationsElapsed);
                }
            }
        }
    }
    
    if ( gatherLDvalues >= 1 ) {
        calculateLDselectedSitesOnly( gatherLDvalues );
        calculateLDneutralSitesOnly( gatherLDvalues );
        calculateLD( gatherLDvalues );
    }
    else
        fprintf(effMigRates, " 0.0 0.0 0.0 0.0 0.0 0.0\n");
    
    
    
    
    if ( (totalGenerationsElapsed > END_PERIOD_ALLOPATRY) && (m > 100) && !PURE_NEUTRAL_MODEL ) {
        if ( approachingSpeciationThreshold ) {
            // we already had some indication of *potentially*
            if ( (maxFit[0]/avgResFit[0] > avgResFit[0]/avgImmFit[0]) && (maxFit[j]/avgResFit[j] > avgResFit[j]/avgImmFit[j]) ) { // magnitude of resident fitness is still closer to that of immigrant than to max possible
                lastTimeResNearImm = totalGenerationsElapsed;
                // fprintf(stderr, "\ncandidate lastTimeResNearImm = %li\n", lastTimeResNearImm);
            }
            if ( (maxFit[0]/avgResFit[0] > avgResFit[0]/randomFit[0]) && (maxFit[j]/avgResFit[j] > avgResFit[j]/randomFit[j]) ) { // magnitude of resident fitness is still closer to that of immigrant than to max possible
                lastTimeResNearRand = totalGenerationsElapsed;
                // fprintf(stderr, "\ncandidate lastTimeResNearImm = %li\n", lastTimeResNearImm);
            }
        }
        else { // have not yet gotten indications of approaching threshold
            if ( (maxFit[0]/avgResFit[0] <= avgResFit[0]/avgImmFit[0]) || (maxFit[j]/avgResFit[j] <= avgResFit[j]/avgImmFit[j]) ) { // first time avg. res closer in magnitude to max than to avg. imm
                approachingSpeciationThreshold = 1;
                firstTimeResNearMaxFit = totalGenerationsElapsed;
                if ( RECORDING_TIMES_IN_FILE ) {
                    if ( recordingTimesCompleted > 0 ) {
                        lastTimeResNearImm = vectorOfRecordingTimes[(recordingTimesCompleted - 1)];
                        lastTimeResNearRand = lastTimeResNearImm;
                    }
                    else {
                        lastTimeResNearImm = 0;
                        lastTimeResNearRand = 0;
                    }
                }
                else {
                    lastTimeResNearImm = totalGenerationsElapsed - TS_SAMPLING_FREQUENCY;
                    lastTimeResNearRand = totalGenerationsElapsed - TS_SAMPLING_FREQUENCY;
                }
                fprintf(stderr, "\nfirstTimeResNearMaxFit = %li\n", firstTimeResNearMaxFit);
            }
        }
    }
    
    
    if ( totalGenerationsElapsed > END_PERIOD_ALLOPATRY && totalMigrants > 0  && !PURE_NEUTRAL_MODEL ) {
        
        stillLooking = 1;
        allPatchesDiverged = 0;
        j = 0;
        while ( (stillLooking) && (j < nPATCHES) ) {
            if ( (eme[j] * N) < RI_THRESH )
                allPatchesDiverged = 1;
            else {
                allPatchesDiverged = 0;
                stillLooking = 0;
            }
            j++;
        }
        
        if ( allPatchesDiverged ) {
            RI_REACHED = 1;
        }
    }
    
    
    
    //	if ( ((eme[(nPATCHES-1)] > (100000.0 * eme[0])) && (eme[0] > 0.0)) || ((eme[0] > (100000.0 * eme[(nPATCHES-1)])) && (eme[(nPATCHES-1)] > 0.0)) ) {
    //		int totalImmigrants, cct = 0;
    //		fprintf(stderr,"\n\n**** Warning: Strange values of effective migration rates ****\n\n");
    //		fprintf(stderr,"%li\teme[0] = %E\t\teme[1] = %E\n",totalGenerationsElapsed, eme[0], eme[1]);
    //		fprintf(stderr,"\tn1 = %i\t\tn2 = %i\n", n_in_each_patch[0], n_in_each_patch[1]);
    //		fvpt = fitpt;
    //		fspt = fitsumpt;
    //		pppt = previousPatches;
    //		for ( i = 0; i < nPATCHES; i++ ) {
    //			fprintf(stderr,"\npatch %i:\n",i);
    //			fprintf(stderr,"Prev patch\timmigrant?\tfitness\tcuml. fit.\txloc\tpatch\n");
    //			nj = n_in_each_patch[i];
    //			immFitSum = 0.0;
    //			totalImmigrants = 0;
    //			for ( j = 0; j < nj; j++ ) {
    //				fprintf(stderr,"%i\t\t",(*pppt));
    //				if ( j > 0 )
    //					focalFitVal = ((*fvpt) - (*(fvpt-1))); // because the fitnesses vector hold cumulative sums
    //				else
    //					focalFitVal = (*fvpt);
    //
    //				if ( (*pppt) != i ) { // this is a migrant
    //					fprintf(stderr,"1\t\t");
    //					totalImmigrants++;
    //					immFitSum += focalFitVal;
    //				}
    //				else
    //					fprintf(stderr,"0\t\t");
    //				fprintf(stderr,"%E\t%E\t%E\t%i\n", focalFitVal, (*fvpt), x_locations[cct], patch_locations[cct] );
    //				fvpt++;
    //				pppt++;
    //				cct++;
    //			}
    //			fprintf(stderr,"\ntotalImmigrants = %i, immFitSum = %E, fitSumForPatch = %E\n\n", totalImmigrants, immFitSum, (*fspt));
    //			fspt++;
    //		}
    //		fprintf(stderr,"\n\n\t**** Exiting program ****\n\n");
    //		testPrints();
    //		exit(1);
    //	}
    
    
}


double calcMeanMagnitude(double *valuesArray, long int n)
{
    long int i, j;
    double mu = 0.0, var = 0.0, *dpt;
    
    // get the mean
    dpt = valuesArray;
    for ( i = 0; i < n; i++ ) {
        mu += fabs( (*dpt) );
        dpt++;
    }
    return ( mu / ((double) n) );
    
}


double calcRandomHWFitness(int patchNum)
{
    int i, j;
    double fitVal = 1.0, p, q, f1, f0, fH;
    _Bool rl;
    long int fitArrayIndex;
    static int giveWarning = 1;
    
    if ( CONSIDER_EPISTASIS ) {
        if ( giveWarning ) {
            giveWarning = 0;
            fprintf(stderr,"\n\tWARNING: calcRandomHWFitness() does not account for\n\tCONSIDER_EPISTASIS\n\n");
        }
    }
    
    for ( i = 0; i < nLOCI; i++ ) {
        if ( IS_SELECTED_LOCUS[i] ) {
            rl = is_reversed_locus[i];
            fitArrayIndex = ( i * nPATCHES ) + patchNum;
            if ( fixed_allele[i] == 1 ) {
                if ( rl )
                    fitVal *= (*(fit0 + fitArrayIndex));
                else
                    fitVal *= (*(fit1 + fitArrayIndex));
            }
            else if ( fixed_allele[i] == 0 ) {
                fprintf(stderr, "\n Shouldn't be here in calcRandomHWFitness()\n");
                exit(1);
            }
            else {
                p = allele_frequencies[i];
                q = 1.0 - p;
                f1 = (*(fit1 + fitArrayIndex));
                f0 = (*(fit0 + fitArrayIndex));
                fH = (*(fitH + fitArrayIndex));
                if ( rl ) {
                    if ( MULTIPLICATIVE_FITNESS )
                        fitVal *= ( (p * p * f0) + (q * q * f1) + ( 2.0 * p * q * fH) );
                    else if ( ADDITIVE_FITNESS )
                        fitVal += ( (p * p * f0) + (q * q * f1) + ( 2.0 * p * q * fH) );
                }
                else {
                    if ( MULTIPLICATIVE_FITNESS )
                        fitVal *= ( (p * p * f1) + (q * q * f0) + ( 2.0 * p * q * fH) );
                    else if ( ADDITIVE_FITNESS )
                        fitVal += ( (p * p * f1) + (q * q * f0) + ( 2.0 * p * q * fH) );
                }
            }
        }
    }
    return fitVal;
}



double calculateMaxPossibleFitness(int patchNum)
{
    int i, j, l, p, locus, locus2, epiInteraction;
    long int fitArrayIndex[nLOCI], SELECTED_LOCI[nSELECTED_LOCI];
    _Bool locusTakenCareOf[nLOCI], derivedFitHighest, ancestralFitHighest;
    short int bestGT;
    double ecmult, ec, fitVal = 1.0, fv1, fv2, resid;
    
    ecmult = epi_patch_multipliers[patchNum];
    
    p = 0;
    for ( i = 0; i < nLOCI; i++ ) {
        fitArrayIndex[i] = 0;
        if ( IS_SELECTED_LOCUS[i] ) {
            SELECTED_LOCI[p] = i;
            p++;
        }
        locusTakenCareOf[i] = 0;
    }
    
    
    if ( ecmult < 0.0 ) {
        derivedFitHighest = 0; // this is the array, fit1, not the alleles per se
        ancestralFitHighest = 1; // in this case, the fit0 array will have higher values than fit1
    }
    else {
        derivedFitHighest = 1; // in this case, the fit1 array will have higher values than fit0
        ancestralFitHighest = 0;
    }
    
    ecmult = fabs(ecmult);
    
    for ( l = 0; l < nSELECTED_LOCI; l++ ) {
        locus = SELECTED_LOCI[l];
        fitArrayIndex[locus] = ( locus * nPATCHES ) + patchNum;
    }
#if ( MULTIPLICATIVE_FITNESS )
    if ( CONSIDER_EPISTASIS ) {
        for ( l = 0; l < nSELECTED_LOCI; l++ ) {
            locus = SELECTED_LOCI[l];
            
            if ( any_epis_for_this_locus[locus] ) {
                
                if ( l < (nSELECTED_LOCI - 1) ) {
                    for ( j = l+1; j < nSELECTED_LOCI; j++ ) { // ensures we don't count an interaction twice
                        locus2 = SELECTED_LOCI[j];
                        epiInteraction = *(epistasisMatrix + (locus * nLOCI) + locus2); // locus-th row, locus2-th column
                        
                        if ( epiInteraction == -1 ) {
                            if ( fixed_allele[locus] == 1 && fixed_allele[locus2] == 1 ) {
                                // negative DMI is forced
                                ec = *(epi_coeff_matrix + (nLOCI * locus) + locus2);
                                fitVal *= (1.0 - ec);  // reduction of 100*ec% in fitness of individual
                                locusTakenCareOf[locus] = 1;
                                locusTakenCareOf[locus2] = 1;
                            }
                        }
                        
                        else if ( epiInteraction == 1 ) { // positive interaction
                            if ( is_reversed_locus[locus] ) {
                                if ( derivedFitHighest )
                                    bestGT = 0;
                                else
                                    bestGT = 2;
                            }
                            else { // regular locus
                                if ( derivedFitHighest )
                                    bestGT = 2;
                                else
                                    bestGT = 0;
                            }
                            
                            if ( bestGT == 2 ) { // epis at this locus can work to increase max fitness
                                //  because derived is best in this patch; does not matter if 1 allele is fixed or not
                                ec = 1.0 + (ecmult * ( *(epi_coeff_matrix + (nLOCI * locus) + locus2) ));
                                // ec is now a multiplier that can be used below
                                
                                if ( is_reversed_locus[locus] ) {
                                    if ( !ancestralFitHighest ) {
                                        fprintf(stderr, "\nError!  fit0 should have been highest if you got to this point.\n");
                                        exit(1);
                                    }
                                    fv1 = (*(fit0 + fitArrayIndex[locus]));
                                    fv2 = (*(fit0 + fitArrayIndex[locus2]));
                                }
                                else {
                                    if ( !derivedFitHighest ) {
                                        fprintf(stderr, "\nError!  fit1 should have been highest if you got to this point.\n");
                                        exit(1);
                                    }
                                    fv1 = (*(fit1 + fitArrayIndex[locus]));
                                    fv2 = (*(fit1 + fitArrayIndex[locus2]));
                                }
                                
                                resid = (fv1 * fv2) - 1.0; // get the increase in fitness from this pair of loci
                                fitVal *= (1.0 + (ec * resid));
                                
                                locusTakenCareOf[locus] = 1;
                                locusTakenCareOf[locus2] = 1;
                            }
                            
                            else if ( fixed_allele[locus] == 1 && fixed_allele[locus2] == 1 ) {
                                // the best allele is 0 but the 1 alleles fixed, so GT == 0 or 1 isn't possible, and epi hurts
                                ec = 1.0 + (ecmult * ( *(epi_coeff_matrix + (nLOCI * locus) + locus2) ));
                                // ec is now a multiplier that can be used below
                                
                                if ( ancestralFitHighest ) {
                                    fv1 = (*(fit0 + fitArrayIndex[locus]));
                                    fv2 = (*(fit0 + fitArrayIndex[locus2]));
                                }
                                else {
                                    fv1 = (*(fit1 + fitArrayIndex[locus]));
                                    fv2 = (*(fit1 + fitArrayIndex[locus2]));
                                }
                                
                                resid = (fv1 * fv2) - 1.0; // get the increase in fitness from this pair of loci
                                fitVal *= (resid + 1.0) / (1.0 + (ec * resid)); // reduction in fitness
                                
                                locusTakenCareOf[locus] = 1;
                                locusTakenCareOf[locus2] = 1;
                            }
                        }
                    }
                }
            }
        }
    }
    // do the following, regardless of epistasis
    for ( l = 0; l < nSELECTED_LOCI; l++ ) {
        
        locus = SELECTED_LOCI[l];
        
        if ( !locusTakenCareOf[locus] ) { // effects of this locus on fitness have NOT yet been accounted for
            
            if ( fixed_allele[locus] == 1 ) {
                if ( is_reversed_locus[locus] )
                    fitVal *= (*(fit0 + fitArrayIndex[locus]));
                else
                    fitVal *= (*(fit1 + fitArrayIndex[locus]));
            }
            else {
                if ( derivedFitHighest ) // if it is a regular locus, then we use the fit1 array; if it were a reversed locus, then 0 would be the best allele, but we would reverse and use the fit1 array anyway
                    fitVal *= (*(fit1 + fitArrayIndex[locus]));
                else
                    fitVal *= (*(fit0 + fitArrayIndex[locus]));
            }
        }
    }
    
#elif ( ADDITIVE_FITNESS )
    if ( CONSIDER_EPISTASIS ) {
        fprintf(stderr, "\nHey dude, you didn't code all the calculations for additive fitness with epistasis yet.\n\tCode it and then compile again.\n\t... Exiting ...\n\n");
        exit(1);
    }
#endif
    
    
    if ( fitVal < 0.0 ) {
        fprintf(stderr, "\nBad max fitness: fitval = %E\n", fitVal);
        exit(1);
    }
    
    return fitVal;
}


double calculateMinPossibleFitness(int patchNum)
{
    int i, j, l, p, locus, locus2, epiInteraction;
    long int fitArrayIndex[nLOCI], SELECTED_LOCI[nSELECTED_LOCI];
    _Bool locusTakenCareOf[nLOCI], derivedFitHighest, ancestralFitHighest;
    short int bestGT;
    double ecmult, ec, fitVal = 1.0, fv1, fv2, resid;
    
    ecmult = epi_patch_multipliers[patchNum];
    
    p = 0;
    for ( i = 0; i < nLOCI; i++ ) {
        fitArrayIndex[i] = 0;
        if ( IS_SELECTED_LOCUS[i] ) {
            SELECTED_LOCI[p] = i;
            p++;
        }
        locusTakenCareOf[i] = 0;
    }
    
    
    if ( ecmult < 0.0 ) {
        derivedFitHighest = 0; // this is the array, not the alleles per se
        ancestralFitHighest = 1; // in this case, fit0 array will have higher values than fit1
    }
    else {
        derivedFitHighest = 1;
        ancestralFitHighest = 0;
    }
    
    ecmult = fabs(ecmult);
    
    for ( l = 0; l < nSELECTED_LOCI; l++ ) {
        locus = SELECTED_LOCI[l];
        fitArrayIndex[locus] = ( locus * nPATCHES ) + patchNum;
    }
#if ( MULTIPLICATIVE_FITNESS )
    if ( CONSIDER_EPISTASIS ) {
        for ( l = 0; l < nSELECTED_LOCI; l++ ) {
            locus = SELECTED_LOCI[l];
            
            if ( any_epis_for_this_locus[locus] ) {
                
                if ( l < (nSELECTED_LOCI - 1) ) {
                    for ( j = l+1; j < nSELECTED_LOCI; j++ ) { // ensures we don't count an interaction twice
                        locus2 = SELECTED_LOCI[j];
                        epiInteraction = *(epistasisMatrix + (locus * nLOCI) + locus2); // locus-th row, locus2-th column
                        
                        if ( epiInteraction == -1 ) {
                            // negative DMI
                            ec = *(epi_coeff_matrix + (nLOCI * locus) + locus2);
                            fitVal *= (1.0 - ec);  // reduction of 100*ec% in fitness of individual
                            locusTakenCareOf[locus] = 1;
                            locusTakenCareOf[locus2] = 1;
                        }
                        
                        else if ( epiInteraction == 1 ) { // positive interaction; only count if it is forced by fixation
                            if ( is_reversed_locus[locus] ) {
                                if ( derivedFitHighest )
                                    bestGT = 0;
                                else
                                    bestGT = 2;
                            }
                            else { // regular locus
                                if ( derivedFitHighest )
                                    bestGT = 2;
                                else
                                    bestGT = 0;
                            }
                            
                            if ( bestGT == 2 ) { // epis at this locus can work to increase max fitness
                                //  because derived is best in this patch; does not matter if 1 allele is fixed or not
                                if ( fixed_allele[locus] == 1 && fixed_allele[locus2] == 1 ) {
                                    ec = 1.0 + (ecmult * ( *(epi_coeff_matrix + (nLOCI * locus) + locus2) ));
                                    // ec is now a multiplier that can be used below
                                    
                                    if ( is_reversed_locus[locus] ) {
                                        if ( !ancestralFitHighest ) {
                                            fprintf(stderr, "\nError!  fit0 should have been highest if you got to this point.\n");
                                            exit(1);
                                        }
                                        fv1 = (*(fit0 + fitArrayIndex[locus]));
                                        fv2 = (*(fit0 + fitArrayIndex[locus2]));
                                    }
                                    else {
                                        if ( !derivedFitHighest ) {
                                            fprintf(stderr, "\nError!  fit1 should have been highest if you got to this point.\n");
                                            exit(1);
                                        }
                                        fv1 = (*(fit1 + fitArrayIndex[locus]));
                                        fv2 = (*(fit1 + fitArrayIndex[locus2]));
                                    }
                                    
                                    resid = (fv1 * fv2) - 1.0; // get the increase in fitness from this pair of loci
                                    fitVal *= (1.0 + (ec * resid));
                                    
                                    locusTakenCareOf[locus] = 1;
                                    locusTakenCareOf[locus2] = 1;
                                }
                            }
                            
                            else {
                                // the best allele is 0, so derived-derived pair can reduce fitness
                                ec = 1.0 + (ecmult * ( *(epi_coeff_matrix + (nLOCI * locus) + locus2) ));
                                // ec is now a multiplier that can be used below
                                
                                if ( ancestralFitHighest ) {
                                    fv1 = (*(fit0 + fitArrayIndex[locus]));
                                    fv2 = (*(fit0 + fitArrayIndex[locus2]));
                                }
                                else {
                                    fv1 = (*(fit1 + fitArrayIndex[locus]));
                                    fv2 = (*(fit1 + fitArrayIndex[locus2]));
                                }
                                
                                resid = (fv1 * fv2) - 1.0; // get the increase in fitness from this pair of loci
                                fitVal *= (resid + 1.0) / (1.0 + (ec * resid)); // reduction in fitness
                                
                                locusTakenCareOf[locus] = 1;
                                locusTakenCareOf[locus2] = 1;
                            }
                        }
                    }
                }
            }
        }
    }
    // do the following, regardless of epistasis
    for ( l = 0; l < nSELECTED_LOCI; l++ ) {
        
        locus = SELECTED_LOCI[l];
        
        if ( !locusTakenCareOf[locus] ) { // effects of this locus on fitness have NOT yet been accounted for
            
            if ( fixed_allele[locus] == 1 ) {
                if ( is_reversed_locus[locus] )
                    fitVal *= (*(fit0 + fitArrayIndex[locus]));
                else
                    fitVal *= (*(fit1 + fitArrayIndex[locus]));
            }
            else {
                if ( derivedFitHighest ) // if it is a regular locus, then we use the fit1 array; if it were a reversed locus, then 0 would be the best allele, but we would reverse and use the fit1 array anyway
                    fitVal *= (*(fit0 + fitArrayIndex[locus]));
                else
                    fitVal *= (*(fit1 + fitArrayIndex[locus]));
            }
        }
    }
    
#elif ( ADDITIVE_FITNESS )
    if ( CONSIDER_EPISTASIS ) {
        fprintf(stderr, "\nHey dude, you didn't code all the calculations for additive fitness with epistasis yet.\n\tCode it and then compile again.\n\t... Exiting ...\n\n");
        exit(1);
    }
#endif
    
    return fitVal;
}


double calcVariance(double *valuesArray, long int n)
{
    long int i, j;
    double mu = 0.0, var = 0.0, *dpt;
    
    // get the mean
    dpt = valuesArray;
    for ( i = 0; i < n; i++ ) {
        mu += (*dpt);
        dpt++;
    }
    mu = mu / ((double) n);
    
    dpt = valuesArray;
    for ( i = 0; i < n; i++ ) {
        var += ( ((*dpt) - mu) * ((*dpt) - mu) );
        dpt++;
    }
    
    return (var / ((double) n));
    
}

void chooseDemeSample(int *sampleArray, int demeNumber, int nSamples)
{
    int i, j, startIndex, endIndex;
    int focalN, dumi;
    _Bool chooseAgain;
    
    focalN = n_in_each_patch[demeNumber];
    
    startIndex = 0;
    if ( demeNumber > 0 ) {
        for ( i = 0; i < demeNumber; i++ ) {
            startIndex += n_in_each_patch[i];
        }
    }
    endIndex = (startIndex + focalN) - 1;
    
    dumi = startIndex + (randI() % focalN);
    if ( dumi < startIndex || dumi > endIndex ) {
        fprintf(stderr, "\nERROR!  Random sample algorithm busted!\n");
        fprintf(stderr, "startIndex = %i, endIndex = %i, dumi = %i\n", startIndex, endIndex, dumi);
        exit(-1);
    }
    sampleArray[0] = dumi;
    for ( i = 1; i < nSamples; i++ ) {
        chooseAgain = 1;
        while ( chooseAgain ) {
            chooseAgain = 0;
            dumi = startIndex + (randI() % focalN);
            for ( j = 0; j < i; j++ ) {
                if ( sampleArray[j] == dumi ) {
                    chooseAgain = 1;
                }
            }
        }
        sampleArray[i] = dumi;
        if ( (dumi < startIndex) || (dumi > endIndex) ) {
            fprintf(stderr, "\nERROR!  Random sample algorithm busted!\n");
            fprintf(stderr, "startIndex = %i, endIndex = %i, dumi = %i\n", startIndex, endIndex, dumi);
            exit(-1);
        }
    }
}


//void collectDetailedTimeSample(void)
//{
//    char dirname[80], counterString[10];
//
//    timeSampleCounter++;
//    sprintf(counterString, "%d", timeSampleCounter); // sample number as a character string
//    sprintf(dirname, "TimeSample%d", timeSampleCounter);
//    mkdir(dirname, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH );
//
//}



int compare_doubles(const void *a, const void *b)
{
    /* code used for qsort.  Code sourced from: http://www.gnu.org/software/libc/manual/html_node/Comparison-Functions.html#Comparison-Functions
     accessed on 4/26/12 */
    const double *da = (const double *) a;
    const double *db = (const double *) b;
    
    return (*da > *db) - (*da < *db);
}


int findNearestSelectedNeighbor(int focalLocus)
{
    int i, nsn = -1, cmemb, testLoc;
    double maplocation, minDist = 50.0, testDist;
    _Bool checkIt;
    
    maplocation = MAP[focalLocus];
    cmemb = chromosomeMembership[focalLocus];
    
    // look "down" first
    if (focalLocus > 0) {
        testLoc = focalLocus - 1;
        checkIt = 1;
        while ( (chromosomeMembership[testLoc] == cmemb) && checkIt ) {
            if ( IS_SELECTED_LOCUS[testLoc] && variable_loci[testLoc] ) {
                testDist = fabs( (maplocation - MAP[testLoc]) );
                if ( testDist < minDist ) {
                    minDist = testDist;
                    nsn = testLoc;
                }
                checkIt = 0; // first one has to be the closest below
            }
            testLoc--;
            if ( testLoc < 0 ) {
                checkIt = 0;
                testLoc = 0; // prevents bad indexing in while condition
            }
        }
    }
    
    // then look "up"
    if (focalLocus < (nLOCI-1)) {
        testLoc = focalLocus + 1;
        checkIt = 1;
        while ( (chromosomeMembership[testLoc] == cmemb) && checkIt ) {
            if ( IS_SELECTED_LOCUS[testLoc] && variable_loci[testLoc] ) {
                testDist = fabs( (MAP[testLoc] - maplocation) );
                if ( testDist < minDist ) {
                    minDist = testDist;
                    nsn = testLoc;
                }
                checkIt = 0; // first one found has to be the closest above
            }
            testLoc++;
            if ( testLoc >= nLOCI ) {
                checkIt = 0;
                testLoc = (nLOCI - 1); // safety to prevent bad indexing in while condition
            }
        }
    }
    
    // if nsn is still -1, that is a flag for none on this chromosome
    return nsn;
}



void growEpistasisMatrix(void)
{
    long int i, j, lasti;
    int newPositiveEpiCount = 0;
    int newDMIcount = 0;
    long int count = 0;
    long int nElements;
    short int *newEpistasisMatrix, *sipt, *oldsipt, *newsipt, dum1, dum2;
    double *newEpiCoeffMatrix, *dpt, *olddpt, *newdpt;
    
    nElements = nLOCI * nLOCI;
    
    newEpistasisMatrix = (short int *) palloc( POOL_EPISTASIS, (sizeof(short int) * nElements) );
    newEpiCoeffMatrix = (double *) palloc( POOL_EPI_COEFFS, (sizeof(double) * nElements) );
    
    if ( nSELECTED_LOCI < 2 ) {
        sipt = newEpistasisMatrix;
        dpt = newEpiCoeffMatrix;
        for ( i = 0; i < nElements; i++ ) {
            //printf("i = %li, *sipt = %i\n", i, *sipt);
            *sipt = 0;
            *dpt = 0.0;
            sipt++;
            dpt++;
            count++;
        }
    }
    else {
        // copy old relationships and assign new
        oldsipt = epistasisMatrix;
        newsipt = newEpistasisMatrix;
        olddpt = epi_coeff_matrix;
        newdpt = newEpiCoeffMatrix;
        for ( i = 0; i < nLOCI; i++ ) {
            if ( i == newestLocus ) {
                dum1 = is_reversed_locus[newestLocus];
                for ( j = 0; j < nLOCI; j++ ) {
                    dum2 = is_reversed_locus[j];
                    if ( (j == newestLocus) || (is_reference_locus[j]) ) { // can't have epistasis here
                        *newsipt = 0;
                        *newdpt = 0.0;
                    }
                    else if ( dum1 == dum2 ) {
                        if ( randU() < POSITIVE_EPI_PROBABILITY )  { // the two loci are eligible for a positive interaction
                            *newsipt = 1;
                            *newdpt = log(1.0 - randU()) * -MEAN_POSITIVE_EPI_COEFF;
                            any_epis_for_this_locus[i] = 1;
                            any_epis_for_this_locus[j] = 1;
                            newPositiveEpiCount++;
                            if ( dum1 )
                                nPositiveEpisReversed++;
                            else
                                nPositiveEpisRegular++;
                        }
                        else {
                            *newsipt = 0;
                            *newdpt = 0.0;
                        }
                    }
                    else if ( randU() < DMI_PROBABILITY ) { // the two loci are eligible for a DMI under our assumptions
                        *newsipt = -1;
                        *newdpt = log(1.0 - randU()) * -DMI_MEAN_EFFECT;
                        any_epis_for_this_locus[i] = 1;
                        any_epis_for_this_locus[j] = 1;
                        newDMIcount++;
                    }
                    else {
                        *newsipt = 0;
                        *newdpt = 0.0;
                    }
                    if ( *newsipt > 1 || *newsipt < -1 ) {
                        fprintf(stderr,"\n*newsipt = %i, i = %li, j = %li, nLOCI = %li, dum1 = %i, dum2 = %i, t = %li\n", *newsipt, i, j, nLOCI, dum1, dum2, totalGenerationsElapsed);
                        exit(1);
                    }
                    newsipt++;
                    newdpt++;
                    count++;
                }
            }
            else {
                for ( j = 0; j < nLOCI; j++ ) {
                    if ( j == newestLocus ) {
                        *newsipt = 0;
                        *newdpt = 0.0;
                        newsipt++; // skip this one for now; it is a spot for an interaction with the newest locus
                        newdpt++;
                    }
                    else { // this is a spot for copying from the old matrix
                        *newsipt = *oldsipt;
                        *newdpt = *olddpt;
                        newsipt++;
                        oldsipt++;
                        newdpt++;
                        olddpt++;
                    }
                    count++;
                }
            }
        }
        
    }
    
    if ( count != nElements ) {
        fprintf(stderr,"\nError!  Arithmetic of new epistasis interaction matrix is off:\n\tCount = %li, nElements = %li\n", count, nElements);
        exit(1);
    }
    
    
    if ( newPositiveEpiCount > 0 || newDMIcount > 0 ) {  /* at least one new epistatic interaction was added
                                                          need to adjust the newestLocus-th "row" of elements to match the corresponding column (symmetry)
                                                          above in the for loop in the part where the new counts are incremented, we have already set elements
                                                          (newestLocus * nLOCI) thru ((newestLocus+1)*nLOCI - 1) to the right values, just need to copy these
                                                          to make the matrix symmetrical */
        long int firstLoc, lastLoc;
        
        firstLoc = newestLocus;
        lastLoc = newestLocus + ((nLOCI - 1) * nLOCI);
        if ( lastLoc != (nElements - (nLOCI - newestLocus)) ) {
            fprintf(stderr,"\nError in arithmetic tidying up incompat. matrix:\n\tlastLoc = %li, but should have been %li\n", lastLoc, (nElements - (nLOCI - newestLocus)) );
            exit(1);
        }
        
        oldsipt = newEpistasisMatrix + (newestLocus * nLOCI); // we already filled in nLOCI new values starting here
        newsipt = newEpistasisMatrix + newestLocus;  // there are several, at nLOCI intervals, that we need to fill in starting here
        olddpt = newEpiCoeffMatrix + (newestLocus * nLOCI);
        newdpt = newEpiCoeffMatrix + newestLocus;
        for ( i = firstLoc; i <= lastLoc; i+=nLOCI ) { // starting at newestLocus, and jumping nLOCI elements each time
            *newsipt = *oldsipt;
            *newdpt = *olddpt;
            newsipt += nLOCI; // jump to next spot to fill
            newdpt += nLOCI;
            oldsipt++; // advance to next element, already filled above
            olddpt++;
            // at some point, these will point to the same address, but that should be fine
            // dbbp
            lasti = i;
        }
        
        if ( lasti != (nElements - (nLOCI - newestLocus)) ) {
            fprintf(stderr,"\nError in arithmetic tidying up incompat. matrix:\n\tlasti = %li, but should have been %li\n", lasti, (nElements - (nLOCI - newestLocus)) );
            exit(1);
        }
    }
    
    pfree( epistasisMatrix );
    epistasisMatrix = newEpistasisMatrix;
    pfree( epi_coeff_matrix );
    epi_coeff_matrix = newEpiCoeffMatrix;
    
    totalDMIs += newDMIcount;
    totalPositiveEpis += newPositiveEpiCount;
}




void initializePopulation(void)
{
    if ( INITIAL_CONDITIONS == 1 ) {
        // set up for secondary contact
    }
    else if ( INITIAL_CONDITIONS == 2 ) {
        loadCustomPopulationFromFiles();
    }
    else {
        if ( !FIXED_N )
            INITIAL_POPULATION_SIZE = ((int) K) * nPATCHES; // population at carrying capacity approximately in each patch; overrides command line option
        // default = random uniformly distributed population, monomorphic at all loci
        int i;
        for ( i = 0; i < nPATCHES; i++ )
            n_in_each_patch[i] = 0;
        // spatial locations
        patch_locations = (int*) palloc(POOL_PATCH_LOCATIONS, (sizeof(int) * INITIAL_POPULATION_SIZE) );
        //n_in_each_patch = (int*) palloc( (sizeof(int) * nPATCHES) );
        x_locations = (double*) palloc(POOL_X_LOCATIONS, (sizeof(double) * INITIAL_POPULATION_SIZE) );
#if (D == 2)
        y_locations = (double*) palloc(POOL_Y_LOCATIONS, (sizeof(double) * INITIAL_POPULATION_SIZE) );
#endif
        for ( i = 0; i < INITIAL_POPULATION_SIZE; i++ ) {
            // set all alleles to zero, locations to uniform random, assign patch numbers
            x_locations[i] = randU();
#if (D == 1)
            patch_locations[i] = (int) (x_locations[i] * ((double) PATCHES));
#elif (D == 2)
            y_locations[i] = randU();
            patch_locations[i] = scalarPatchNumber(x_locations[i],y_locations[i]);
#endif
            n_in_each_patch[(patch_locations[i])] = n_in_each_patch[(patch_locations[i])] + 1;
        }
        
        // patch numbers
        //		double incr = 1.0 / ((double) PATCHES);
        //		for ( i = 0; i < (PATCHES-1); i++) {
        //			patch_division_points[i] = (i+1) * incr;
        //		}
        //		patch_division_points[(PATCHES-1)] = 1.0;
        
        // intial genotypes
        genotypes = (short int*) palloc(POOL_GENOTYPES, (sizeof(short int) * INITIAL_POPULATION_SIZE * nLOCI * 2) );
        if ( EARLY_TEST_KILL > 0 )
            fprintf(stderr, "\nSize of genotypes in memory = %li\n", (sizeof(short int) * INITIAL_POPULATION_SIZE * nLOCI * 2));
        for ( i = 0; i < (INITIAL_POPULATION_SIZE * nLOCI * 2); i++ )
            genotypes[i] = 0;
        
        for ( i = 0; i < nLOCI; i++ )
            allele_frequencies[i] = 0.0;
        
        
        // sort the population; not a true sort, but rather a binning
        // no need to sort genotypes here bcause they are monomorphic.  Instead, just sort coordinates (x and y locations) and patch numbers (patch locations)
        
        double xsorted[INITIAL_POPULATION_SIZE], ysorted[INITIAL_POPULATION_SIZE];
        int patchSorted[INITIAL_POPULATION_SIZE], nextOpenSpot[nPATCHES], putMeHere, patch;
        nextOpenSpot[0] = 0;
        for ( i = 1; i < nPATCHES; i++ )
            nextOpenSpot[i] = n_in_each_patch[(i-1)] + nextOpenSpot[(i-1)];
        for ( i = 0; i < INITIAL_POPULATION_SIZE; i++ ) {
            patch = patch_locations[i];
            putMeHere = nextOpenSpot[patch];
            xsorted[putMeHere] = x_locations[i];
            patchSorted[putMeHere] = patch;
#if (D == 2)
            ysorted[putMeHere] = y_locations[i];
#endif
            nextOpenSpot[patch] = nextOpenSpot[patch] + 1;
        }
        // copy back to global
        for ( i = 0; i < INITIAL_POPULATION_SIZE; i++ ) {
            x_locations[i] = xsorted[i];
            patch_locations[i] = patchSorted[i];
#if (D == 2)
            y_locations[i] = ysorted[i];
#endif
        }
        N = INITIAL_POPULATION_SIZE;
    }
    
}


void loadCustomPopulationFromFiles(void)
{
    fprintf(stderr,"\nError! \n\tloadCustomPopulationFromFiles(): \n\tfunction called but not written yet!\n");
    exit(-1);
}

void loadMapFromFile(void)
{
    fprintf(stderr,"\nError! \n\tloadMapFromFile(): \n\tfunction called but not written yet!\n");
    exit(-1);
}


int maxInt(int i1, int i2)
{
    if ( i1 > i2 )
        return i1;
    else
        return i2;
}

int minInt(int i1, int i2)
{
    if ( i1 < i2 )
        return i1;
    else
        return i2;
}


void openDataFilesForRecording( int gatherLDvalues ) {
    int i, j;
    
    // open data files for recording
    effMigRates = fopen("EffectiveMigrationRates.txt","w");
    FSTtimeSeries = fopen("FSTtimeSeries.txt","w");
    AFtimeSeries = fopen("AlleleFreqTimeSeries.txt","w");
    grossMigRates = fopen("GrossMigrationRates.txt","w");
    fixationLog = fopen("FixationLog.txt","w");
    selectedFrequencies = fopen("SelectedAlleleFrequencies.txt","w");
    neutralFrequencies = fopen("NeutralAlleleFrequencies.txt","w");
    if ( RECORD_FIT_TS ) {
        fitTSdeme0 = fopen("FitnessTSdeme0.txt","w");
        fitTSdeme1 = fopen("FitnessTSdeme1.txt","w");
    }
    
    DXYTS = fopen("DXYtimeSeries.csv","w");
    fprintf(DXYTS, "Time");
    for ( i = 0; i < (nPATCHES - 1); i++ ) {
        for ( j = (i + 1); j < nPATCHES; j++ ) {
            fprintf(DXYTS, ",DxyDeme%iDeme%i", i, j);
        }
    }
    fprintf(DXYTS, "\n");
    
    
    nVarTS = fopen("NumberVariableLoci.txt","w");
    if ( CONSIDER_EPISTASIS )
        epistasisTS = fopen("EpistasisOverTime.txt","w");
    mutationLog = fopen("MutationLog.txt","w");
    fprintf(mutationLog,"time mut_num map_location chromosome chrom_loc S_MAX0 S_MAX1 H patch FitnessReversal MutationType\n");
    logOfRemovedLoci = fopen("LociRemoved.txt","w");
    fprintf(logOfRemovedLoci, "time locusID map_location chromosome S_MAX0 S_MAX1 H FitnessReversal MutationType\n");
    if ( gatherLDvalues >= 1 ) {
        LDselSitesAvg = fopen("LDselSitesAvg.txt","w");
        LDneutSitesAvg = fopen("LDneutSitesAvg.txt","w");
    }
    if ( gatherLDvalues >= 3 ) {
        LDselSitesSame = fopen("LDselSitesSame.txt","w");
        LDselSitesDiff = fopen("LDselSitesDiff.txt", "w");
        LDneutSitesSame = fopen("LDneutSitesSame.txt","w");
        LDneutSitesDiff = fopen("LDneutSitesDiff.txt", "w");
        LDfpt = fopen("LDvalues.txt","w");
        fprintf(LDfpt,"time locus1 locus2 LDcoef Dprime Delta MapDist S_MAX1loc1 S_MAX1loc2 locus1Type locus2type\n");
    }
    
    effPopSizeData = fopen("EffPopSizeData.txt","w");
    
    AFSTS = fopen("AlleleFrequencySpectrum.csv", "w");
    selAFSTS = fopen("SelectedAlleleFrequencySpectrum.csv", "w");
    neutAFSTS = fopen("NeutralAlleleFrequencySpectrum.csv", "w");
    fprintf( AFSTS, "Time,AlleleCount,NumSites\n" );
    fprintf( selAFSTS, "Time,AlleleCount,NumSites\n" );
    fprintf( neutAFSTS, "Time,AlleleCount,NumSites\n" );
    
    // for JAFS data:
    const char* datadir;
    datadir = "JAFSdata";
    struct stat sb;
    if (stat(datadir, &sb) == 0 && S_ISDIR(sb.st_mode)) {
        fprintf(stderr, "\n\tWARNING: JAFSdata directory already exists, so I'm using it ...\n");
        fprintf(stderr, "\tAs a result, data already in that directory may be overwritten, \n");
        fprintf(stderr, "\tdepending upon recording times and file names.\n");
    }
    else {
        printf("\n\tJAFSdata directory does not exist.\n\t\t... Creating it now ...\n");
        mkdir(datadir, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH );
    }
    
}



long int Poisson(double mm)
/* ==================================================
 * Returns a Poisson distributed non-negative integer.
 * NOTE: use mm > 0
 * code modified from http://www.cs.wm.edu/~va/software/park/rvgs.c
 * accessed 4/12/12
 * ==================================================
 */
{
    double tt = 0.0;
    long int events = 0;
    
    while (tt < mm) {  // cumulative "time" between events; mm is the expected events per unit of time
        tt += -(log(1.0 - randU())); /* adding an amount of "time" from a random exponential draw with mean 1.0 unit
                                      hence, this winds up treating mm as the expected number of units of time to achieve
                                      mm events.  the number of times we have to increment tt to reach mm is the actual number of
                                      events that occurred in mm units of time + 1. x counts the number of increments. s*/
        events++;  // one more "event"
    }
    return (events-1); // corrects for counting one past the last "event" in the while loop
}







void printParameters(long int estRunLength)
{
    int i, j;
    FILE *pfpt;
    
    pfpt = fopen("parameters.m","w");
    
    fprintf(pfpt,"CodeVersion = '%s';\n\n",version);
    
    fprintf(pfpt,"D = %i;\n",D);
    fprintf(pfpt,"PATCHES = %i;\n",PATCHES);
    fprintf(pfpt,"nPATCHES = %i;\n",nPATCHES);
    
    fprintf(pfpt,"USE_MUTATIONS_FROM_FILES = %i;\n",USE_MUTATIONS_FROM_FILES);
    fprintf(pfpt,"ADDITIVE_FITNESS = %i;\n",ADDITIVE_FITNESS);
    fprintf(pfpt,"MULTIPLICATIVE_FITNESS = %i;\n",MULTIPLICATIVE_FITNESS);
    fprintf(pfpt,"TWO_DEME = %i;\n",TWO_DEME);
    fprintf(pfpt,"DETERMINISTIC = %i;\n",DETERMINISTIC);
    fprintf(pfpt,"nGENERATIONS_MAX = %li;\n",nGENERATIONS_MAX);
    fprintf(pfpt,"nMUTATIONS = %i;\n",nMUTATIONS);
    fprintf(pfpt,"INITIAL_CONDITIONS = %i;\n",INITIAL_CONDITIONS);
    fprintf(pfpt,"INITIAL_POPULATION_SIZE = %i;\n",INITIAL_POPULATION_SIZE);
    fprintf(pfpt,"MEAN_S = %G;\n",MEAN_S);
    fprintf(pfpt,"SYMMETRIC_MUTATIONS = %i;\n",SYMMETRIC_MUTATIONS);
    fprintf(pfpt,"MAP_TYPE = %i;\n",MAP_TYPE);
    fprintf(pfpt,"GAMETE_PRODUCTION_MODE = %i;\n",GAMETE_PRODUCTION_MODE);
    fprintf(pfpt,"FIXED_N = %i;\n",FIXED_N);
    fprintf(pfpt,"K = %G;\n",K);
    fprintf(pfpt,"MOSAIC = %i;\n",MOSAIC);
    fprintf(pfpt,"SD_MOVE = %G;\n",SD_MOVE);
    fprintf(pfpt,"OFFSPRING_IN_RANDOM_LOCATIONS = %i;\n",OFFSPRING_IN_RANDOM_LOCATIONS);
    fprintf(pfpt,"nCHROMOSOMES = %i;\n",nCHROMOSOMES);
    fprintf(pfpt,"TOTAL_MAP_LENGTH = %G;\n",TOTAL_MAP_LENGTH);
    fprintf(pfpt,"MUTATION_DISTRIBUTION = %i;\n",MUTATION_DISTRIBUTION);
    fprintf(pfpt,"START_WITH_BIG_IN_PLACE = %i;\n",START_WITH_BIG_IN_PLACE);
    fprintf(pfpt,"BIG_START_S_VAL = %G;\n",BIG_START_S_VAL);
    fprintf(pfpt,"FATTEN_TAIL_PROPORTION = %G;\n",FATTEN_TAIL_PROPORTION);
    fprintf(pfpt,"FATTEN_TAIL_MAX = %G;\n",FATTEN_TAIL_MAX);
    fprintf(pfpt,"NUMBER_BIG_TO_ESTABLISH = %i;\n",NUMBER_BIG_TO_ESTABLISH);
    fprintf(pfpt,"END_PERIOD_ALLOPATRY = %li;\n",END_PERIOD_ALLOPATRY);
    fprintf(pfpt,"START_PERIOD_ALLOPATRY = %li;\n", START_PERIOD_ALLOPATRY);
    fprintf(pfpt,"PERIOD_ALLOPATRY = %i;\n", PERIOD_ALLOPATRY);
    fprintf(pfpt,"DMI_MEAN_EFFECT = %E;\n", DMI_MEAN_EFFECT);
    fprintf(pfpt,"DMI_PROBABILITY = %E;\n", DMI_PROBABILITY);
    fprintf(pfpt,"PROBABILITY_DERIVED_REVERSED = %E;\n", PROBABILITY_DERIVED_REVERSED);
    fprintf(pfpt,"POSITIVE_EPI_PROBABILITY = %E;\n", POSITIVE_EPI_PROBABILITY);
    fprintf(pfpt,"MEAN_POSITIVE_EPI_COEFF = %E;\n", MEAN_POSITIVE_EPI_COEFF);
    fprintf(pfpt,"START_THRESH_FOR_LD = %E;\n", START_THRESH_FOR_LD);
    fprintf(pfpt,"END_THRESH_FOR_LD = %E;\n", END_THRESH_FOR_LD);
    fprintf(pfpt,"MUTATIONS_PER_GENERATION = %i;\n", MUTATIONS_PER_GENERATION);
    fprintf(pfpt,"DEME0_CONSTANT_S = %E;\n", DEME0_CONSTANT_S);
    fprintf(pfpt,"DEME1_CONSTANT_S = %E;\n", DEME1_CONSTANT_S);
    fprintf(pfpt,"MU = %E;\n", ( 0.5 * ((double) N) / ((double) MUTATIONS_PER_GENERATION) ) );
    fprintf(pfpt,"FRACTION_SELECTED_LOCI = %E;\n", FRACTION_SELECTED_LOCI);
    fprintf(pfpt,"PURE_NEUTRAL_MODEL = %i;\n", PURE_NEUTRAL_MODEL);
    fprintf(pfpt,"PURE_ALLOPATRY_MODEL = %i;\n", PURE_ALLOPATRY_MODEL);
    
    
    fprintf(pfpt,"epi_patch_multipliers = [");
    for ( i = 0; i < nPATCHES; i++ )
        fprintf(pfpt,"%f ", epi_patch_multipliers[i]);
    fprintf(pfpt,"];\n");
    
    fprintf(pfpt,"H = [");
    for ( i = 0; i < nLOCI; i++ )
        fprintf(pfpt,"%E ",H[i]);
    fprintf(pfpt,"];\n");
    
    if ( MOSAIC  && nPATCHES > 2 ) {
        fprintf(pfpt,"BEST_ALLELE_IN_PATCH = [");
        for ( i = 0; i < nPATCHES; i++ )
            fprintf(pfpt,"%i ",BEST_ALLELE_IN_PATCH[i]);
        fprintf(pfpt,"];\n");
    }
    
    
    fprintf(pfpt,"\n%% Parameters and state variables that affected data recording or length of simulation (but not results themselves):\n");
    
    fprintf(pfpt,"TS_SAMPLING_FREQUENCY = %li;\n",TS_SAMPLING_FREQUENCY);
    fprintf(pfpt, "RECORDING_TIMES_IN_FILE = %i;\n", RECORDING_TIMES_IN_FILE);
    fprintf(pfpt,"RI_THRESH = %E;\n", RI_THRESH);
    fprintf(pfpt,"RI_REACHED = %i;\n\n", RI_REACHED);
    fprintf(pfpt,"RECORD_LD_VALUES = %i;\n", RECORD_LD_VALUES);
    fprintf(pfpt,"BeginRecordingLD = %i;\n", BeginRecordingLD);
    //    fprintf(pfpt,"LD_SAMPLING_FREQUENCY = %i;\n", LD_SAMPLING_FREQUENCY);
    fprintf(pfpt,"LD_LOCI_SUBSAMPLE = %i;\n", LD_LOCI_SUBSAMPLE);
    //    fprintf(pfpt,"EM_THRESH_FOR_LD = %G;\n", EM_THRESH_FOR_LD);
    fprintf(pfpt,"LD_LowerBound = %E;\n", LD_LowerBound);
    //    fprintf(pfpt,"LD_UpperBound = %G;\n\n", LD_UpperBound);
    fprintf(pfpt,"FST_MIN_RECORDING_THRESH = %E;\n", FST_MIN_RECORDING_THRESH);
    
    
    fprintf(pfpt,"%% Some useful variable states at the end of the simulation\n");
    fprintf(pfpt,"approachingSpeciationThreshold = %i;\n", approachingSpeciationThreshold);
    fprintf(pfpt,"firstTimeResNearMaxFit = %li;\n", firstTimeResNearMaxFit);
    fprintf(pfpt,"lastTimeResNearImm = %li;\n", lastTimeResNearImm);
    fprintf(pfpt,"lastTimeResNearRand = %li;\n", lastTimeResNearRand);
    fprintf(pfpt,"totalGenerationsElapsed = %li;\n",totalGenerationsElapsed);
    fprintf(pfpt,"estRunLength = %li;\n", estRunLength);
    fprintf(pfpt,"totalMutationsIntroduced = %li;\n",m);
    fprintf(pfpt,"N = %i;\n",N);
    fprintf(pfpt,"nLOCI = %li;\n",nLOCI);
    fprintf(pfpt,"nVariableLoci = %i;\n",nVariableLoci);
    fprintf(pfpt,"totalDMIs = %li;\n", totalDMIs);
    fprintf(pfpt,"totalPositiveEpis = %li;\n", totalPositiveEpis);
    fprintf(pfpt,"CONSIDER_EPISTASIS = %i;\n", CONSIDER_EPISTASIS);
    fprintf(pfpt, "LOCI_PER_CHROMOSOME = [");
    for ( i = 0; i < nCHROMOSOMES; i++ ) {
        fprintf(pfpt, "%i  ", LOCI_PER_CHROMOSOME[i]);
    }
    fprintf(pfpt, "];\n");
    
    
    
    fclose(pfpt);
}


int scalarPatchNumber(double x, double y)
{
    int xpatch, ypatch;
    
    if ( x >= 1.0 ) xpatch = PATCHES - 1;
    else xpatch = (int) ( x * ((double) PATCHES) ); // truncated, round-down multiplier
    
    if ( y >= 1.0 ) ypatch = PATCHES - 1;
    else ypatch = (int) ( y * ((double) PATCHES) );
    
    return ( ( xpatch * PATCHES ) + ypatch );
}





int seed_gen(void)
{
    /* use calendar time to seed random number generator.  Code adopted from Schildt's textbook */
    
    int stime;
    long ltime;
    FILE *rseed;
    
    /* get the calendar time */
    ltime=time(NULL);
    stime=(unsigned) ltime/2;
    
    
    // generate and store random number seed
    rseed = fopen("RnumSeed.txt","w");
    fprintf(rseed,"%i\n",stime);
    fclose(rseed);
    
    return stime;
}


long int setTSfreq(void)
{
    double tsf, maxGen;
    long int et;
    int i;
    long int incr;
    FILE *recTimesFile;
    
    // expected total generations to barrier of 200 = m / me from matlab multiple linear regression runs 200002 - 201700
    // log10(T) = 5.1411*SD_MOVE - 77.3223*MEAN_S - 3.6414e-06*INITIAL_POPULATION_SIZE - 0.1116*MUTATIONS_PER_GENERATION + 6.8657
    
    tsf = 5.1411*SD_MOVE - 77.3223*MEAN_S - 3.6414E-06*((double) INITIAL_POPULATION_SIZE) - 0.1116*((double) MUTATIONS_PER_GENERATION) + 6.8657;
    
    tsf = pow(10.0, tsf);
    et = (long int) tsf;
    
    tsf = tsf / ((double) nTIME_SAMPLES);
    
    TS_SAMPLING_FREQUENCY = (long int) tsf;
    
    if ( TS_SAMPLING_FREQUENCY < MIN_TS_SAMP_FREQ )
        TS_SAMPLING_FREQUENCY = MIN_TS_SAMP_FREQ;
    
    if ( TS_SAMPLING_FREQUENCY > MAX_TS_SAMP_FREQ )
        TS_SAMPLING_FREQUENCY = MAX_TS_SAMP_FREQ;
    
    if ( et < 1 )
        et = 1;
    
    if ( RECORDING_TIMES_IN_FILE ) {
        
        recTimesFile = fopen("RecordingTimes.txt", "r");
        if ( recTimesFile == NULL ) {
            //if ( PURE_NEUTRAL_MODEL || PURE_ALLOPATRY_MODEL ) {
            // make the file for use and for metadata purposes:
            recTimesFile = fopen("RecordingTimes.txt", "w");
            nRecordingTimes = 5;
            fprintf(recTimesFile, "%i\n", nRecordingTimes); // number of recording times
            maxGen = ((double) nMUTATIONS) / ((double) MUTATIONS_PER_GENERATION);
            incr = (long int) ( maxGen / ((double) nRecordingTimes) );
            for ( i = 1; i < nRecordingTimes; i++ )
                fprintf(recTimesFile, "%li ", (incr * i));
            fprintf(recTimesFile, "%li", ((long int) maxGen) - 1);
            
            fclose(recTimesFile);
            recTimesFile = fopen("RecordingTimes.txt", "r");
            //}
            //            else {
            //                fprintf(stderr, "\nError in setTSfreq():\n\tRecordingTimes.txt not found!\n");
            //                exit(-1);
            //            }
        }
        fscanf(recTimesFile, "%i", &nRecordingTimes);
        vectorOfRecordingTimes = (long int *) malloc( nRecordingTimes * sizeof(long int) );
        for ( i = 0; i < nRecordingTimes; i++ ) {
            if ( feof(recTimesFile) ) {
                printf("\n\tError in setTSfreq():\n\t\tWrong number of rec times in file!\n");
                exit(-1);
            }
            else
                fscanf( recTimesFile, "%li", (vectorOfRecordingTimes+i) );
            
            if ( vectorOfRecordingTimes[i] < 0 ) {
                printf("\n\tError in setTSfreq():\n\t\tNegative time in rec times in file!\n");
                exit(-1);
            }
            if ( i > 0 ) {
                if ( vectorOfRecordingTimes[(i - 1)] >= vectorOfRecordingTimes[i] ) {
                    printf("\n\tError in setTSfreq():\n\t\tRecording times in file are not strictly increasing!\n");
                    exit(-1);
                }
            }
        }
        
        fclose(recTimesFile);
        nextRecordingTime = vectorOfRecordingTimes[0];
        
        printf("\n\tData will be recorded at the following %i times:\n\t\t", nRecordingTimes);
        for ( i = 0; i < nRecordingTimes; i++ ) {
            printf("%li\t", vectorOfRecordingTimes[i]);
        }
        printf("\n\n");
        
        // make sure no other times are recorded:
        TS_SAMPLING_FREQUENCY = 100 * nMUTATIONS;
        
        BeginRecordingLD = 1;
        //exit(0);
    }
    
    return et;
}



void shrinkEpistasisMatrix(int locusToRemove)
{
    long int i, j;
    short int *oldsipt, *newsipt;
    int nDMIsLost = 0, nPositiveEpisLost = 0;
    double *olddpt, *newdpt;
    
    oldsipt = epistasisMatrix;
    newsipt = epistasisMatrix;
    olddpt = epi_coeff_matrix;
    newdpt = epi_coeff_matrix;
    
    //nLOCI has NOT yet been decremented
    
    for ( i = 0; i < nLOCI; i++ ) {
        if ( i == locusToRemove ) {
            // skip over all these entries
            for ( j = 0; j < nLOCI; j++ ) {
                if ( (*oldsipt) == -1 )
                    nDMIsLost++;
                else if ( (*oldsipt) == 1 ) {
                    nPositiveEpisLost++;
                    if ( is_reversed_locus[i] != is_reversed_locus[j] ) {
                        fprintf(stderr, "algorithm busted.  rev i = %i, rev j = %i\n", is_reversed_locus[i], is_reversed_locus[j]);
                    }
                    if ( is_reversed_locus[j] )
                        nPositiveEpisReversed--;
                    else {
                        nPositiveEpisRegular--;
                        if ( nPositiveEpisRegular < 0 ) {
                            fprintf(stderr,"t = %li, locus lost = %i, j = %li, posreg = %li < 0, rev = %i\n", totalGenerationsElapsed, locusToRemove, j, nPositiveEpisRegular, is_reversed_locus[i]);
                        }
                    }
                }
                oldsipt++;
                olddpt++;
            }
        }
        else {
            for ( j = 0; j < nLOCI; j++ ) {
                if ( j == locusToRemove ) {
                    oldsipt++; // skip it
                    olddpt++;
                }
                else {
                    *newsipt = *oldsipt;
                    *newdpt = *olddpt;
                    newsipt++;
                    oldsipt++;
                    newdpt++;
                    olddpt++;
                }
                
            }
        }
    }
    
    totalDMIs -= nDMIsLost;
    totalPositiveEpis -= nPositiveEpisLost;
    
}







void sortPopulation(void)
{
    int i, j, p, spot, *ript;
    short int *sortedGenotypes = (short int*) palloc(POOL_GENOTYPES, sizeof(short int) * N * 2 * nLOCI );
    short int *opt, *npt, *rspt;
    double *sorted_x = (double*) palloc(POOL_X_LOCATIONS, sizeof(double) * N );
    double *rdpt;
    int *sortedPatchLocations = (int*) palloc(POOL_PATCH_LOCATIONS, sizeof(int) * N );
    int *sortedPrevPatLoc = (int*) palloc(POOL_PATCH_LOCATIONS, sizeof(int) * N );
    int putMeHere[nPATCHES];
#if ( D == 2 )
    double *sorted_y = (double*) palloc(POOL_Y_LOCATIONS, sizeof(double) * N );
#endif
    
    int lastp = -1;
    int samePatchCnt = 0;
    
    printf("sortPopulation\n");
    putMeHere[0] = 0;
    for ( i = 1; i < nPATCHES; i++ )
        putMeHere[i] = putMeHere[(i-1)] + n_in_each_patch[(i-1)];
    
    for ( i = 0; i < N; i++ ) {
        p = patch_locations[i];
        //if (p != 1)
        //printf("p = %d ", p);
        
        if (p == lastp)
            samePatchCnt++;
        else {
            printf("%d\n", samePatchCnt);
            samePatchCnt = 0;
            printf("%d: ", p);
            lastp = p;
        }
        
        spot = putMeHere[p];
        
        sortedPatchLocations[spot] = p;
        sortedPrevPatLoc[spot] = previousPatches[i];
        
        sorted_x[spot] = x_locations[i];
#if ( D == 2 )
        sorted_y[spot] = y_locations[i];
#endif
        putMeHere[p] = putMeHere[p] + 1;
        
        opt = genotypes + ( i * 2 * nLOCI );
        npt = sortedGenotypes + ( spot * 2 * nLOCI );
        memcpy(npt, opt, sizeof(short int) * 2 * nLOCI);
    }
    
    pfree(x_locations);
    x_locations = sorted_x;
    pfree(previousPatches);
    previousPatches = sortedPrevPatLoc;
    
#if ( D == 2)
    pfree(y_locations);
    y_locations = sorted_y;
#endif
    
    pfree(genotypes);
    genotypes = sortedGenotypes;
    
    pfree(patch_locations);
    patch_locations = sortedPatchLocations;
    printf("%d\n", samePatchCnt);
    //printf("same patch %d / %d times\n", samePatchCnt, N);
    
    lastp = -1;
    samePatchCnt = 0;
    
    printf("\n");
    for ( i = 0; i < N; i++ ) {
        p = patch_locations[i];
        if (p == lastp)
            samePatchCnt++;
        else {
            printf("%d\n", samePatchCnt);
            samePatchCnt = 0;
            printf("%d: ", p);
            lastp = p;
        }
    }
    printf("%d\n", samePatchCnt);
    
}


void sortPopulation2(void)
{
    int i, j, p, spot, *ript;
    short int *sortedGenotypes = (short int*) palloc(POOL_GENOTYPES, sizeof(short int) * N * 2 * nLOCI );
    short int *opt, *npt, *rspt;
    double *sorted_x = (double*) palloc(POOL_X_LOCATIONS, sizeof(double) * N );
    double *rdpt;
    int *sortedPatchLocations = (int*) palloc(POOL_PATCH_LOCATIONS, sizeof(int) * N );
    int *sortedPrevPatLoc = (int*) palloc(POOL_PATCH_LOCATIONS, sizeof(int) * N );
    int putMeHere[nPATCHES];
    int adjacentPatchCount, lastPatch;
#if ( D == 2 )
    double *sorted_y = (double*) palloc(POOL_Y_LOCATIONS, sizeof(double) * N );
#endif
    
    putMeHere[0] = 0;
    for ( i = 1; i < nPATCHES; i++ )
        putMeHere[i] = putMeHere[(i-1)] + n_in_each_patch[(i-1)];
    
    
    //
    // Set up initial conditions for batch copying of loci
    //
    lastPatch = patch_locations[0];
    adjacentPatchCount = 0;
    opt = genotypes;
    p = patch_locations[0];
    if ( p > 0 )
        npt = sortedGenotypes + ( 2 * nLOCI * (putMeHere[p]) );
    else
        npt = sortedGenotypes; // when p == 0
    
    //
    // Go through all members, placing them in order within their patches
    //
    for ( i = 0; i < N; i++ ) {
        p = patch_locations[i];
        
        spot = putMeHere[p];
        sortedPatchLocations[spot] = p;
        sorted_x[spot] = x_locations[i];
        sortedPrevPatLoc[spot] = previousPatches[i];
#if ( D == 2 )
        sorted_y[spot] = y_locations[i];
#endif
        putMeHere[p] = putMeHere[p] + 1;
        
        if (p == lastPatch)
            adjacentPatchCount++;
        else {
            //
            // Copy all loci that were adjacent
            //
            memcpy(npt, opt, sizeof(short int) * 2 * nLOCI * adjacentPatchCount);
            
            //
            // Reset all parameters for next batch of loci to copy
            //
            lastPatch = p;
            adjacentPatchCount = 1;	// always have 1 to copy after initial condition
            opt = genotypes + ( i * 2 * nLOCI );
            npt = sortedGenotypes + ( spot * 2 * nLOCI );
        }
    }
    
    if (adjacentPatchCount)	// Copy any remaining loci
        memcpy(npt, opt, sizeof(short int) * 2 * nLOCI * adjacentPatchCount);
    
    pfree(x_locations);
    x_locations = sorted_x;
    
    pfree(previousPatches);
    previousPatches = sortedPrevPatLoc;
    
#if ( D == 2)
    pfree(y_locations);
    y_locations = sorted_y;
#endif
    
    pfree(genotypes);
    genotypes = sortedGenotypes;
    
    pfree(patch_locations);
    patch_locations = sortedPatchLocations;
}


void testPrints(void)
{
    if ( EARLY_TEST_KILL > 0 ) {
        fprintf(stderr, "\n");
        
        fprintf(stderr, "Version = %s\n\n",version);
    }
    //fprintf(stderr, "nGENERATIONS_MAX_DEFAULT = %li\n", nGENERATIONS_MAX);
    
    int i,j,count;
    FILE *fpt;
    fpt = fopen("FinalLocations.txt","w");
    for ( i = 0; i < N; i++ ) {
#if (D == 1)
        fprintf(fpt,"%G,  %i\n", x_locations[i], patch_locations[i]);
#elif (D == 2)
        fprintf(fpt,"%G,  %G,  %i\n", x_locations[i], y_locations[i], patch_locations[i]);
#endif
    }
    fclose(fpt);
    
    
    //	fprintf(stderr, "Mutation order:\n");
    //	for ( i = 0; i < nMUTATIONS; i++ ) {
    //		fprintf(stderr, "%i, ",MUTATION_ORDER[i]);
    //		if ( (i-9) % 10 == 0 ) fprintf(stderr, "\n");
    //	}
    //
    //	fprintf(stderr, "\n\nS_MAX1\t S_MAX0:\n");
    //	for ( i = 0; i < nLOCI; i++ ) {
    //		fprintf(stderr, "%G\t%G\n", S_MAX1[i],S_MAX0[i]);
    //	}
    
    /* // test
     fpt = fopen("TestExp.txt","w");
     for ( i = 0; i < 100000; i++ ) fprintf(fpt,"%G\n",( log(1 - randU()) * -MEAN_S ) );
     fclose(fpt); // */
    
    
    if ( EARLY_TEST_KILL > 1 ) {
        fprintf(stderr, "\n\n");
        fprintf(stderr, "\n\nMap:\n");
        count = 0;
        for (i = 0; i < nCHROMOSOMES; i++) {
            for ( j = 0; j < LOCI_PER_CHROMOSOME[i]; j++ ) {
                fprintf(stderr, "%G, ", MAP[count]);
                count++;
            }
            fprintf(stderr, "\n");
        }
    }
    if ( EARLY_TEST_KILL > 0 ) {
        fprintf(stderr, "\n\nNumber of chromosomes = %i\n\n", nCHROMOSOMES);
        
        fprintf(stderr, "\nLOCI_PER_CHROMOSOME:\n");
        for (i=0; i<nCHROMOSOMES; i++) {
            fprintf(stderr, "%i  ", LOCI_PER_CHROMOSOME[i]);
        }
        fprintf(stderr, "\n\n");
        
        fprintf(stderr, "\nMAP_LENGTHS:\n");
        for ( i = 0; i < nCHROMOSOMES; i++ )
            fprintf(stderr, "%G  ",MAP_LENGTHS[i]);
        fprintf(stderr, "\n\n");
    }
    
    /*
     fpt = fopen("FinalGenotypes.txt","w");
     count = 0;
     for ( i = 0; i < N; i++ ) {
     for ( j = 0; j < (2*nLOCI); j++ ) {
     fprintf(fpt,"%i ",genotypes[count++]);
     }
     fprintf(fpt,"\n");
     }
     fclose(fpt);
     */
    
    fpt = fopen("AllelesEstablished.m","w");
    fprintf(fpt,"variable_loci = [");
    for ( i = 0; i < nLOCI; i++ )
        fprintf(fpt,"%i ", variable_loci[i]);
    fprintf(fpt,"];\n");
    
    fprintf(fpt,"fixed_allele = [");
    for ( i = 0; i < nLOCI; i++ )
        fprintf(fpt,"%i ", fixed_allele[i]);
    fprintf(fpt,"];\n");
    fclose(fpt);
    
    
    FILE *pfpt;
    
    if ( MUTATION_DISTRIBUTION == 2 ) { // could have changed first few values
        fpt = fopen("Mutations.m","w");  // overwrite original file set up at initialization
        fprintf(fpt,"MUTATION_ORDER = [");
        for ( i = 0; i < nMUTATIONS; i++ )
            fprintf(fpt, "%E ", MUTATION_ORDER[i]);
        fprintf(fpt,"];\n");
        fprintf(fpt,"MUTATION_TYPE_SEQUENCE = [");
        for ( i = 0; i < nMUTATIONS; i++ )
            fprintf(fpt, "%i ", MUTATION_TYPE_SEQUENCE[i]);
        fprintf(fpt,"];\n");
        fprintf(fpt,"MUTATION_LOCATIONS = [");
        for ( i = 0; i < nMUTATIONS; i++ )
            fprintf(fpt, "%i ", MUTATION_LOCATIONS[i]);
        fprintf(fpt,"];\n");
        fprintf(fpt,"S_MAX0_SEQUENCE = [");
        for ( i = 0; i < nMUTATIONS; i++ )
            fprintf(fpt, "%E ", S_MAX0_SEQUENCE[i]);
        fprintf(fpt,"];\n");
        fprintf(fpt,"S_MAX1_SEQUENCE = [");
        for ( i = 0; i < nMUTATIONS; i++ )
            fprintf(fpt, "%E ", S_MAX1_SEQUENCE[i]);
        fprintf(fpt,"];\n");
        fclose(fpt);
    }
    
    fpt = fopen("lociSvalues.m","w");  // overwrite
    fprintf(fpt,"S_MAX0 = [");
    for ( i = 0; i < nLOCI; i++ )
        fprintf(fpt, "%E ", S_MAX0[i]);
    fprintf(fpt,"];\n");
    
    fprintf(fpt,"S_MAX1 = [");
    for ( i = 0; i < nLOCI; i++ )
        fprintf(fpt, "%E ", S_MAX1[i]);
    fprintf(fpt,"];\n");
    fclose(fpt);
    
    pfpt = fopen("FitnessMatrices.m","w"); // overwrite
    
    fprintf(pfpt, "fit0 = [");
    for ( i = 0; i < nLOCI; i++ ) {
        for ( j = 0; j < nPATCHES; j++ )
            fprintf(pfpt,"%E ",fit0[((i * nPATCHES) + j)]);
        if ( i < (nLOCI-1) )
            fprintf(pfpt,";\n");
    }
    fprintf(pfpt,"];\n\n");
    
    fprintf(pfpt, "fit1 = [");
    for ( i = 0; i < nLOCI; i++ ) {
        for ( j = 0; j < nPATCHES; j++ )
            fprintf(pfpt,"%E ",fit1[((i * nPATCHES) + j)]);
        if ( i < (nLOCI-1) )
            fprintf(pfpt,";\n");
    }
    fprintf(pfpt,"];\n\n");
    
    fprintf(pfpt, "fitH = [");
    for ( i = 0; i < nLOCI; i++ ) {
        for ( j = 0; j < nPATCHES; j++ )
            fprintf(pfpt,"%E ",fitH[((i * nPATCHES) + j)]);
        if ( i < (nLOCI-1) )
            fprintf(pfpt,";\n");
    }
    fprintf(pfpt,"];\n");
    
    fclose(pfpt);
    
    
    FILE *mpt;
    mpt = fopen("Map.m","w");
    fprintf(mpt,"MAP = [");
    for ( i = 0; i < nLOCI; i++ )
        fprintf(mpt,"%.10E ",MAP[i]);
    fprintf(mpt,"];\n\n");
    
    int locus = 0;
    for ( i = 0; i < nCHROMOSOMES; i++ ) {
        fprintf(mpt,"Chromosome%i = [",i);
        for ( j = 0; j < LOCI_PER_CHROMOSOME[i]; j++ ) {
            fprintf(mpt,"%.10E ", MAP[locus]);
            locus++;
        }
        fprintf(mpt,"];\n");
    }
    fclose(mpt);
    if (locus != nLOCI ) {
        fprintf(stderr, "\n\nERROR!! locus = %i, nLOCI = %li \n\n", locus, nLOCI);
        exit(1);
    }
    if ( EARLY_TEST_KILL > 1 ) {
        for ( i = 1; i < nLOCI; i++ ) {
            if ( MAP[i] < MAP[(i-1)] ) {
                fprintf(stderr, "\nChromosome break at: MAP[%i] = %E, MAP[%i] = %E",i,MAP[i],(i-1),MAP[(i-1)]);
            }
            
        }
        fprintf(stderr, "\n");
    }
    
}


void
usage(char *s)
{
    fprintf(stderr,  "Usage: %s [options]\n\n", s);
    
    fprintf(stderr,  "\tNOTE: There is no error checking on user inputs.\n");
    fprintf(stderr,  "\tProgram behavior is unpredictable (not defined) for\n");
    fprintf(stderr,  "\tinputs that do not conform to guidelines below.\n\n");
    
    fprintf(stderr,  "\t[-A <integer>] \tLength of period of ALLOPATRY, expressed in\n");
    fprintf(stderr,  "\t\t\tterms of numbers of mutations introduced (not necessarily\n");
    fprintf(stderr,  "\t\t\taccumulated). Set to 0 for primary contact, i.e., divergence\n");
    fprintf(stderr,  "\t\t\twith constant gene flow, i.e., sympatry/parapatry.\n");
    fprintf(stderr,  "\t\t\tDefault is %i mutations\n\n", END_PERIOD_ALLOPATRY_DEFAULT);
    
    fprintf(stderr,  "\t[-B <integer>] \tNumber of large-effect mutations to establish\n");
    fprintf(stderr,  "\t\t\tearly in a run.  Used with '-f 2' or '-b 1' (see below).\n");
    fprintf(stderr,  "\t\t\tDefault is %i\n\n",NUMBER_BIG_TO_ESTABLISH_DEFAULT);
    
    fprintf(stderr,  "\t[-b <0 or 1>] \tStart with -B large effect mutations established\n");
    fprintf(stderr,  "\t\t\tat approximate migration-selection balance.\n");
    fprintf(stderr,  "\t\t\tDefault is %i\n\n",START_WITH_BIG_IN_PLACE_DEFAULT);
    
    fprintf(stderr,  "\t[-C <integer>] \tNumber of chromosomes.\n");
    fprintf(stderr,  "\t\t\tDefault is %i\n\n",nCHROMOSOMES_DEFAULT);
    
    fprintf(stderr,  "\t[-D]\t\tRun in deterministic mode.\n");
    fprintf(stderr,  "\t\t\tWith -D, the random number seed\n");
    fprintf(stderr,  "\t\t\tis read from the file 'RnumSeed.txt'.\n");
    fprintf(stderr,  "\t\t\tWith -D, if 'RnumSeed.txt' is not\n");
    fprintf(stderr,  "\t\t\tin the current wd, program will seg fault.\n");
    fprintf(stderr,  "\t\t\tWithout -D, program will seed RNG with system time\n\n");
    
    fprintf(stderr,  "\t[-d <sigma>] \tStandard deviation of dispersal distances.\n");
    fprintf(stderr,  "\t\t\tDefault is %G\n",SD_MOVE_DEFAULT);
    fprintf(stderr,  "\t\t\tRecall that the habitat is the unit line\n");
    fprintf(stderr,  "\t\t\tor unit square, so <sigma> is also the fraction\n");
    fprintf(stderr,  "\t\t\tof the habitat's width covered by a movement\n");
    fprintf(stderr,  "\t\t\tof length <sigma> in one dimension.\n\n");
    
    fprintf(stderr,  "\t[-e <0 or 1>]\tRun in two deme mode with -e 1.\n");
    fprintf(stderr,  "\t\t\tDefault is %i\n",TWO_DEME_DEFAULT);
    fprintf(stderr,  "\t\t\t-e 1 option should ONLY be called if the \n");
    fprintf(stderr,  "\t\t\tdimensionality of the habitat is 1 and \n");
    fprintf(stderr,  "\t\t\tthe total number of patches is 2.\n");
    fprintf(stderr,  "\t\t\tThis option eliminates effects of\n");
    fprintf(stderr,  "\t\t\tcontinuous space and essentially creates\n");
    fprintf(stderr,  "\t\t\ta two-patch island model.\n\n");
    
    
    fprintf(stderr,  "\t[-F <0,1>] \tFixed population size\n");
    fprintf(stderr,  "\t\t\t'-F 1': population size is fixed\n");
    fprintf(stderr,  "\t\t\t'-F 0': logistic (density-dependent) growth\n");
    fprintf(stderr,  "\t\t\tDefault is %i; see also '-K' below.\n\n",FIXED_N_DEFAULT);
    
    fprintf(stderr,  "\t[-f <0,1,2>] \tDistribution of mutation S values\n");
    fprintf(stderr,  "\t\t\t'-f 0': exponential distribution with mean set by '-S'\n");
    fprintf(stderr,  "\t\t\t'-f 1': 'fattened' tail (see also usage of '-P' below)\n");
    fprintf(stderr,  "\t\t\t'-f 2': 'flattened' distribution with early mutations\n");
    fprintf(stderr,  "\t\t\t\tbeing drawn from uniform distribution on [0,%G]\n",FATTEN_TAIL_MAX);
    fprintf(stderr,  "\t\t\t\tand later mutations following exponential\n\t\t\t\t(see '-B' above).\n");
    fprintf(stderr,  "\t\t\tDefault is %i\n\n",MUTATION_DISTRIBUTION_DEFAULT);
    
    fprintf(stderr,  "\t[-G <0,1,2>] \tGamete production mode.\n");
    fprintf(stderr,  "\t\t\t'-G 0': divergence and genome hitchhiking\n");
    fprintf(stderr,  "\t\t\t'-G 1': genome hitchhiking\n");
    fprintf(stderr,  "\t\t\t'-G 2': direct selection only (loci inherited\n\t\t\t\tindependently)\n");
    fprintf(stderr,  "\t\t\tDefault is %i\n\n",GAMETE_PRODUCTION_MODE_DEFAULT);
    
    fprintf(stderr,  "\t[-H <value>] \tDominance coefficient of '1' alleles over '0' alleles.\n");
    fprintf(stderr,  "\t\t\tDefault is %G; must lie between 0 and 1.\n\n",H_DEFAULT);
    
    fprintf(stderr,  "\t[-I <value>] \tProportion of pairs of loci with an incompatibility.\n");
    fprintf(stderr,  "\t\t\tDefault is %f.\n\n",DMI_PROB_DEFAULT);
    
    fprintf(stderr,  "\t[-i <value>] \tProportion of pairs of loci positive epistasis.\n");
    fprintf(stderr,  "\t\t\tDefault is %f.\n\n",POSITIVE_EPI_PROB_DEFAULT);
    
    fprintf(stderr,  "\t[-K <value>] \tCarrying capacity of each patch.\n");
    fprintf(stderr,  "\t\t\t<value> is only used if population is NOT fixed\n");
    fprintf(stderr,  "\t\t\t(see '-F' above)\n");
    fprintf(stderr,  "\t\t\tDefault is %G\n\n",K_DEFAULT);
    
    fprintf(stderr,  "\t[-L <0,1>] \tOffspring born in random locations\n");
    fprintf(stderr,  "\t\t\t'-L 0': each offspring occupies the same position as\n");
    fprintf(stderr,  "\t\t\t\tone of its parents (randomly chosen).\n");
    fprintf(stderr,  "\t\t\t'-L 1': each offspring occupies a position in its\n");
    fprintf(stderr,  "\t\t\t\tnatal patch that is independent of its parents.\n");
    fprintf(stderr,  "\t\t\tDefault is %i\n\n",OIRL_DEFAULT);	
    
    fprintf(stderr,  "\t[-l <value>] \tTOTAL length of genetic map (all chromosomes combined).\n");
    fprintf(stderr,  "\t\t\tUnless a custom map is specified, each of the C\n");	
    fprintf(stderr,  "\t\t\tchromosomes will have length <value>/C\n");	
    fprintf(stderr,  "\t\t\t(see also '-C' above).\n");	
    fprintf(stderr,  "\t\t\tDefault is %G\n\n",TOTAL_MAP_LENGTH_DEFAULT);
    
    fprintf(stderr,  "\t[-m <integer>] \tNumber of mutations to introduce.\n");
    fprintf(stderr,  "\t\t\tDefault is %i\n\n",nMUTATIONS_DEFAULT);
    
    fprintf(stderr,  "\t[-N <integer>] \tPopulation size (or initial population size\n");
    fprintf(stderr,  "\t\t\tif used with '-F 0'; see above).\n");
    fprintf(stderr,  "\t\t\tDefault is %i\n\n",INITIAL_POPULATION_SIZE_DEFAULT);
    
    fprintf(stderr,  "\t[-O <0,1>] \tMosaic or gradient environmental variation\n");
    fprintf(stderr,  "\t\t\t'-O 0': habitat with a discrete 'gradient' of\n\t\t\t\tfitness variation.\n");
    fprintf(stderr,  "\t\t\t'-O 1': habitat with a mosaic of patches.\n");
    fprintf(stderr,  "\t\t\tIf the total number of patches is only 2, this\n\t\t\toption does not matter.\n");
    fprintf(stderr,  "\t\t\tDefault is %i\n\n",MOSAIC_DEFAULT);	
    
    fprintf(stderr,  "\t[-P <value>] \tPercent of S distribution tail 'fattened'.\n");
    fprintf(stderr,  "\t\t\tSee usage of '-f' above.\n\t\t\t<value> set here is only used with '-f 1'.\n");	
    fprintf(stderr,  "\t\t\tDefault is %G\n\n",FATTEN_TAIL_PROPORTION_DEFAULT);
    
    fprintf(stderr,  "\t[-R <thresh>] \tThreshold value for determining RI.\n");
    fprintf(stderr,  "\t\t\tIf avg. resident fitness > thresh * avg. immigrant\n");
    fprintf(stderr,  "\t\t\tfitness, then the program is exited.\n");
    fprintf(stderr,  "\t\t\tDefault is %E\n\n",RI_THRESH);
    
    fprintf(stderr,  "\t[-S <mean>] \tMean value of selection coefficients.\n");
    fprintf(stderr,  "\t\t\tCoefficients are drawn from an exponential\n");
    fprintf(stderr,  "\t\t\tdistribution with this as the distribution's\n");
    fprintf(stderr,  "\t\t\tparameter (but see -f above).\n");
    fprintf(stderr,  "\t\t\tDefault is %G\n\n",MEAN_S_DEFAULT);
    
    
    fprintf(stderr,  "\t[-s <0,1>] \tSymmetric S values for 0 and 1 alleles\n");
    fprintf(stderr,  "\t\t\t'-s 0': values chosen for 0 and 1 alleles at a locus\n");
    fprintf(stderr,  "\t\t\t\tare independent of each other.\n");
    fprintf(stderr,  "\t\t\t'-s 1': allele '0' is just as good in patch #0 as\n");
    fprintf(stderr,  "\t\t\t\tallele '1' is in patch #nPATCHES-1.\n");
    fprintf(stderr,  "\t\t\tDefault is %i\n\n",SYMMETRIC_MUTATIONS_DEFAULT);	
    
    fprintf(stderr,  "\t[-T <integer>] \tdesired number of time samples.  Approximate!\n");
    fprintf(stderr,  "\t\t\tDefault is %i \n\n",nTIME_SAMPLES_DEFAULT);
    
    fprintf(stderr,  "\t[-t <integer>] \tAmount of time (number of generations)\n");
    fprintf(stderr,  "\t\t\tbetween mutation introductions.\n");
    fprintf(stderr,  "\t\t\tDefault is %i generations\n\n",nGENERATIONS_MAX_DEFAULT);
    
    fprintf(stderr,  "\t[-V <0,1>] \tRecord linkage disequilibrium time series.\n");
    fprintf(stderr,  "\t\t\tBoolean variable:  0 = do not record, 1 = record.\n");
    fprintf(stderr,  "\t\t\tLinkage disequilibrium data file ('LDvalues.txt')\n");
    fprintf(stderr,  "\t\t\tproduced with -V 1 may be very large (10s of GB).\n");
    fprintf(stderr,  "\t\t\tDefault is %i\n\n",RECORD_LD_VALUES_DEFAULT);
    
}



void warmUpRNG(void)
{
    long int i;
    for ( i = 0; i < 1000000; i++ ) randU();
}

//
// printTime	-- display an elapsed time in hours, minutes, seconds
//
void
printTime(long t)
{
    int days, hours, min, sec;
    
    sec = (int)(t / CLOCKS_PER_SEC);
    days = sec / (24 * 60 * 60);
    sec  -= days * (24 * 60 * 60);
    hours = sec / (60 * 60);
    sec  -= hours * (60 * 60);
    min   = sec / 60;
    sec  -= min * 60;
    
    if (days)
        printf("%dd ", days);
    if (hours)
        printf("%dh ", hours);
    if (min)
        printf("%dm ", min);
    
    //
    // If it's a short run, print fractions of seconds
    //
    if ((days + hours + min) == 0)
        printf("%.2f sec\n", (float)t / CLOCKS_PER_SEC);
    else
        printf("%d sec\n", sec);
}


//
// Palloc/Pfree	-- specialized memory allocation. Keeps previously
//	allocated memory on a free list for quick reuse, no matter
//	what the size of the allocation.  
//
//	Must pass in a non-zero pool tag to help palloc() reuse previously allocated/freed
//	memory.  Allocations from the same pool should be of similar size/use.
//
void *
palloc(int tag, size_t size)
{
    palloc_hdr *buf, *prev;
    static int initialized = 0;
    int i;
    size_t realloc_size;
    char *palloc_env;
    
    palloc_count++;
    
    if (!initialized) {
        initialized = 1;
        
        palloc_env = getenv("PALLOC");
        if (palloc_env) {
            printf("PALLOC = %s\n", palloc_env);
            palloc_debug = atoi(palloc_env);
            printf("palloc_debug = %d\n", palloc_debug);
        }
        
        for (i = 0; i < PALLOC_TRACK_LIMIT; i++) {
            palloc_pool_track[i].tag = 0;
            palloc_pool_track[i].size = 0L;
            palloc_pool_track[i].requested = 0;
            palloc_pool_track[i].allocated = 0;
            palloc_pool_track[i].reallocated = 0;
        }
    }
    
    if (tag == 0) {				 // make sure tag is legal
        fprintf(stderr, "Must use pool tag != 0\n");
        exit(-1);
    }
    
    if (palloc_debug > 1)
        printf("palloc%d(%ld)\n", palloc_count, (long)size);
    
    //
    // See if we've had an allocation of this pool previously
    //
    for (i = 0; i < palloc_track_count; i++)
        if (tag == palloc_pool_track[i].tag)
            break;
    
    //
    // if it's a new pool, start tracking it
    //
    if (i == palloc_track_count) {
        palloc_track_count++;
        palloc_pool_track[i].tag = tag;
        palloc_pool_track[i].size = size;
    }
    
    palloc_pool_track[i].requested++;
    
    //
    // Search our free list for a previously allocated block
    // of the correct size.
    //
    if (palloc_debug > 4) {
        printf("FREELIST: ");
        buf = palloc_free_list;
        while (buf) {
            printf("%p %d	%d ", buf, (int)buf->size, buf->tag);
            buf = buf->next;
        }
        printf("\n");
    }
    
    if (palloc_debug > 3) 
        printf("Looking for reclaim for size %d\n", (int)size);
    
    prev = NULL;
    buf = palloc_free_list;
    while (buf != NULL) {
        if (buf->tag == tag)
            break;			// found one
        prev = buf;
        buf = buf->next;
    }
    
    if (buf) {				// if we found one, udpate
        if (palloc_debug > 3) 
            printf("success: %p %d   %d\n", buf, (int)buf->size, buf->tag);
        palloc_free--;			// stats
        palloc_reclaimed++;
        
        if (prev) 			// pull it out of the list
            prev->next = buf->next;
        else
            palloc_free_list = buf->next;
        
        if (buf->size < size) {		// make sure it's big enough
            //
            // make the new buffer size 10% larger than what was
            // requested so we avoid excess reallocations
            //
            realloc_size = (size_t)((size * 11) / 10);
            
            //printf("before realloc: buf = %p ", buf);
            buf = (palloc_hdr *)realloc(buf, realloc_size + PALLOC_HDR_SIZE);
            //printf("after: buf = %p\n", buf);
            //fflush(NULL);
            buf->size = realloc_size;
            buf->tag = tag;
            
            palloc_realloced++;
            palloc_pool_track[i].reallocated++;
        }
        
        buf->next = NULL;
    }
    else {
        if (palloc_debug > 3) 
            printf("failed\n");
        
        //
        // need to allocate a new buffer with our headed prepended
        //
        if (++palloc_active > palloc_max_active)
            palloc_max_active = palloc_active;
        
        palloc_pool_track[i].allocated++;
        buf = (palloc_hdr *)malloc(size+PALLOC_HDR_SIZE);
        //printf("palloc(%d): buf = %p\n", size+PALLOC_HDR_SIZE, buf);
        buf->size = size;
        buf->tag = tag;
        buf->next = NULL;
    }
    
    return (void *)((unsigned char *)buf + PALLOC_HDR_SIZE);
}


void 
pfree(void *ptr)
{
    int i;
    palloc_hdr *buf;
    
    palloc_active--;
    if (++palloc_free > palloc_max_free)
        palloc_max_free = palloc_free;
    
    //
    // Adjust the pointer to get back to our header.
    //
    buf = (palloc_hdr *)((unsigned char *)ptr - PALLOC_HDR_SIZE);
    
    //
    // Add to head of our free list for later reuse.
    //
    buf->next = palloc_free_list;
    palloc_free_list = buf;
    
    if (palloc_debug > 2) 
        printf("freeing size %d\n", (int)buf->size);
    
    if (palloc_debug > 4) {
        printf("FREELIST: ");
        buf = palloc_free_list;
        while (buf) {
            printf("%d ", (int)buf->size);
            buf = buf->next;
        }
        printf("\n");
    }
    
    //
    // Don't actually free the buffer, but this is what we'd do
    // if we did.
    //
    //	free((void *)((unsigned char *)ptr - PALLOC_HDR_SIZE));
}


void
palloc_stats()
{
    int i;
    
    if (palloc_debug) {
        printf("pool	size	alloc	/ request	realloc\n");
        for (i = 0; i < palloc_track_count; i++)
            printf("%d	%ld	%d	/ %d		%d\n", i, 
                   palloc_pool_track[i].size,
                   palloc_pool_track[i].allocated,
                   palloc_pool_track[i].requested,
                   palloc_pool_track[i].reallocated);
        
        printf("palloc_count = %d\n", palloc_count);
        printf("palloc_max_active = %d\n", palloc_max_active);
        printf("palloc_max_free = %d\n", palloc_max_free);
        printf("palloc_reclaimed = %d\n", palloc_reclaimed);
        printf("palloc_realloced = %d\n", palloc_realloced);
    }
}
