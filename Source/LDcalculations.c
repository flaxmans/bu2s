// File for calculations of LD involved in BU2S model.  Older calculations may
// still exist in bu2s.c
// This set of functions is intended to eventually supercede those older functions

// includes:
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "bu2s.h"

// magic number codes for LD pair types:
// magic number codes for LD pairs:
#define LD_PAIR_CODE_BOTH_SEL_SAME_CHROM 0
#define LD_PAIR_CODE_BOTH_SEL_DIFF_CHROM 1
#define LD_PAIR_CODE_BOTH_NEUT_SAME_CHROM 2
#define LD_PAIR_CODE_BOTH_NEUT_DIFF_CHROM 3
#define LD_PAIR_CODE_SEL_NEUT_SAME_CHROM 4
#define LD_PAIR_CODE_SEL_NEUT_DIFF_CHROM 5
#define NUM_LD_PAIR_CATEGORIES 6 // number different kinds of pairs
// the following names should be in exactly the order given by the magic number codes above
const char *LD_VALS_FILE_NAMES[NUM_LD_PAIR_CATEGORIES] = {"LD_vals_bothSelSameChrom.csv","LD_vals_bothSelDiffChrom.csv","LD_vals_bothNeutSameChrom.csv","LD_vals_bothNeutDiffChrom.csv", "LD_vals_selNeutSameChrom.csv", "LD_vals_selNeutDiffChrom.csv"};
const char *LD_AVG_FILE_NAMES[(NUM_LD_PAIR_CATEGORIES + 1)] = {"LD_avg_bothSelSameChrom.csv","LD_avg_bothSelDiffChrom.csv","LD_avg_bothNeutSameChrom.csv","LD_avg_bothNeutDiffChrom.csv", "LD_avg_selNeutSameChrom.csv", "LD_avg_selNeutDiffChrom.csv", "LD_avg_GenomeWide.csv"};

// file pointers
FILE *LD_filePointers[NUM_LD_PAIR_CATEGORIES];
FILE *LD_averagesFilePointers[(NUM_LD_PAIR_CATEGORIES + 1)];

void calculateLDwithinDemes( int gatherLDvalues, double *alleleFrequenciesByPatch, int *patchAbundancePtr )
{
    int i, j, patch, haplo, iLocusType, jLocusType, chromi, chromj;
    long int locusID1, locusID2, *lociIDs, *lipt;
    int *pairCategories, *ipt, *locusIndexes, *locuspt;
    double *afbppti, *afbpptj, p1, p2, *individualLDvals, *dpt, *distPtr, *distances;
    long int pairCount = 0;
    int categoryCounts[NUM_LD_PAIR_CATEGORIES];
    double averageLDvalues[(NUM_LD_PAIR_CATEGORIES+1)][(nPATCHES+1)];
    
    individualLDvals = (double *) malloc( ((nLOCI * (nLOCI - 1) / 2) * (nPATCHES + 1) * sizeof(double)) );
    // to store results for each pair of loci in each patch AND an overall value (between deme)
    pairCategories = (int *) malloc( ((nLOCI * (nLOCI - 1) / 2) * (nPATCHES + 1) * sizeof(int)) );
    lociIDs = (long int *) malloc( ((nLOCI * (nLOCI - 1) ) * (nPATCHES + 1) * sizeof(long int)) ); // no divided by 2 because we need two entries for each
    locusIndexes = (int *) malloc( ((nLOCI * (nLOCI - 1) ) * (nPATCHES + 1) * sizeof(int)) ); // no divided by 2 because we need two entries for each
    distances = (double *) malloc( ((nLOCI * (nLOCI - 1) / 2) * (nPATCHES + 1) * sizeof(double)) );
    
    
    for ( i = 0; i < NUM_LD_PAIR_CATEGORIES; i++ )
        categoryCounts[i] = 0;
    
    lipt = lociIDs;
    dpt = individualLDvals;
    ipt = pairCategories;
    distPtr = distances;
    locuspt = locusIndexes;
    
    // average LD values WITHIN demes for neutral-neutral, selected-selected
    // neutral-selected; all within chromosome and between chromosome
    
    for ( i = 0; i < (nLOCI - 1); i++ ) {
        locusID1 = locusID[i];
        if ( variable_loci[i] && (locusID1 % LD_LOCI_SUBSAMPLE == 0) ) {
            iLocusType = IS_SELECTED_LOCUS[i];
            chromi = chromosomeMembership[i];
            afbppti = alleleFrequenciesByPatch + ( i * nPATCHES ); // pointer to frequency of allele 1 of locus i in patch 0
            for ( j = (i + 1); j < nLOCI; j++ ) {
                locusID2 = locusID[j];
                if ( variable_loci[j] && (locusID2 % LD_LOCI_SUBSAMPLE == 0) ) {
                    
                    jLocusType = IS_SELECTED_LOCUS[j];
                    chromj = chromosomeMembership[j];
                    
                    afbpptj = alleleFrequenciesByPatch + ( j * nPATCHES ); // pointer to frequency of allele 1 of locus i in patch 0
                    
                    // now calculate:
                    calculateLDonePairInDemes( i, j, afbppti, afbpptj, dpt, patchAbundancePtr );
                    
                    // now classify:
                    if ( chromi == chromj ) { // on same chrom
                        if ( iLocusType && jLocusType ) {
                            // they are both selected loci
                            *ipt = LD_PAIR_CODE_BOTH_SEL_SAME_CHROM;
                        }
                        else if ( iLocusType || jLocusType ) {
                            // same chrom, one sel, one not (note previous if)
                            *ipt = LD_PAIR_CODE_SEL_NEUT_SAME_CHROM;
                        }
                        else {
                            // both neut
                            *ipt = LD_PAIR_CODE_BOTH_NEUT_SAME_CHROM;
                        }
                        *distPtr = fabs(MAP[j] - MAP[i]);
                    }
                    else {  // on different chroms
                        if ( iLocusType && jLocusType ) {
                            // they are both selected loci
                            *ipt = LD_PAIR_CODE_BOTH_SEL_DIFF_CHROM;
                        }
                        else if ( iLocusType || jLocusType ) {
                            // same chrom, one sel, one not (note previous if)
                            *ipt = LD_PAIR_CODE_SEL_NEUT_DIFF_CHROM;
                        }
                        else {
                            // both neut
                            *ipt = LD_PAIR_CODE_BOTH_NEUT_DIFF_CHROM;
                        }
                        *distPtr = NAN;
                    }
                    categoryCounts[ (*ipt) ] += 1; // keep track of how many pairs of each type exist
                    pairCount++; // one more pair taken care of; will store total number of pairs by the end
                    
                    ipt++; // increment pairCategory pointer
                    *lipt = locusID1; // Locus ID pointer
                    lipt++;
                    *lipt = locusID2;
                    lipt++;
                    *locuspt = i;  // locus index pointer
                    locuspt++;
                    *locuspt = j;
                    locuspt++;
                    distPtr++; // distance array pointer
                    dpt += (nPATCHES + 1); // increment pointer to LD values storage array
                    
                }  // end on if check on conditional use of second locus; nothing done otherwise
                
            } // end of loop over second member of pair
            
        } // end of conditional checks on validity of first locus for use
        
    } // end of for loop over first locus
    
    // calculate averages and write to data files:
    calculateLDaverageVals( pairCount, categoryCounts, pairCategories, individualLDvals, &averageLDvalues[0][0] );
    
    // write individual values to data files:
    ipt = pairCategories;
    dpt = individualLDvals;
    distPtr = distances;
    lipt = lociIDs;
    locuspt = locusIndexes;
    int dumi;
    if ( gatherLDvalues >= 3 ) {
        for ( i = 0; i < pairCount; i++ ) {
            dumi = *ipt;
            fprintf(LD_filePointers[dumi], "%li,%li,%li", totalGenerationsElapsed, *lipt, *(lipt + 1));
            for ( j = 0; j < (nPATCHES + 1); j++ ) {
                fprintf(LD_filePointers[dumi], ",%E", *dpt);
                dpt++;
            }
            if ( dumi == LD_PAIR_CODE_BOTH_SEL_SAME_CHROM || dumi == LD_PAIR_CODE_BOTH_NEUT_SAME_CHROM || dumi == LD_PAIR_CODE_SEL_NEUT_SAME_CHROM ) {
                fprintf(LD_filePointers[dumi], ",%E", *distPtr);
            }
            if ( dumi != LD_PAIR_CODE_BOTH_NEUT_SAME_CHROM && dumi != LD_PAIR_CODE_BOTH_NEUT_DIFF_CHROM ) {
                fprintf(LD_filePointers[dumi], ",%E,%E", S_MAX1[(*locuspt)], S_MAX1[(*(locuspt + 1))]);
            }
            fprintf(LD_filePointers[dumi], "\n");
            ipt++;
            distPtr++;
            lipt += 2;
            locuspt += 2;
        }
    }
    FILE *dumFpt;
    if ( gatherLDvalues >= 1 ) {
        for ( i = 0; i < (NUM_LD_PAIR_CATEGORIES + 1); i++ ) {
            dumFpt = LD_averagesFilePointers[i];
            fprintf(dumFpt, "%li", totalGenerationsElapsed);
            for ( j = 0; j < (nPATCHES + 1); j++ ) {
                fprintf(dumFpt, ",%E", averageLDvalues[i][j]);
            }
            fprintf(dumFpt, "\n");
        }
    }
    
    
    // free up memory grabbed using malloc:
    free(individualLDvals);
    free(pairCategories);
    free(lociIDs);
    free(distances);
    free(locusIndexes);
}

void calculateLDonePairInDemes( int l1, int l2, double *afbppti, double *afbpptj, double *LDstorageArray, int *patchAbundancePtr ) {
    
    double p1, p2, *dpt = LDstorageArray, Ninv;
    int patch, oo[(nPATCHES + 1)], oz[(nPATCHES + 1)], zo[(nPATCHES + 1)], zz[(nPATCHES + 1)];
    int haplo, n_count = 0, ii;
    short int *gpt1, *gpt2;
    double obs10, obs11, obs00, obs01, DD, corrAlleleStates;
    
    // initialize counts:
    for ( ii = 0; ii < (nPATCHES + 1); ii++ ) {
        oo[ii] = 0;
        oz[ii] = 0;
        zo[ii] = 0;
        zz[ii] = 0;
    }
    
    // loop over patches
    for ( patch = 0; patch < nPATCHES; patch++ ) {
        p1 = *(afbppti + patch); // locus 1, "1" allele frequency in patch
        p2 = *(afbpptj + patch); // locus 2, "1" allele frequency in patch
        if ( p1 > 0.0 && p1 < 1.0 && p2 > 0.0 && p2 < 1.0 ) {
            // both loci are variable in this patch
            for ( haplo = 0; haplo < 2; haplo++ ) {  // loop over each haplotype
                
                gpt1 = genotypes + ( 2 * nLOCI * n_count ) + (2 * l1) + haplo; // allele of first locus on chromosome for first individual in patch
                gpt2 = genotypes + ( 2 * nLOCI * n_count ) + (2 * l2) + haplo; // allele on second locus on chromosome
                
                for ( ii = 0; ii < patchAbundancePtr[patch]; ii++ ) {
                    if ( *gpt1 ) {
                        if ( *gpt2 ) {
                            oo[patch]++;
                        }
                        else {
                            oz[patch]++;
                        }
                    }
                    else {
                        if ( *gpt2 ) {
                            zo[patch]++;
                        }
                        else {
                            zz[patch]++;
                        }
                    }
                    gpt1 += (2 * nLOCI); // next individual
                    gpt2 += (2 * nLOCI);
                }
                
            }
            
            Ninv = 0.5 / ((double) patchAbundancePtr[patch]);
            // haplotype frequencies
            obs11 = ((double) oo[patch]) * Ninv; // observed 11
            obs10 = ((double) oz[patch]) * Ninv; // observed 10
            obs01 = ((double) zo[patch]) * Ninv; // observed 01
            obs00 = ((double) zz[patch]) * Ninv; // observed 00
            
            // LD coeff:
            DD = (( obs11 * obs00 ) - ( obs01 * obs10 ));
            // correlation of allelic states:
            *dpt = DD / sqrt( p1 * p2 * (1.0 - p1) * (1.0 - p2) );
            
            // add to total count over ALL patches:
            oo[nPATCHES] += oo[patch];
            oz[nPATCHES] += oz[patch];
            zo[nPATCHES] += zo[patch];
            zz[nPATCHES] += zz[patch];
        }
        else {
            // LD is zero in this patch
            *dpt = 0.0;
        }
        n_count += patchAbundancePtr[patch]; // add in those already done
        dpt++; // increment pointer to storage array for next deme
    }
    
    // overall LD:
    Ninv = 0.5 / ((double) N);
    obs11 = ((double) oo[nPATCHES]) * Ninv; // observed 11
    obs10 = ((double) oz[nPATCHES]) * Ninv; // observed 10
    obs01 = ((double) zo[nPATCHES]) * Ninv; // observed 01
    obs00 = ((double) zz[nPATCHES]) * Ninv; // observed 00
    
    // LD coeff:
    DD = (( obs11 * obs00 ) - ( obs01 * obs10 ));
    
    p1 = allele_frequencies[l1];
    p2 = allele_frequencies[l2];
    // correlation of allelic states:
    *dpt = DD / sqrt( p1 * p2 * (1.0 - p1) * (1.0 - p2) );
    
}


void calculateLDaverageVals( long int pairCount, int *categoryCounts, int *pairCategories, double *individualLDvals, double *averageLDvalues )
{
    // note averageLDvalues[(NUM_LD_PAIR_CATEGORIES + 1)][(nPATCHES+1)];
    int j, row, *ipt = pairCategories;
    long int i;
    double categorySums[(NUM_LD_PAIR_CATEGORIES + 1)][(nPATCHES + 1)];
    double *dpt, denom;
    int testCountSum[(NUM_LD_PAIR_CATEGORIES + 1)];
    // the +1's on NUM_LD_PAIR_CATEGORIES are for genome-wide averages across all locus-pair types
    
    // initialize arrays
    for ( i = 0; i < (NUM_LD_PAIR_CATEGORIES + 1); i++ ) {
        testCountSum[i] = 0; // testing piece
        for ( j = 0; j < (nPATCHES + 1); j++ ) {
            categorySums[i][j] = 0.0;
        }
    }
    
    // calculate sums for each category and each patch, plus the aggregate
    dpt = individualLDvals;
    for ( i = 0; i < pairCount; i++ ) {
        // loop over categories first
        row = *ipt;  // storage row is category number code
        for ( j = 0; j < (nPATCHES + 1); j++ ) {
            // loop over patches
            categorySums[row][j] += *dpt;
            categorySums[NUM_LD_PAIR_CATEGORIES][j] += *dpt; // genome wide numbers in last row
            dpt++; // increment pointer to LD values
        }
        testCountSum[(*ipt)] += 1;
        ipt++; // increment pointer to pair categories
    }
    
    //testing piece:
    for ( i = 0; i < NUM_LD_PAIR_CATEGORIES; i++ ) {
        if ( categoryCounts[i] != testCountSum[i] ) {
            fprintf(stderr, "\nError in calculateLDaverageVals():\n\tcategoryCounts[%li] (%i) != testCountSum[%li] (%i)\n\n", i, categoryCounts[i], i, testCountSum[i]  );
            exit(-1);
        }
    }
    
    // calculate the averages
    dpt = averageLDvalues;
    for ( i = 0; i < NUM_LD_PAIR_CATEGORIES; i++ ) {
        denom = 1.0 / ((double) categoryCounts[i]); // denominator for calculating means
        for ( j = 0; j < (nPATCHES + 1); j++ ) {
            *dpt = categorySums[i][j] * denom;
            dpt++;
        }
    }
    // genome wide uses different denominator, but continues from same pointer
    denom = 1.0 / ((double) pairCount);
    for ( i = 0; i < (nPATCHES + 1); i++ ) {
        *dpt = categorySums[NUM_LD_PAIR_CATEGORIES][i] * denom;
        dpt++;
    }
}


void openLDdataFiles( int gatherLDvalues )
{
    int i, j;
    
    if ( gatherLDvalues >= 3 ) {
        for ( i = 0; i < NUM_LD_PAIR_CATEGORIES; i++ ) {
            LD_filePointers[i] = fopen( LD_VALS_FILE_NAMES[i], "w" );
            
            // make headers:
            fprintf(LD_filePointers[i], "totalGenerationsElapsed,locusID1,locusID2");
            for  ( j = 0; j < nPATCHES; j++ )
                fprintf( LD_filePointers[i], ",LDinDeme%i", j);
            fprintf(LD_filePointers[i], ",globalLD)");
            if ( i == LD_PAIR_CODE_BOTH_SEL_SAME_CHROM || i == LD_PAIR_CODE_BOTH_NEUT_SAME_CHROM || i == LD_PAIR_CODE_SEL_NEUT_SAME_CHROM ) {
                fprintf(LD_filePointers[i], ",distance");
            }
            if ( i != LD_PAIR_CODE_BOTH_NEUT_SAME_CHROM && i != LD_PAIR_CODE_BOTH_NEUT_DIFF_CHROM ) {
                fprintf(LD_filePointers[i], ",Scoeff1,Scoeff2" );
            }
            fprintf(LD_filePointers[i], "\n");
        }
    }
    if ( gatherLDvalues >= 1 ) {
        for ( i = 0; i < (NUM_LD_PAIR_CATEGORIES + 1); i++ ) {
            LD_averagesFilePointers[i] = fopen( LD_AVG_FILE_NAMES[i], "w" );
            
            // make headers
            fprintf(LD_averagesFilePointers[i], "totalGenerationsElapsed");
            for ( j = 0; j < nPATCHES; j++ )
                fprintf(LD_averagesFilePointers[i], ",avgInDeme%i", j);
            fprintf(LD_averagesFilePointers[i], ",globalAvg\n");
        }
    }
}

void closeLDdataFiles( int gatherLDvalues )
{
    int i;
    
    if ( gatherLDvalues >= 3 ) {
        for ( i = 0; i < NUM_LD_PAIR_CATEGORIES; i++ )
            fclose( LD_filePointers[i] );
    }
    if ( gatherLDvalues >= 1 ) {
        for ( i = 0; i < (NUM_LD_PAIR_CATEGORIES + 1); i++ )
            fclose( LD_averagesFilePointers[i] );
    }
}





