# Linkage disequilibrium
is calculated in several ways and places in the code.  There are also a number of controls on it.  This creates a bit of a mess.  The following is a description of relevant functions and the controls on data collection.

## I. Summary of command line control
The most important control variable is `gatherLDvalues`, which is an int that can be set from the command line with -g
* if `-g = 0`, there won't be much LD data
* if `-g = 1`, there will be a lot of average LD data
* if `-g = 3`, there will be tons of site specific data (potentially gigabytes).  Use with caution.

`RECORD_LD_VALUES` is a legacy control that still plays some role.  `gatherLDvalues` was added because finer control was desired.  `RECORD_LD_VALUES` is automatically set to a value of "1" (i.e., true) in `main()` if `gatherLDvalues >= 3`.



## II. Where LD calculations occur in the code:

### II. A. The bulk of the calculations are called from `calcExpectedME()` and depend upon `gatherLDvalues`

Note that `calcExpectedME()` is called from `reproduce()`, and `reproduce()` is called from `main()`

Within `calcExpectedME()`, we find:
* Line 4108-4123: if `RECORD_LD_VALUES == 1`, then we can either start or stop recording some values with `BeginRecordingLD` being changed to a value of 0 or 1.
* Separately, on lines 4125-4129, we have:
```
if ( gatherLDvalues >= 1 ) {
    calculateLDselectedSitesOnly( gatherLDvalues );
    calculateLDneutralSitesOnly( gatherLDvalues );
    calculateLD( gatherLDvalues );
}
```
* `calculateLDselectedSitesOnly()` and `calculateLDneutralSitesOnly()` are nearly identical in their operations, except the former works on selected sites and the latter on neutral sites.
    * When either of these is called, they make a call to `calculateLDpairOneOff()`
    * returned value from that function is used to calculate an average value.  These average values are written to `LDneutSitesAvg` and `LDselSitesAvg`
    * the returned value is ALSO individually written to a data file if `gatherLDvalues >= 3`.  These individual values are written to `LDselSitesSame`, `LDselSitesDiff`, `LDneutSitesSame`, and `LDneutSitesDiff`

* `calculateLD()` calls `calculateLDpair(l1,l2,dist, dpt1, dpt2, dpt3, gatherLDvalues)`
    * `calculateLD()` writes averages and variances to `effMigRates`
    * `calculateLDpair()` writes individual values to `LDfpt` IF `gatherLDvalues >= 3` AND `BeginRecordingLD` AND `fabs(DD) > LD_LowerBound`


### II. B. Additional calculations of LD with nearest selected neighbor are called from `calculateAlleleFrequencies`

* These calls do NOT depend on `gatherLDvalues`.
* `calculateAlleleFrequencies()` is called from `main()`
* nsnLD values are obtained from `calculateLDpairOneOff()`
* nsnLD values are written to `selectedFrequencies` or `neutralFrequencies`


## III. Which files have which data?
| File Pointer | File Name | Condition of writing data | Description of data |
| ------------- | ----------- | ----------------------------- | --------------------- |
| `selectedFrequencies` | "SelectedAlleleFrequencies.txt" | always when `calculateAlleleFrequencies()` is called | nsnLD |
| `neutralFrequencies` | "SelectedAlleleFrequencies.txt" | always when `calculateAlleleFrequencies()` is called | nsnLD |




## Additional random notes

`alleleFrequenciesByPatch[nLOCI][nPATCHES]` is a local variable created in `calculateAlleleFrequencies()`

Control via `gatherLDvalues ( -g <int>)`
+ local variable created in Main --> sent to reproduce --> sent to calcExpectedME
+ LINE 560:
```
if ( gatherLDvalues >= 3 )
    RECORD_LD_VALUES = 1;
```
+ LINE 763:
``` if ( gatherLDvalues >= 1 ) {
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
```
+ used by:
openDataFilesForRecording( gatherLDvalues );
reproduce(gatherLDvalues);
calculateLD( gatherLDvalues );
+ LINE 1721:  if ( gatherLDvalues >= 1 ) in calculateLD()
then D, Dprime, and Delta are calculated and written to effMigRates
for LD_LOCI_SUBSAMPLE
line 1744: sent to calculateLDpair(l1,l2,dist, dpt1, dpt2, dpt3, gatherLDvalues)

+ in calculateLDneutralSitesOnly()
line 1816: if ( gatherLDvalues >= 3 )
then individual site data is recorded
same for calculateLDselectedSitesOnly() line 1890
same for calculateLDpair, line 1983






Control via RECORD_LD_VALUES
+ set with -V
+ automatically set to 1 if gatherLDvalues >= 3 (line 517)
+ can cause BeginRecordingLD = 1 (line 4118) in conjunction with START_THRESH_FOR_LD

Control via BeginRecordingLD
affects calculateLDpair(), line 1984; that is ONLY effect: whether individual site data are written to LDfpt
set at line 4115, but this is unnecessary if setTSfreq() has been called:
else { // check to see if we should start recording
if ( migrationCount[0] && migrationCount[j] ) {
if ( (avgImmFit[0]/avgResFit[0] <= START_THRESH_FOR_LD) || (avgImmFit[j]/avgResFit[j] <= START_THRESH_FOR_LD) ) {
BeginRecordingLD = 1;
fprintf(stderr, "\nStarted recording LD at generation %li\n", totalGenerationsElapsed);
}
}
}
IF RECORDING_TIMES_IN_FILE: then set in setTSfreq():  line 5398: BeginRecordingLD = 1;

Control via LD_LowerBound
affects calculateLDpair(), line 1985
if ( fabs(DD) > LD_LowerBound )
hard coded as 0.001
only in play if gatherLDvalues >= 3
in which case LD values are printed to LDfpt

Control via END_THRESH_FOR_LD
line 4108:
if ( RECORD_LD_VALUES && (totalGenerationsElapsed > END_PERIOD_ALLOPATRY) ) {
if ( BeginRecordingLD ) { // we already started recorded before now
if ( (avgImmFit[0]/avgResFit[0] <= END_THRESH_FOR_LD) && migrationCount[0] && (avgImmFit[j]/avgResFit[j] <= END_THRESH_FOR_LD) && migrationCount[j] && !RECORDING_TIMES_IN_FILE ) {
RECORD_LD_VALUES = 0;
fprintf(stderr, "\nStopped recording LD at generation %li\n", totalGenerationsElapsed);
}
}


Control via LD_LOCI_SUBSAMPLE:
calculateLDneutralSitesOnly
calculateLDselectedSitesOnly
calculateLD

