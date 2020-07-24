function [missingDirs, missingParams, missingFinals] = make_parameter_summary(indexes, workingdir, prefix)

orig = pwd;
cd(workingdir);
if nargin < 3
    prefix = 'Run';
end



[indexes, errs, missingDirs, missingParams, missingFinals] = filterDirectories(indexes, prefix);

n = numel(indexes);

fname = ['RunSummary_' num2str(indexes(1)) '-' num2str(indexes(n)) '.csv'];
fpt = fopen(fname,'w');

fprintf(fpt,['run_index,RnumSeed,DH_GH_or_direct,D,PATCHES,nPATCHES,' ...
    'USE_MUTATIONS_FROM_FILES,ADDITIVE_FITNESS,MULTIPLICATIVE_FITNESS,' ...
    'TWO_DEME,DETERMINISTIC,nGENERATIONS_MAX,nMUTATIONS,' ...
    'INITIAL_CONDITIONS,INITIAL_POPULAION_SIZE,MEAN_S,' ...
    'SYMMETRIC_MUTATIONS,MAP_TYPE,GAMETE_PRODUCTION_MODE,FIXED_N,' ...
    'K,H,MOSAIC,SD_MOVE,OFFSPRING_IN_RANDOM_LOCATIONS,' ...
    'nCHROMOSOMES,TOTAL_MAP_LENGTH,' ...
    'MUTATION_DISTRIBUTION,FATTEN_TAIL_PROPORTION,FATTEN_TAIL_MAX,' ...
    'NUMBER_BIG_TO_ESTABLISH, START_WITH_BIG_IN_PLACE,BIG_START_S_VAL,' ...
    'END_PERIOD_ALLOPATRY,START_PERIOD_ALLOPATRY,PERIOD_ALLOPATRY,' ...
    'MUTATIONS_PER_GENERATION,' ...
    'DEME0_CONSTANT_S,DEME1_CONSTANT_S,DMI_MEAN_EFFECT,DMI_PROBABILITY,' ...
    'PROBABILITY_DERIVED_REVERSED,POSITIVE_EPI_PROBABILITY,' ...
    'MEAN_POSITIVE_EPI_COEFF,TS_SAMPLING_FREQUENCY,RECORD_LD_VALUES,' ...
    'START_THRESH_FOR_LD,END_THRESH_FOR_LD,LD_LOCI_SUBSAMPLE,' ...
    'LD_LowerBound,RI_THRESH,FST_MIN_RECORDING_THRESH,' ...
    'FRACTION_SELECTED_LOCI,CodeVersion,RI_REACHED,LDfilesPresent,TotalFileCount\n']);
for i = 1:n
    wi = indexes(i);
    wd = [prefix num2str(wi)];
    cd(wd);
    % some defaults;
    END_PERIOD_ALLOPATRY = nan;
    START_PERIOD_ALLOPATRY = nan;
    PERIOD_ALLOPATRY = nan;
    MUTATIONS_PER_GENERATION = nan;
    DEME0_CONSTANT_S = nan;
    DEME1_CONSTANT_S = nan;
    DMI_MEAN_EFFECT = nan;
    FST_MIN_RECORDING_THRESH = nan;
    FRACTION_SELECTED_LOCI = nan;
    LD_LowerBound = nan;
    RI_REACHED = nan;
        parameters;
        if exist('RnumSeed.txt.bz2','file')
            system('bunzip2 RnumSeed.txt.bz2');
        end
        rnumseed = importdata('RnumSeed.txt');
        if GAMETE_PRODUCTION_MODE == 2
            gpm = 'direct';
        elseif GAMETE_PRODUCTION_MODE == 1
            gpm = 'GH';
        else
            gpm = 'DH';
        end
        
        
        printString = ['%i,%i,%s,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,' ... % through nMUTATIONS
            '%i,%i,%f,%i,%i,%i,%i,%f,%f,%i,%f,%i,%i,%f,%i,%f,%f,' ... % through FATTEN_TAIL_MAX
            '%i,%i,%f,%i,%i,%i,%i,%f,%f,%f,%f,%f,%f,%f,%i,%i,' ... % through RECORD_LD_VALUES
            '%f,%f,%i,%f,%f,%f,%f,%s,%i,'];
        fprintf(fpt,printString, ... 
            wi,rnumseed,gpm,D,PATCHES, ...
            nPATCHES,USE_MUTATIONS_FROM_FILES,ADDITIVE_FITNESS, ...
            MULTIPLICATIVE_FITNESS,TWO_DEME,DETERMINISTIC, ...
            nGENERATIONS_MAX,nMUTATIONS,INITIAL_CONDITIONS, ...
            INITIAL_POPULATION_SIZE,MEAN_S,SYMMETRIC_MUTATIONS, ...
            MAP_TYPE,GAMETE_PRODUCTION_MODE,FIXED_N,K,H(1),MOSAIC, ...
            SD_MOVE,OFFSPRING_IN_RANDOM_LOCATIONS, ...
            nCHROMOSOMES,TOTAL_MAP_LENGTH,MUTATION_DISTRIBUTION, ...
            FATTEN_TAIL_PROPORTION,FATTEN_TAIL_MAX,NUMBER_BIG_TO_ESTABLISH, ...
            START_WITH_BIG_IN_PLACE,BIG_START_S_VAL, ...
            END_PERIOD_ALLOPATRY,START_PERIOD_ALLOPATRY,PERIOD_ALLOPATRY, ...
            MUTATIONS_PER_GENERATION, ...
            DEME0_CONSTANT_S,DEME1_CONSTANT_S,DMI_MEAN_EFFECT, ... 
            DMI_PROBABILITY,PROBABILITY_DERIVED_REVERSED, ...
            POSITIVE_EPI_PROBABILITY,MEAN_POSITIVE_EPI_COEFF, ...
            TS_SAMPLING_FREQUENCY,RECORD_LD_VALUES,START_THRESH_FOR_LD, ...
            END_THRESH_FOR_LD,LD_LOCI_SUBSAMPLE,LD_LowerBound,RI_THRESH, ...
            FST_MIN_RECORDING_THRESH,FRACTION_SELECTED_LOCI,CodeVersion,RI_REACHED);
        
        if exist('LDneutSitesDiff.txt','file') || exist('LDneutSitesDiff.txt.bz2','file')
            fprintf(fpt,'1,');
        else
            fprintf(fpt,'0,');
        end
        
        [~,fileCount] = system('ls | wc -l');
        fileCount = str2double(fileCount);
        fprintf(fpt, '%i\n', fileCount);
        
    cd ..;
end



fclose(errs);
fclose(fpt);

cd(orig);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [indexes, errs, missingDirs, missingParams, missingFinals] = filterDirectories(indexes, prefix)

if nargin < 2
    prefix = 'Run';
end

ninit = numel(indexes);
firsti = indexes(1);
lasti = indexes(ninit);
filename = ['Msgs_ParamSummary_Runs' num2str(firsti) '-' num2str(lasti) '.txt'];
errs = fopen(filename,'w');

missingDirs = zeros(1,ninit);
mdcount = 0;
missingParams = missingDirs;
mpcount = 0;
missingFinals = missingDirs;
mfcount = 0;

% filter out nonexistent runs
for i = ninit:-1:1
    dirname = [prefix num2str(indexes(i))];
    if ~exist(dirname,'dir')
        fprintf(errs,'Directory %s does not exist. Removing...\n', dirname);
        mdcount = mdcount + 1;
        missingDirs(mdcount) = indexes(i);        
        indexes(i) = [];
    else
        pname = [dirname '/parameters.m'];
        pzname = [pname '.bz2'];
        finname = [dirname '/FinalLocations.txt'];
        zfinname = [dirname '/FinalLocations.txt.bz2'];
        if exist(pzname,'file')
            syscall = ['bunzip2 ' pzname];
            system(syscall);
        elseif ~exist(pname,'file')
            fprintf(errs,'Parameters file %s does not exist!  Removing...\n', pname);
            mpcount = mpcount + 1;
            missingParams(mpcount) = indexes(i);
            indexes(i) = [];
        elseif ~exist(finname,'file') && ~exist(zfinname,'file')
            fprintf(errs,'FinalLocations file %s does not exist!  Removing...\n', finname);
            mfcount = mfcount + 1;
            missingFinals(mfcount) = indexes(i);
            indexes(i) = [];
        end
    end
end
if mdcount < 1
    missingDirs = 0;
else
    missingDirs = missingDirs(1:mdcount);
end
if mpcount < 1
    missingParams = 0;
else
    missingParams = missingParams(1:mpcount);
end
if mfcount < 1
    missingFinals = 0;
else
    missingFinals = missingFinals(1:mfcount);
end

if mdcount > 0
    rdf = fopen('AbsentDirs.txt','w');
    for i = 1:mdcount
        fprintf(rdf, '%i ', missingDirs(i));
    end
    fprintf(rdf,'\n');
    fclose(rdf);
end
if mpcount > 0
    rdf = fopen('RedoAbsentParams.txt','w');
    for i = 1:mpcount
        fprintf(rdf, '%i ', missingParams(i));
    end
    fprintf(rdf,'\n');
    fclose(rdf);
end
if mfcount > 0
    rdf = fopen('RedoAbsentFinals.txt','w');
    for i = 1:mfcount
        fprintf(rdf, '%i ', missingFinals(i));
    end
    fprintf(rdf,'\n');
    fclose(rdf);
end

n = numel(indexes);
%report on step one
msg = ['gathering info for ' num2str(n) ' directories with parameters.m and FinalLocations.txt files'];
disp(msg);
fprintf(errs, '%s\n', msg);

