%%%%%
%%%%% PURPOSE: Once ran /PSID/prepData2023.m, now need to create
%%%%%       pseudopopulations of neurons/trials to run PSID training/test on 
%%%%% Allows you to balance # neural features across cortical regions
%%%%%
%%%%% This will load '[region]_allTrials.mat', 
%%%%%
%%%%% OUTPUT: 'pp_NN.mat'; 1 per pseudopop iteration 
%%%%%           Contains, in addition to adjuvent info:
%%%%%           PSIDpseudopop -- 1 field per subject (B, T)
%%%%%               Each contains 5 structures -- 1 per region, plus 'All'
%%%%%               regions together 
%%%%%               Each region field contains:
%%%%%                    'neural' -- nTrials x 1 cellarray, each containing
%%%%%                    nTimepts x nFeatures; NOTE all trials from a given
%%%%%                    recording date contain same sampling of neurons,
%%%%%                    but across recording dates (and across
%%%%%                    pseduopopulations/matfiles) the sampled neurons
%%%%%                    differ. Said differently -- across
%%%%%                    pseudopopulations, the behavioral trials are all
%%%%%                    the same (you use every trial across every day,
%%%%%                    that fits criteria, it's just the NEURON RESPONSES that
%%%%%                    differ)
%%%%%               For the combinedSubjs data, a pseudopop will have n=4
%%%%%               "cellID" fields -- one per day (two barney, two thor).
%%%%%               Those tell you for a given pseudopop, a given day,
%%%%%                   which cells responses made it into the pseudopop
%%%%% Sampling with replacement
%%% Each trial is nTimepoints x nCells array

clear all; close all;
workdir = '/Users/geena/Dropbox/PhD/SUAinfo/PSIDResults2023/';
subj2analyze = {'Barney','Thor'};
regionList = {'S1','M1','PMv','M3','All'};

collapseSubjectsFlag = true;
nIter = 25; % number of pseudopopulations /cell groupings to make
nUnitsPerRegion = 50; % number of cells per region
useAllUnitsFlag = false;
sampleNeuronsWithReplacement = true;

bhv2use = 'ThrOnly'; % LSOnly, LSThrOnly, allBhv, etc.
saveFlag = true;

% set up output directories
if saveFlag && ~collapseSubjectsFlag && ~useAllUnitsFlag 
    outDir = ['/Users/geena/Dropbox/PhD/SUAinfo/PSIDResults2023/Pseudopopulations/' bhv2use '/balanced/separateSubjs/N_' num2str(nUnitsPerRegion) 'cells/'];
elseif saveFlag && ~collapseSubjectsFlag && useAllUnitsFlag 
    outDir = ['/Users/geena/Dropbox/PhD/SUAinfo/PSIDResults2023/Pseudopopulations/' bhv2use '/unbalanced/separateSubjs/N_' num2str(nUnitsPerRegion) 'cells/'];
elseif saveFlag && collapseSubjectsFlag && ~useAllUnitsFlag 
    outDir = ['/Users/geena/Dropbox/PhD/SUAinfo/PSIDResults2023/Pseudopopulations/' bhv2use '/balanced/combinedSubjs/N_' num2str(nUnitsPerRegion) 'cells/'];
elseif saveFlag && collapseSubjectsFlag && useAllUnitsFlag 
    outDir = ['/Users/geena/Dropbox/PhD/SUAinfo/PSIDResults2023/Pseudopopulations/' bhv2use '/unbalanced/combinedSubjs/N_' num2str(nUnitsPerRegion) 'cells/'];
end

% To combine subjects
if collapseSubjectsFlag
    subj2analyze = {'combined'};
end

for iter = 1: nIter
    for ss = 1:length(subj2analyze)

        thisSubj = subj2analyze{ss};

        for rr = 1:length(regionList)

            allTrialsNeural = [];
            allTrialsBhv = [];
            allTrialsBhvLabels = [];
            allTrialsCellIDs = []; 

            allTrialData = load([workdir regionList{rr} '_' bhv2use 'Trials.mat']);
            
            % if collapsing subjects, take all dates, otherwise loop thru
            if collapseSubjectsFlag
                dates = unique(allTrialData.allLabels(:,4),'stable');
            else
                dates = unique(allTrialData.allLabels(strcmpi(allTrialData.allLabels(:,3),thisSubj),4),'stable');
            end

            % for each recording date (ie., one "set" of neurons), grab
            % neural & bhv trials 
            for dd = 1:length(dates)
                
                theseTrialsPseudopop = [];

                % grab trials from that date
                dateIdx = strcmpi(allTrialData.allLabels(:,4), dates{dd});
                theseNeuralTrials = allTrialData.allNeuralTrials(dateIdx);
                theseBhvTrials = allTrialData.allBhvTrials(dateIdx);
                theseBhvTrialLabels = allTrialData.allBhvLabels(dateIdx);

                % grab cellIDs from that date
                if collapseSubjectsFlag
                    if dd == 1
                        cellIDs = cellstr(allTrialData.cellListOut{1,1}{1,1});
                    elseif dd == 2
                        cellIDs = cellstr(allTrialData.cellListOut{1,1}{1,2});
                    elseif dd == 3
                        cellIDs = cellstr(allTrialData.cellListOut{1,2}{1,1});
                    elseif dd == 4
                        cellIDs = cellstr(allTrialData.cellListOut{1,2}{1,2});
                    end
                else
                    cellIDs = cellstr(allTrialData.cellListOut{1,ss}{1,dd});
                end

                if ~useAllUnitsFlag % if selecting subset of cells
                    if sampleNeuronsWithReplacement %  with replacement
                        randIdx = randi(length(cellIDs),nUnitsPerRegion,1);
                    else %  without replacement
                        randIdx = randperm(length(cellIDs),nUnitsPerRegion);
                    end
                else % use all cells
                    disp('not set up to do this yet!')
                    if sampleNeuronsWithReplacement
                        randIdx = randi(length(cellIDs),length(cellIDs),1);
                    else
                        randIdx = randperm(length(cellIDs),length(cellIDs));
                    end
                end

                randUnits = cellIDs(randIdx);

                % pull out trials (neuraldata & bhv) & stash
                for tt = 1:length(theseNeuralTrials)
                    theseTrialsPseudopop{tt,1} = theseNeuralTrials{tt,1}(:,randIdx);

                end

                allTrialsNeural = [allTrialsNeural ;theseTrialsPseudopop];
                allTrialsBhv = [allTrialsBhv; theseBhvTrials];
                allTrialsBhvLabels = [allTrialsBhvLabels; theseBhvTrialLabels];
                allTrialsCellIDs = [allTrialsCellIDs; {randUnits}];

            end % end dates

            PSIDpseudopop.(thisSubj).(regionList{rr}).neural = allTrialsNeural;
            PSIDpseudopop.(thisSubj).(regionList{rr}).cellIDs = allTrialsCellIDs;

        end % end regions

        PSIDpseudopop.(thisSubj).bhv = allTrialsBhv;
        PSIDpseudopop.(thisSubj).bhvLabels = allTrialsBhvLabels;

    end % end subjects

    % grab a couple variables 
    taxis = -allTrialData.tmin:allTrialData.analysisWin:allTrialData.tmax;
    bhvPCs = allTrialData.bhvPCs; bhvs2train = allTrialData.bhvs2train; bhvs2test = allTrialData.bhvs2test;
    markerPCA = allTrialData.markerPCA; markers2use = allTrialData.markers2use;

    if ~exist(outDir, 'dir')
        mkdir(outDir)
    end

    % save it 
    save([outDir '/pp_N' num2str(iter) '.mat'],...
        'PSIDpseudopop','collapseSubjectsFlag','iter','nUnitsPerRegion','regionList','sampleNeuronsWithReplacement','useAllUnitsFlag',...
        'taxis','bhvPCs','bhvs2test','bhvs2train','markerPCA','markers2use', '-v7.3');
    
    disp(['saving N= ' num2str(iter) ' ...']);

end % end iterations
