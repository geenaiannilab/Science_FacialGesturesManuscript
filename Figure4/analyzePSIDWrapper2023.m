clear all; close all; 
set(0,'defaultAxesFontSize', 16); % bc im blind 


% define the dataset you are running
thisSubj = 'combined';
pseduoPop = 'LSonly/balanced'; % vs 'unbalanced'
subj = 'combinedSubjs'; % vs 'separateSubjs' v 'combinedSubjs'
nullModelFlag = true;  nullModelType = 'shuffleBhvTimepoints';
nNullShuffles = 50;

nPseudopops = 25; % pseduopopulations
popSize = '50'; % nCells
outcomeMetric = 'CC'; % vs R2 vs RSME etc
saveFlag = true; 

dataDir = (['/Users/geena/Dropbox/PhD/SUAinfo/PSIDResults2023/Pseudopopulations/' pseduoPop '/' subj '/N_' popSize 'cells/']);
outDir = ['/Users/geena/Dropbox/PhD/SUAinfo/PSIDResults2023/Pseudopopulations/' pseduoPop '/' subj '/N_' popSize 'cells/results'];

for pp = 4:nPseudopops

   % load one pseudopop data
   trueData = load([dataDir 'pp_N' num2str(pp) '.mat']);
   regionList = trueData.regionList;

   tic 

   % run true results
    trueResults = runCorePSIDAnalysis(trueData, thisSubj, regionList, outcomeMetric);

    % run null shuffles
    for nn = 1:nNullShuffles
        shuffledData = makeDataShuffled(trueData, nullModelType, regionList);
    	shuffledResult(nn) = runCorePSIDAnalysis(shuffledData, thisSubj, regionList, outcomeMetric);
        disp(['SHUFFLE # ' num2str(nn)])
    end 

    for rr  = 1:length(regionList)
        
        % pull out decoding mean/std, neural pred mean/std across shuffles
        meanDecodingNull = arrayfun(@(x) x.combined.(regionList{rr}).kinematicDecoding, shuffledResult,'UniformOutput',false);
        meanDecodingNull = cellfun(@mean , meanDecodingNull,'UniformOutput',false);
        [stdCompareDecodingNull, meanCompareDecodingNull] = std(cell2mat(meanDecodingNull'),[],1); % nBhvVars long; averaged over nOuterFolds per shuffle, and then over shuffles

        meanSelfPred1Null = arrayfun(@(x) x.combined.(regionList{rr}).neuralSelfPred1stage, shuffledResult,'UniformOutput',false);
        meanSelfPred1Null = cellfun(@mean, meanSelfPred1Null, 'UniformOutput',false);
        [stdCompareSelfPred1Null, meanCompareSelfPred1Null] = std(cell2mat(meanSelfPred1Null'),[],1); % nCells long; averaged over nOuterFolds per shuffle, and then over shuffles

        meanSelfPredFullNull = arrayfun(@(x) x.combined.(regionList{rr}).neuralSelfPred2stage, shuffledResult,'UniformOutput',false);
        meanSelfPredFullNull = cellfun(@mean, meanSelfPredFullNull, 'UniformOutput',false);
        [stdCompareSelfPredFullNull, meanCompareSelfPredFullNull] = std(cell2mat(meanSelfPredFullNull'),[],1); % nCells long; averaged over nOuterFolds per shuffle, and then over shuffles

        % stash null results
        shuffledResultsAll.(thisSubj).(regionList{rr}).kinematicDecoding(:,pp) = meanCompareDecodingNull'; 
        shuffledResultsAll.(thisSubj).(regionList{rr}).kinematicDecodingStd(:,pp) = stdCompareDecodingNull';
        shuffledResultsAll.(thisSubj).(regionList{rr}).outcomeMetric = outcomeMetric;
        shuffledResultsAll.(thisSubj).(regionList{rr}).neuralSelfPred2stage(:,pp) =  meanCompareSelfPredFullNull';
        shuffledResultsAll.(thisSubj).(regionList{rr}).neuralSelfPred2stageStd(:,pp) = stdCompareSelfPredFullNull';
        shuffledResultsAll.(thisSubj).(regionList{rr}).neuralSelfPred1stage(:,pp) = meanCompareSelfPred1Null';
        shuffledResultsAll.(thisSubj).(regionList{rr}).neuralSelfPred1stageStd(:,pp) = stdCompareSelfPred1Null';
        shuffledResultsAll.(thisSubj).(regionList{rr}).N1(pp) = mean(arrayfun(@(x) x.combined.(regionList{rr}).N1, shuffledResult));
        shuffledResultsAll.(thisSubj).(regionList{rr}).NX(pp) = mean(arrayfun(@(x) x.combined.(regionList{rr}).NX, shuffledResult));
        shuffledResultsAll.(thisSubj).(regionList{rr}).horz(pp) = mean(arrayfun(@(x) x.combined.(regionList{rr}).horz, shuffledResult));

        % stash true results
%         trueResultsAll.(thisSubj).(regionList{rr}).kinematicDecoding(:,:,pp) = trueResults.(thisSubj).(regionList{rr}).kinematicDecoding;
%         trueResultsAll.(thisSubj).(regionList{rr}).outcomeMetric = trueResults.(thisSubj).(regionList{rr}).outcomeMetric;
%         trueResultsAll.(thisSubj).(regionList{rr}).neuralSelfPred2stage(:,:,pp) = trueResults.(thisSubj).(regionList{rr}).neuralSelfPred2stage;
%         trueResultsAll.(thisSubj).(regionList{rr}).neuralSelfPred1stage(:,:,pp) = trueResults.(thisSubj).(regionList{rr}).neuralSelfPred1stage;
%         trueResultsAll.(thisSubj).(regionList{rr}).N1(pp) = trueResults.(thisSubj).(regionList{rr}).N1;
%         trueResultsAll.(thisSubj).(regionList{rr}).NX(pp) = trueResults.(thisSubj).(regionList{rr}).NX;
%         trueResultsAll.(thisSubj).(regionList{rr}).horz(pp) = trueResults.(thisSubj).(regionList{rr}).horz;
%         trueResultsAll.(thisSubj).(regionList{rr}).bhvTest{pp}= trueResults.(thisSubj).(regionList{rr}).bhvTest;
%         trueResultsAll.(thisSubj).(regionList{rr}).bhvPred{pp}= trueResults.(thisSubj).(regionList{rr}).bhvPred;
%         trueResultsAll.(thisSubj).(regionList{rr}).bhvLabels{pp} = trueResults.(thisSubj).(regionList{rr}).bhvLabels;
%         trueResultsAll.(thisSubj).(regionList{rr}).latent{pp} = trueResults.(thisSubj).(regionList{rr}).latent;
%         trueResultsAll.(thisSubj).(regionList{rr}).finalModel(pp) = trueResults.(thisSubj).(regionList{rr}).finalModel;
    end

    if saveFlag
        if ~exist(outDir, 'dir')
            mkdir(outDir)
        end
       % save([outDir '/decodingResults.mat' ],'trueResultsAll','-v7.3');
        save([outDir '/decodingNullResults_pp_N' num2str(num2str(pp))  '.mat' ],'shuffledResultsAll','nullModelFlag', 'nullModelType', 'nNullShuffles','-v7.3');
    end
    clear shuffledResult; clear shuffledData;
    toc
end




% shuffle Data function
function data = makeDataShuffled(data, shuffleMethod, regionList)

% shuffle method

if strcmpi(shuffleMethod,'shuffleNeuralTimepoints')

    for rr = 1:length(regionList)

        nTrials = length(data.PSIDpseudopop.combined.(regionList{rr}).neural);

        for tt = 1:nTrials

            nCells = size(data.PSIDpseudopop.combined.(regionList{rr}).neural{tt},2);
            nTimepoints = size(data.PSIDpseudopop.combined.(regionList{rr}).neural{tt},1);

            for cc = 1:nCells

                shuffledTimepoints = randperm(nTimepoints);
                data.PSIDpseudopop.combined.(regionList{rr}).neural{tt}(:,cc) = data.PSIDpseudopop.combined.(regionList{rr}).neural{tt}(shuffledTimepoints,cc);

            end % end nCells

        end % end trials 

    end % end regions

elseif strcmpi(shuffleMethod,'shuffleBhvTimepoints')
   
    nTrials = length(data.PSIDpseudopop.combined.bhv);

    for tt = 1:nTrials
        nBhvVar = size(data.PSIDpseudopop.combined.bhv{tt},2);
        nTimepoints = size(data.PSIDpseudopop.combined.bhv{tt},1);

        for bb = 1:nBhvVar
            shuffledTimepoints = randperm(nTimepoints);
            data.PSIDpseudopop.combined.bhv{tt}(:,bb) = data.PSIDpseudopop.combined.bhv{tt}(shuffledTimepoints,bb);

        end % end bhvVar

    end % endTrials 

end % end shuffleFlagType

end % end function




