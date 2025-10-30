function resultsOut = runCorePSIDAnalysis(data,thisSubj, regionList, outcomeMetric)

% assumes you've run prepData2023.m --> makePseudopopPSID2023.m
% latentStatePredictions --  behaviorally-relevant (first n1 columns) or irrelevant states (last n2 columns). 
%
checkPlotRawDataFlag = false; 
plotN1TuningResultsFlag = false; plotNXTuningResultsFlag = false;
zScoreNeural = true;  zScoreBhv = true;

% Set the number of outer and inner folds
outerFolds = 4;
innerFolds = 8;

% Set the hyperparameter values to evaluate
statesBhvRelv = [2 3 4 8 10 12 20]; %  N1
additionalStates = [2 4 8 10]; % ** this is NOT Nx (Nx = statesBhvRelv + additional) **
horizonValues = [5 10 15  ];  % max N1 <= nz(bhv) * i ; max Nx <= ny(neural) * i
nz = 5;                       % bhvDims       
%taxis = -0.5:0.05:1.5;

% define behavior (kinematic & categorical labels)
bhvPCs2Predict = 1:nz;
nBhvVars = length(bhvPCs2Predict);
bhvs2test = data.bhvs2test;
allBhvTrialLabels = cell2mat(data.PSIDpseudopop.(thisSubj).bhvLabels);
allBhvTrials = data.PSIDpseudopop.(thisSubj).bhv; fun = @(x) x(:,bhvPCs2Predict);
allBhvTrials = cellfun(fun,allBhvTrials,'UniformOutput',false);

for rr = 1:length(regionList)

    % define neural
    taxis = data.taxis;
    allNeuralTrials = data.PSIDpseudopop.(thisSubj).(regionList{rr}).neural;

    % Outer cross-validation loop
    outerCV = cvpartition(allBhvTrialLabels,'KFold',outerFolds,'Stratify',true);

    % Initialize variables to store results
    outerKinematicAccuracy = zeros(outerFolds, nBhvVars);

    for outerFold = 1:outerFolds

        % Extract the training and testing indices for the current outer fold
        outerTrainIdx = training(outerCV, outerFold);
        outerTestIdx = test(outerCV, outerFold);

        % Split the data into training and testing sets for the outer fold
        outerNeuralTrain = allNeuralTrials(outerTrainIdx, :);
        outerBhvTrain = allBhvTrials(outerTrainIdx);
        outerBhvTrainLabels = allBhvTrialLabels(outerTrainIdx);

        outerNeuralTest = allNeuralTrials(outerTestIdx, :);
        outerBhvTest = allBhvTrials(outerTestIdx);
        outerBhvTestLabels = allBhvTrialLabels(outerTestIdx);

        % Inner cross-validation loop for hyperparameter tuning
        innerCV = cvpartition(outerBhvTrainLabels, 'KFold', innerFolds,'Stratify',true);
        nCombos = length(statesBhvRelv)*length(horizonValues);
        combo = 1;

        for statesBhvRelvIdx = 1:length(statesBhvRelv)

            for horizonIdx = 1:length(horizonValues)

                % Set current hyperparameters
                currentHyperparam.statesBhvRelv = statesBhvRelv(statesBhvRelvIdx);
                currentHyperparam.horizon = horizonValues(horizonIdx);

                % Initialize to store validation performance across inner folds
                innerKinematicAccuracy = zeros(innerFolds, nBhvVars);

                for innerFold = 1:innerFolds

                    % Extract the training and validation indices for the current inner fold
                    innerTrainIdx = training(innerCV, innerFold);
                    innerValIdx = test(innerCV, innerFold);

                    % Split the data into training and validation sets for the inner fold
                    innerNeuralTrain = outerNeuralTrain(innerTrainIdx, :);
                    innerBhvTrain = outerBhvTrain(innerTrainIdx);
                    innerBhvTrainLabels = outerBhvTrainLabels(innerTrainIdx);

                    innerNeuralVal = outerNeuralTrain(innerValIdx, :);
                    innerBhvVal = outerBhvTrain(innerValIdx); innerBhvValLabels = outerBhvTrainLabels(innerValIdx);

                    if checkPlotRawDataFlag
                        figure;
                        for i = 1:length(innerNeuralVal)
                            subplot 211
                            imagesc(taxis, [], innerNeuralVal{i,:}'); colormap(jet); colorbar; hold on;
                            xlabel('time, s'); ylabel('cells'); title('Neural Activity')
                            subplot 212
                            plot(taxis, innerBhvVal{i,1}(:,1),'linew',6)
                            xlabel('time, s'); ylabel('face PC1, au'); title('Behavior')

                            pause; clf;
                        end
                    end

                    % train the model, using one combo of hyperparms
                    %  ONLY using first stage, ie., n1=nX
                    idSys = PSID(innerNeuralTrain, innerBhvTrain, currentHyperparam.statesBhvRelv, currentHyperparam.statesBhvRelv, currentHyperparam.horizon,[],[],[],...
                        'zscore_Y',zScoreNeural,'zscore_Z',zScoreBhv);

                    % Predict behavior using the learned model
                    %  ONLY using first stage, ie., n1=nX
                    [innerBhvPred, ~, ~] = PSIDPredict(idSys, innerNeuralVal);
                    innerBhvPred = cell2mat(innerBhvPred); innerBhvVal = cell2mat(innerBhvVal);

                    % Evaluate the classifier on the validation set; prediction is nFolds x nPredictors  (calc btwn trials)
                    innerKinematicAccuracy(innerFold,:) = evalPrediction(innerBhvVal, innerBhvPred, outcomeMetric);


                end % end innerFolds

                % Compute the mean accuracy across inner folds
                optimN1(combo).meanInnerKinematicAccuracy = mean(innerKinematicAccuracy,'omitnan')';
                optimN1(combo).statesBhvRelv = currentHyperparam.statesBhvRelv;
                optimN1(combo).horz = currentHyperparam.horizon;

                combo = combo + 1;
            end % end 'horizon' gridsearch
        end % end 'statesBhvRelv' gridsearch

        % In inner CV of N1, pick lowest N1 within 1 SEM of peak decoding
        meanInnerPerformanceAll = mean([optimN1(:).meanInnerKinematicAccuracy],1, 'omitnan'); % averaged across nBhvsPredicted & nFolds
        meanInnerPerformancePerBhv = [optimN1(:).meanInnerKinematicAccuracy];
        searchStatesBhvRelv = [optimN1(:).statesBhvRelv]; searchHorz = [optimN1(:).horz];

        %  optim N1
        [peakPerf, ~] = max(meanInnerPerformanceAll);
        semAcc = std(meanInnerPerformanceAll,'omitnan') / sqrt(length(meanInnerPerformanceAll));
        rangePeakPerf = [(peakPerf - semAcc), (peakPerf + semAcc)];
        if diff(rangePeakPerf) >= 0.001
            rangePeakPerf = round(rangePeakPerf,3);
        end
        optBhvRelevantStates = min(searchStatesBhvRelv(meanInnerPerformanceAll >= rangePeakPerf(1) & meanInnerPerformanceAll <= rangePeakPerf(2)));

        %  optim horiz
        [optHorzIdx, ~] = max(meanInnerPerformanceAll(searchStatesBhvRelv == optBhvRelevantStates));
        optHorz = searchHorz(find(meanInnerPerformanceAll == optHorzIdx));

        % plot the hyperparm tuning results for N1/horz
        if plotN1TuningResultsFlag
            plotN1TuningResults(meanInnerPerformanceAll, meanInnerPerformancePerBhv, peakPerf, rangePeakPerf,searchStatesBhvRelv, optBhvRelevantStates, horizonValues, optHorz, bhvPCs2Predict)
        end

        %%%
        %%%
        %%% NOW THAT you have optim N1/horiz,
        %%% find optim NX w/ those values (ie., re-enter innercross val)
        statesTotal = optBhvRelevantStates + additionalStates;

        for statesTotalIdx = 1:length(statesTotal)

            % Set current hyperparameters
            currentHyperparam.statesTotal = statesTotal(statesTotalIdx);
            innerNeuralSelfPredAccuracy = zeros(innerFolds, size(outerNeuralTest{1},2));

            for innerFold = 1:innerFolds

                % Extract the training and validation indices for the current inner fold
                innerTrainIdx = training(innerCV, innerFold);
                innerValIdx = test(innerCV, innerFold);

                % Split the data into training and validation sets for the inner fold
                innerNeuralTrain = outerNeuralTrain(innerTrainIdx, :);
                innerBhvTrain = outerBhvTrain(innerTrainIdx); innerBhvTrainLabels = outerBhvTrainLabels(innerTrainIdx);

                innerNeuralVal = outerNeuralTrain(innerValIdx, :);
                innerBhvVal = outerBhvTrain(innerValIdx); innerBhvValLabels = outerBhvTrainLabels(innerValIdx);

                % train the model, using optim N1/horiz & ONE Nx
                idSys2 = PSID(innerNeuralTrain, innerBhvTrain, currentHyperparam.statesTotal, optBhvRelevantStates, optHorz,[],[],[],...
                    'zscore_Y',zScoreNeural,'zscore_Z',zScoreBhv);

                % Predict NEURAL activity this time (self-pred), using both
                % stages N1 & NX
                [~, innerNeuralSelfPred, ~] = PSIDPredict(idSys2, innerNeuralVal);
                innerNeuralSelfPred = cell2mat(innerNeuralSelfPred); innerNeuralVal = cell2mat(innerNeuralVal);

                % Evaluate the classifier on the validation set; prediction is nFolds x nPredictors  (calc btwn trials)
                innerNeuralSelfPredAccuracy(innerFold,:) = evalPrediction(innerNeuralVal,innerNeuralSelfPred,outcomeMetric);

            end % end innerFolds (again)

            optimNx(statesTotalIdx).meanInnerSelfPred = mean(innerNeuralSelfPredAccuracy,'omitnan')';
            optimNx(statesTotalIdx).statesTotal = currentHyperparam.statesTotal;

        end % end gridsearch, statesTotal

        % In inner CV pick lowest NX within 1 SEM of peak neural SELF prediction
        meanInnerSelfPredAll = mean([optimNx(:).meanInnerSelfPred],1); % averaged across nFeatures & nFolds
        searchStatesTotal = [optimNx(:).statesTotal];

        %  optim NX
        [peakSelf, ~] = max(meanInnerSelfPredAll);
        semSelf = std(meanInnerSelfPredAll) / sqrt(length(meanInnerSelfPredAll));
        rangePeakSelf = [(peakSelf - semSelf), (peakSelf + semSelf)];
        if diff(rangePeakSelf) >= 0.001
            rangePeakSelf = round(rangePeakSelf,3);
        end
        optTotalStates = min(searchStatesTotal(meanInnerSelfPredAll >= rangePeakSelf(1) & meanInnerSelfPredAll <= rangePeakSelf(2)));

        % plot the hyperparm tuning results for NX
        if plotNXTuningResultsFlag
            plotNXTuningResults(meanInnerSelfPredAll, statesTotal, peakSelf,rangePeakSelf, optimNx, optBhvRelevantStates, optHorz)
        end

        %%%%% NOW you have all three hyperparms optimized
        disp(['optBhvRelevantStates = ' num2str(optBhvRelevantStates)]);
        disp(['optTotalStates = ' num2str(optTotalStates)]);
        disp(['optHorz = ' num2str(optHorz)])
        disp('\\\\\\\\\\\\')

        %%%%%%%
        % Re-Train the final classifier on the entire outer training set (using
        % 2 stages)
        idSysFinal_2Stage = PSID(outerNeuralTrain, outerBhvTrain, optTotalStates, optBhvRelevantStates, optHorz,[],[],[],...
            'zscore_Y',zScoreNeural,'zscore_Z',zScoreBhv);

        % Evaluate the final, 2-stage classifier on the outer test set
        % R2/CC is nFolds x nPredictors  (calc btwn trials)
        [outerBhvPred, outerNeuralSelfPred_2stage, outerLatentPred] = PSIDPredict(idSysFinal_2Stage, outerNeuralTest);
        outerBhvPredC = cell2mat(outerBhvPred); outerBhvTestC = cell2mat(outerBhvTest); outerNeuralTestC = cell2mat(outerNeuralTest);
        outerNeuralSelfPred_2stageC = cell2mat(outerNeuralSelfPred_2stage);outerLatentPredC = cell2mat(outerLatentPred);

        outerKinematicAccuracy(outerFold,:) = evalPrediction(outerBhvTestC, outerBhvPredC, outcomeMetric);
        outerNeuralSelfPred_2stageAcc(outerFold, :) = evalPrediction(outerNeuralTestC,outerNeuralSelfPred_2stageC,outcomeMetric);

        % Re-Train a separate, final classifier on the entire outer training set (using
        % only FIRST stage, ie. NX = N1)
        idSysFinal_1Stage = PSID(outerNeuralTrain, outerBhvTrain, optBhvRelevantStates, optBhvRelevantStates, optHorz,[],[],[],...
            'zscore_Y',zScoreNeural,'zscore_Z',zScoreBhv);

        % Evaluate the final, 1-stage classifier on the outer test set
        % R2/CC is nFolds x nPredictors  (calc btwn trials)
        [~, outerNeuralSelfPred_1stage, ~] = PSIDPredict(idSysFinal_1Stage, outerNeuralTest);
        outerNeuralSelfPred_1stageC = cell2mat(outerNeuralSelfPred_1stage);
        outerNeuralSelfPred_1stageAcc(outerFold, :) = evalPrediction(outerNeuralTestC,outerNeuralSelfPred_1stageC,outcomeMetric);

    end

    [stdKinematicAcc, meanKinematicAcc] = std(outerKinematicAccuracy);
    [stdSelfPred_2stageAcc, meanSelfPred_2stageAcc] = std(outerNeuralSelfPred_2stageAcc,[],'all');
    [stdSelfPred_1stageAcc, meanSelfPred_1stageAcc] = std(outerNeuralSelfPred_1stageAcc,[],'all');

    % Display the average kinematic & neural pred accuracy across outer folds
    disp(['REGION = ' regionList{rr}])
    disp(['Average Kinematic PC1, ' outcomeMetric ' : ' num2str(meanKinematicAcc(1)) ' +/- ' num2str(stdKinematicAcc(1))]);
    disp(['Average Kinematic PC2, ' outcomeMetric ' : ' num2str(meanKinematicAcc(2)) ' +/- ' num2str(stdKinematicAcc(2))]);
    disp(['Average Kinematic PC3, ' outcomeMetric ' : ' num2str(meanKinematicAcc(3)) ' +/- ' num2str(stdKinematicAcc(3))]);

    disp(['Average Neural, Full Model ' outcomeMetric ' : ' num2str(meanSelfPred_2stageAcc) ' +/- ' num2str(stdSelfPred_2stageAcc)]);
    disp(['Average Neural First-stage only ' outcomeMetric ' : ' num2str(meanSelfPred_1stageAcc) ' +/- ' num2str(stdSelfPred_1stageAcc)]);
    disp(['Kinematic : AllOther Info Ratio = ' num2str((meanSelfPred_1stageAcc / meanSelfPred_2stageAcc))])

    % stash results
    resultsOut.(thisSubj).(regionList{rr}).kinematicDecoding = outerKinematicAccuracy;
    resultsOut.(thisSubj).(regionList{rr}).outcomeMetric = outcomeMetric;
    resultsOut.(thisSubj).(regionList{rr}).neuralSelfPred2stage = outerNeuralSelfPred_2stageAcc;
    resultsOut.(thisSubj).(regionList{rr}).neuralSelfPred1stage = outerNeuralSelfPred_1stageAcc;
    resultsOut.(thisSubj).(regionList{rr}).N1 = optBhvRelevantStates;
    resultsOut.(thisSubj).(regionList{rr}).NX = optTotalStates;
    resultsOut.(thisSubj).(regionList{rr}).horz = optHorz;
    resultsOut.(thisSubj).(regionList{rr}).bhvTest = outerBhvTest;
    resultsOut.(thisSubj).(regionList{rr}).bhvPred= outerBhvPred;
    resultsOut.(thisSubj).(regionList{rr}).bhvLabels = outerBhvTestLabels;
    resultsOut.(thisSubj).(regionList{rr}).latent = outerLatentPred;
    resultsOut.(thisSubj).(regionList{rr}).finalModel = idSysFinal_2Stage;

end % end regions


end

function [] = plotN1TuningResults(meanInnerPerformanceAll, meanInnerPerfomancePerBhv, peakPerf, rangePeakPerf,searchStatesBhvRelv, optBhvRelevantStates, horizonValues, optHorz, bhvPCs2Predict)

nCombos = length(unique(searchStatesBhvRelv)) * length(horizonValues);

legendStrings = "i = " + string(horizonValues);
x = [searchStatesBhvRelv(1) searchStatesBhvRelv(1) searchStatesBhvRelv(end) searchStatesBhvRelv(end)];
y = [rangePeakPerf(end) rangePeakPerf(1) rangePeakPerf(1) rangePeakPerf(end)];

figure;
plot(searchStatesBhvRelv(1:length(horizonValues):nCombos),meanInnerPerformanceAll(1,1:length(horizonValues):nCombos),'o','linew',4); hold on
plot(searchStatesBhvRelv(2:length(horizonValues):nCombos),meanInnerPerformanceAll(1,2:length(horizonValues):nCombos),'o','linew',4);
plot(searchStatesBhvRelv(3:length(horizonValues):nCombos),meanInnerPerformanceAll(1,3:length(horizonValues):nCombos),'o','linew',4);
fill(x,y,'b','FaceAlpha',0.2); yline(peakPerf,'b--','linew',4); hold off

xlabel('# rel states'); ylabel('Decoding Performance');
title(['Avg All BhvVars; Opt BhvRelv States = ' num2str(optBhvRelevantStates) ' Horiz = ' num2str(optHorz)]); legend(legendStrings)

figure;
for bb = 1:length(bhvPCs2Predict) % for each bhvVar
    subplot(1,length(bhvPCs2Predict),bb)
    plot(searchStatesBhvRelv(1:length(horizonValues):nCombos),meanInnerPerfomancePerBhv(bb,1:length(horizonValues):nCombos),'o','linew',4);hold on
    plot(searchStatesBhvRelv(2:length(horizonValues):nCombos),meanInnerPerfomancePerBhv(bb,2:length(horizonValues):nCombos),'o','linew',4);
    plot(searchStatesBhvRelv(3:length(horizonValues):nCombos),meanInnerPerfomancePerBhv(bb,3:length(horizonValues):nCombos),'o','linew',4); hold off
    xlabel('# rel states'); ylabel('performance'); title(['Bhv Var = ' num2str(bb)]); legend(legendStrings)
    sgtitle(['Per BhvVars; Opt BhvRelv States = ' num2str(optBhvRelevantStates) ' Horiz = ' num2str(optHorz)])
end


end

function [] = plotNXTuningResults(meanInnerSelfPredAll, statesTotal, peakSelf, rangePeakSelf, optimNx, optBhvRelevantStates, optHorz)

x = [statesTotal(1) statesTotal(1) statesTotal(end) statesTotal(end)];
y = [rangePeakSelf(end) rangePeakSelf(1) rangePeakSelf(1) rangePeakSelf(end)];

figure;
plot(statesTotal, meanInnerSelfPredAll,'o','linew',6); hold on
fill(x,y,'b','FaceAlpha',0.2); yline(peakSelf,'b--','linew',4); hold off
xlim([statesTotal(1) statesTotal(end)])
xlabel('# total states'); ylabel('Decoding Performance');
title(['Avg NeuralSelfPred; Opt BhvRelv States = ' num2str(optBhvRelevantStates) ' Horiz = ' num2str(optHorz)]);



end