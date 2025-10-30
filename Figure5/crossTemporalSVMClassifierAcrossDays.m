%%% written GI 211215ish
%%%
%%% Implements SVM classifier (fitceoc) to classify individual trials of facial
%%% expression behavior, on basis of single-trial PCA neural data.
%%% Current implementation uses 
%%% Crossvalidated by k=20 folds of 80/20 splits. 
%%% Saves validation scores & confusion matrices  



function [ meanValidationScores,SEMValidationScores, meanBetaWeights ] = crossTemporalSVMClassifierAcrossDays(varargin)

%% Yes, you're getting an input parser
p = inputParser;

% Optional params
p.addParameter('trainNeural','trainNeural');
p.addParameter('testNeural','testNeural');
p.addParameter('subject','Barney',@ischar);
p.addParameter('array','M1',@ischar);
p.addParameter('nFolds',8,@isscalar);
p.addParameter('trialLabels','trialLabels');

p.addParameter('denoiseFlag',false,@islogical); % PCA-based denoising 
p.addParameter('zScoreFlag', 1,@isscalar); 
p.addParameter('threshlowFRFlag', 1,@isscalar); % 0/1; whether or not to remove low FR neurons 
p.addParameter('threshlowFR', 0.01,@isscalar); %  sp/s threshold for removing low FR neurons 

p.addParameter('bhvs2classify', [1 2 4]);
p.addParameter('socialVNotDecoder',false)

p.addParameter('minRestFlag', true);
p.addParameter('minRest', 1,@isscalar);

p.addParameter('trainTmin', 0,@isscalar);
p.addParameter('trainTmax',1,@isscalar);
p.addParameter('testTmin', 0,@isscalar);
p.addParameter('testTmax',1,@isscalar);
p.addParameter('win', 0.02,@isscalar); % in sec

p.addParameter('plotFlag', false); 
p.addParameter('saveFig', false);
p.addParameter('verboseFlag',true);

% parse
p.parse(varargin{:})

% because we're paranoid
rng('shuffle');
close all; 

subject = p.Results.subject;
trainNeural = p.Results.trainNeural;
testNeural = p.Results.testNeural;
denoiseFlag = p.Results.denoiseFlag;
bhvs2classify = p.Results.bhvs2classify;
socialVNotDecoder = p.Results.socialVNotDecoder;
trialLabels = p.Results.trialLabels; 

trainTmin = p.Results.trainTmin; 
trainTmax = p.Results.trainTmax; 
testTmin = p.Results.testTmin; 
testTmax = p.Results.testTmax; 

win = p.Results.win; % in sec
trainAx = -trainTmin:win:trainTmax;
testAx = -testTmin:win:testTmax;

plotFlag = p.Results.plotFlag; 
nFolds = p.Results.nFolds;
zScoreFlag = p.Results.zScoreFlag;  % zscoring matters ALOT here bc FRs are so different 
threshlowFRFlag = p.Results.threshlowFRFlag; % 0/1; whether or not to remove low FR neurons 
threshlowFR = p.Results.threshlowFR;
minRestFlag = p.Results.minRestFlag;
minRest = p.Results.minRest;

% remove low FR neurons & normalize 
tmp = [trainNeural; testNeural];
meanFRs = mean(tmp,1) ./ abs(trainTmax - trainTmin);

% remove low FR neurons
if threshlowFRFlag
    tmp = tmp(:,meanFRs > threshlowFR);
end

% take top N PCs accounting for >98% variance, & recompose data 
if denoiseFlag
    mu = mean(tmp);
    [eigenvectors, scores,latent] = pca(tmp);
    varEx = cumsum(latent/sum(latent));
    nComp = find(varEx >= 0.98,1, 'first'); % take top N PCs >= 95% variance
    Xhat = scores(:,1:nComp) * eigenvectors(:,1:nComp)'; 
    tmp = bsxfun(@plus, Xhat, mu);
end

% zScore 
if zScoreFlag
    tmpZ = zscore(tmp);
    trainNeural = tmpZ(1:size(trainNeural,1),:);
    testNeural = tmpZ(size(trainNeural,1)+1:end,:);

else
    trainNeural = tmp(1:size(trainNeural,1),:);
    testNeural = tmp(size(trainNeural,1)+1:end,:);
end

% data entry is nTrials x nCells , 
%    where each entry is # of spikes within a trial

%% split the data 
cv = cvpartition(trialLabels,'KFold',nFolds,'Stratify',true);

for fold = 1:nFolds

    % define training & validation indices based on cvpartitions 
    trainingNeural = trainNeural(training(cv,fold),:);
    trainingBhv = trialLabels(training(cv,fold));
    
    validationNeural = testNeural(test(cv,fold),:);
    validationBhv = trialLabels(test(cv,fold));
  

    % fit SVM on the training data
    t = templateSVM('Standardize',true,'SaveSupportVectors','on');
    SVMmodel = fitcecoc(trainingNeural,trainingBhv,'Learners',t);
    
    % predict for the validation data 
    % predBhvScore is an testTrials-by-nClasses matrix of soft scores; 
    [predictedBhv,~] = predict(SVMmodel,validationNeural);
    
    % compute accuracy 
    validationAcc = sum(validationBhv == predictedBhv) / numel(validationBhv);
    %fprintf('Testing accuracy was: %.2f\n', validationAcc)
    
    % extract beta weights 
    L = size(SVMmodel.CodingMatrix,2); % Number of SVMs
    betaWx = zeros(length(SVMmodel.PredictorNames),L); 
    for j = 1:L
        SVM = SVMmodel.BinaryLearners{j};
        betaWx(:,j) = SVM.Beta;
    end
    
    % keep track per fold 
    validationScores(fold) = validationAcc;
    confusionMatrix(:,:,fold) = confusionmat(validationBhv, predictedBhv);
    betaWxMatrix(:,:,fold) = betaWx;

    clear svm; clear SVM; 
    
    clear trainingNeural; clear trainingBhv; clear betaWx;
    clear validationNeural; clear validationBhv; clear predictedBhv;
    clear SVMmodel;

end


meanValidationScores = mean(validationScores);
SEMValidationScores = std(validationScores)/sqrt(length(validationScores));
meanBetaWeights = mean(betaWxMatrix,3);

fprintf('Mean validation accuracy: %2f\n', mean(validationScores))

confusionPrct = zeros(3,3,nFolds);
for i = 1:nFolds
    confusionPrct(:,:,i) = diag(1./sum(confusionMatrix(:,:,i),2))*confusionMatrix(:,:,i);
end
confusionPrct = mean(confusionPrct,3);


function binnedSpikes = binSpikes(evtime, sptimes, win, tmin, tmax)
%%%%%%%%%
trainAx = tmin:win:tmax;
binnedSpikes = zeros(length(evtime),length(trainAx)-1); %initialize

for i=1:length(evtime)
    curT = evtime(i);
    spt = sptimes(find(sptimes >curT+tmin & sptimes < curT + tmax))-curT;%find spikes times after tmin & before tmax
                                                                         % subtracting curT converts spiketimes to -tmin:tmax scale 
                                                                         % (ie, first possible spktime is -tmin)
                                                                          
    binnedSpikes(i,:) = histcounts(spt, trainAx); %returns integer of # of spikes occurring in that bin


end
end

end
