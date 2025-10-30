%%% written GI 211215ish
%%%
%%% Implements SVM classifier (fitceoc) to classify individual trials of facial
%%% expression behavior, on basis of neural data.
%%% Crossvalidated by k=20 folds of 80/20 splits. 
%%% Saves validation scores & confusion matrices  


function [ meanValidationScores,SEMValidationScores,meanBetaWeights ] = SVMClassifierPermutationTestingAcrossDays(varargin)

%% Yes, you're getting an input parser
p = inputParser;

% Optional params
p.addParameter('dataIn','dataIn');
p.addParameter('subject','Barney',@ischar);
p.addParameter('array','M1',@ischar);
p.addParameter('nFolds',8,@isscalar);
p.addParameter('trialLabels','trialLabels');

p.addParameter('nPerms',10, @isscalar);

p.addParameter('denoiseFlag',false,@islogical); % PCA-based denoising 
p.addParameter('zScoreFlag', 1,@isscalar); 
p.addParameter('threshlowFRFlag', 1,@isscalar); % 0/1; whether or not to remove low FR neurons 
p.addParameter('threshlowFR', 0.01,@isscalar); %  sp/s threshold for removing low FR neurons 

p.addParameter('bhvs2classify', [1 2 4]);
p.addParameter('socialVNotDecoder',false)
p.addParameter('tmin', 0,@isscalar);
p.addParameter('tmax',1,@isscalar); % 1.5 to -0.5; 1 to 0, 0.5 to 0.5; 0 to 1 
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
dataIn = p.Results.dataIn;
nPerms = p.Results.nPerms;
denoiseFlag = p.Results.denoiseFlag;
bhvs2classify = p.Results.bhvs2classify;
socialVNotDecoder = p.Results.socialVNotDecoder;
trialLabels = p.Results.trialLabels; 

tmin = p.Results.tmin; 
tmax = p.Results.tmax; % 1.5 to -0.5; 1 to 0, 0.5 to 0.5; 0 to 1 
win = p.Results.win; % in sec
ax = -tmin:win:tmax;

plotFlag = p.Results.plotFlag; 
nFolds = p.Results.nFolds;
zScoreFlag = p.Results.zScoreFlag;  % zscoring matters ALOT here bc FRs are so different 
threshlowFRFlag = p.Results.threshlowFRFlag; % 0/1; whether or not to remove low FR neurons 
threshlowFR = p.Results.threshlowFR;

verboseFlag = p.Results.verboseFlag;

%% PREPROCESSING 
% sum spike counts over bins (now ntrials x ncells) 
meanFRs = mean(dataIn,1) / win ; % now in sp/s

% remove low FR neurons
if threshlowFRFlag
    dataIn = dataIn(:,(meanFRs >= threshlowFR));
end

% take top N PCs accounting for >98% variance, & recompose data 
if denoiseFlag 
     mu = mean(dataIn);
    [eigenvectors, scores,latent] = pca(dataIn);
    varEx = cumsum(latent/sum(latent));
    nComp = find(varEx >= 0.98,1, 'first'); % take top N PCs >= 95% variance
    Xhat = scores(:,1:nComp) * eigenvectors(:,1:nComp)'; 
    dataIn = bsxfun(@plus, Xhat, mu);
end

% zScore 
if zScoreFlag
    neuralInput = zscore(dataIn);
elseif ~zScoreFlag
    neuralInput = (dataIn);
end


%% permutate the outcome labels
shuffledTrialExpressions = nan(length(trialLabels),nPerms);
for p = 1:nPerms
    shuffledTrialExpressions(:,p) = trialLabels(randperm(length(trialLabels)));
end

%% fit & predict the permuted data 
for p = 1:nPerms

    permTrialData = shuffledTrialExpressions(:,p);

    %% split the data (data entry is nTrials x nCells, where each entry is # of spikes within a trial

    cv = cvpartition(permTrialData,'KFold',nFolds,'Stratify',true);
    for fold = 1:nFolds

        % define training & validation indices based on cvpartitions 
        trainingNeural = neuralInput(training(cv,fold),:);
        trainingBhv = permTrialData(training(cv,fold));

        validationNeural = neuralInput(test(cv,fold),:);
        validationBhv = permTrialData(test(cv,fold));
        
        % fit SVM on the training data
        t = templateSVM('Standardize',true);
        SVMmodel = fitcecoc(trainingNeural,trainingBhv,'Learners',t);
        codingMatrix = SVMmodel.CodingMatrix; 

        % predict for the validation data 
        [predictedBhv,~] = predict(SVMmodel,validationNeural);

         %compute accuracy 
        validationAcc = sum(validationBhv == predictedBhv) / numel(validationBhv);

         % extract beta weights 
        L = size(SVMmodel.CodingMatrix,2); % Number of SVMs
        betaWx = zeros(length(SVMmodel.PredictorNames),L); 
        for j = 1:L
            SVM = SVMmodel.BinaryLearners{j};
            betaWx(:,j) = SVM.Beta;
        end
   
        % keep track per fold 
        validationScores(p,fold) = validationAcc; 
        betaWxMatrix(:,:,fold) = betaWx;

    end
    clear trainingNeural; clear trainingBhv;
    clear validationNeural; clear validationBhv; clear predictedBhv;
    clear SVM; clear betaWx; clear SVMmodel;


    if verboseFlag 
        fprintf('Permuted mean validation accuracy: %2f\n', mean(validationScores,2))
    end

end

meanValidationScores = mean(validationScores,2);
SEMValidationScores = std(validationScores,[],2)/sqrt(length(validationScores));
meanBetaWeights = mean(betaWxMatrix,3);

function binnedSpikes = binSpikes(evtime, sptimes, win, tmin, tmax)
%%%%%%%%%
ax = tmin:win:tmax;
binnedSpikes = zeros(length(evtime),length(ax)-1); %initialize

for i=1:length(evtime)
    curT = evtime(i);
    spt = sptimes(find(sptimes >curT+tmin & sptimes < curT + tmax))-curT;%find spikes times after tmin & before tmax
                                                                         % subtracting curT converts spiketimes to -tmin:tmax scale 
                                                                         % (ie, first possible spktime is -tmin)
                                                                          
    binnedSpikes(i,:) = histcounts(spt, ax); %returns integer of # of spikes occurring in that bin


end
end

end
