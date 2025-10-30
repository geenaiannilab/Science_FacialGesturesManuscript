%%% written GI 211215ish
%%%
%%% Implements SVM classifier (fitceoc) to classify individual trials of facial
%%% expression behavior, on basis of single-trial PCA neural data.
%%% Current implementation uses 
%%% Crossvalidated by k=20 folds of 80/20 splits. 
%%% Saves validation scores & confusion matrices  


function [ meanValidationScores,SEMValidationScores, meanBetaWeights, codingMatrix ] = SVMClassifierAcrossDays(varargin)

%% Yes, you're getting an input parser
p = inputParser;

% Optional params
p.addParameter('dataIn','dataIn');
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
p.addParameter('tmin', 0,@isscalar);
p.addParameter('tmax',1,@isscalar); % 1.5 to -0.5; 1 to 0, 0.5 to 0.5; 0 to 1 
p.addParameter('win', 0.02,@isscalar); % in sec

p.addParameter('plotFlag', false); 
p.addParameter('saveFig', false);
p.addParameter('nClusters', 3,@isscalar);

% parse
p.parse(varargin{:})

close all; 

subject = p.Results.subject;
denoiseFlag = p.Results.denoiseFlag;
dataIn = p.Results.dataIn;
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


%% PREPROCESSING 
% sum spike counts over bins (now ntrials x ncells) 
meanFRs = mean(dataIn,1) / win ; % now in sp/s

% remove low FR neurons
if threshlowFRFlag
    dataIn = dataIn(:,(meanFRs >= threshlowFR));
end

% take top N PCs accounting for >95% variance, & recompose data 
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

% data entry is nTrials x nCells , 
%    where each entry is # of spikes within a trial

%% split the data 
cv = cvpartition(trialLabels,'KFold',nFolds,'Stratify',true);

for fold = 1:nFolds

    % define training & validation indices based on cvpartitions 
    trainingNeural = neuralInput(training(cv,fold),:);
    trainingBhv = trialLabels(training(cv,fold));

    validationNeural = neuralInput(test(cv,fold),:);
    validationBhv = trialLabels(test(cv,fold));

    % fit SVM on the training data
    t = templateSVM('Standardize',true,'SaveSupportVectors','on');
    SVMmodel = fitcecoc(trainingNeural,trainingBhv,'Learners',t);
    codingMatrix = SVMmodel.CodingMatrix; 

    % predict for the validation data 
    % predBhvScore is an testTrials-by-nClasses matrix of soft scores; 
     [predictedBhv,~] = predict(SVMmodel,validationNeural);

    % compute accuracy 
    validationAcc = sum(validationBhv == predictedBhv) / numel(validationBhv);
    fprintf('Testing accuracy was: %.2f\n', validationAcc)
   
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

    % plot the training data, clusters, & predicted/validated data 
    if plotFlag 
    
        % extract each learner 
        L = size(SVMmodel.CodingMatrix,2); % Number of SVMs
        sv = cell(L,1); % Preallocate for support vector indices
        for j = 1:L
            SVM = SVMmodel.BinaryLearners{j};
            sv{j} = SVM.SupportVectors;
            sv{j} = sv{j}.*SVM.Sigma + SVM.Mu;
        end
      
        figure; 
        colors = 'rbgk';
        grid on; hold on; 
    
        %
        bhvs2classify = unique(trainingBhv)'; % training data 
        for i = 1:length(bhvs2classify)
            thisBhv = bhvs2classify(i);
            idx = thisBhv == trainingBhv;
            scatter3(trainingNeural(idx,1),trainingNeural(idx,2),trainingNeural(idx,3), 30,colors(i),'filled')
        end
        hold on
    
        markers = {'ko','mo','yo'}; % support vectors 
        for j = 1:L
            svs = sv{j};
            plot3(svs(:,1),svs(:,2),svs(:,3), markers{j},'LineWidth',3,...
                'MarkerSize',10 + (j - 1)*3);
        end
        
        colors = 'kmy';
        for i = 1:length(SVMmodel.BinaryLearners)
            svm = SVMmodel.BinaryLearners{i};
            
            % Step 1: Choose hyperplane's x1 values
            x1vals = linspace(min(trainingNeural(:,1)),max(trainingNeural(:,1)));
             
            % Step 2: Solve for ((x2-svm.Mu(2))/svm.Sigma(2)) in terms of x1.
            x2fun =@(x1) -svm.Beta(1)/svm.Beta(2) * ((x1-svm.Mu(1))./svm.Sigma(1)) - svm.Bias/svm.Beta(2);
            
            % Step 3: Calculate corresponding x2 values
            x2vals = x2fun(x1vals) * svm.Sigma(2) + svm.Mu(2);
           
            l{i} = mat2str(SVMmodel.ClassNames(SVMmodel.CodingMatrix(:,i)~=0)');
    
            % plot the hyperplane
            plot(x1vals,x2vals,colors(i),'LineWidth',6); hold on; 
        end
    
        hold off 
    
        view(3)
        title('Face Exp -- MultiSVM Hyperplanes &  Support Vectors')
        legend; %legend(['bhv1','bhv2','bhv4' l l])
        pause; close all ; 
        
        clear svm; clear SVM; 
    end

    clear trainingNeural; clear trainingBhv; clear betaWx;
    clear validationNeural; clear validationBhv; clear predictedBhv;
    clear SVMmodel;

end


meanValidationScores = mean(validationScores);
SEMValidationScores = std(validationScores)/sqrt(length(validationScores));
meanBetaWeights = mean(betaWxMatrix,3);

fprintf('Mean validation accuracy: %2f\n', mean(validationScores))

for i = 1:nFolds
    confusionPrct(:,:,i) = diag(1./sum(confusionMatrix(:,:,i),2))*confusionMatrix(:,:,i);
end

confusionPrct = mean(confusionPrct,3);


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
