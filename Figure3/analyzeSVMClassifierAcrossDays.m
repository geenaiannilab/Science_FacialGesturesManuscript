
%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%
clear all;

%% intial set-up 
subject2analyze = {'combined'};
denoiseFlag = true; % PCA denoising 
nFolds = 8; 
popType = 'balanced';
nUnitsPerRegion = '75';

workdir = ['/Users/geena/Dropbox/PhD/SUAinfo/PseudoPopulations/' popType '/N_' num2str(nUnitsPerRegion) 'cells/'];
outdir = ['/Users/geena/Dropbox/PhD/SUAinfo/PseudoPopulations/' popType '/N_' num2str(nUnitsPerRegion) 'cells/SVM_' subject2analyze{:} '/'];
saveResults = true;

classifyFlag = 1;
permuteFlag = 1;
nPerms = 100;
iterations = 50;
threshlowFRFlag = false; 
threshlowFR = 0.1;

% Timebase parameters for SVM 
desiredStart = -1.0;
desiredEnd = 1.0;
windowLength = 0.4; % 0.1, 0.2, 0.3, 0.4, 0.5
windowShiftSize = 0.05; % 0.05, 0.1

% repeat for each pseudopopulation
for iter = 1:iterations
 
    % load the pseudopopulation 
    pseudoPop = load([workdir 'twoSubjPseudopopulation_N' num2str(iter) '.mat']);
    
    % define some variables  
    arrayList = pseudoPop.regionList;
    bhvs2classify = unique(pseudoPop.pseudoPopulation.(subject2analyze{:}).bhvTrials);  
    originalTaxis = pseudoPop.tmin:pseudoPop.win:pseudoPop.tmax;
    
    % create time base
    firstWindowCenter = desiredStart + windowLength/2;
    lastWindowCenter = desiredEnd - windowLength/2;
    windowsCenter = firstWindowCenter:windowShiftSize:lastWindowCenter;
    
    
    %% Run the classifier
    % Loop over time windows
    if classifyFlag
        for ss = 1:length(subject2analyze)
    
            thisSubj = subject2analyze{ss}; 
            trialLabels = pseudoPop.pseudoPopulation.(thisSubj).bhvTrials;
    
            for aa = 1:length(arrayList)
    
                for tt = 1:length(windowsCenter)
                    tmin = windowsCenter(tt) - windowLength/2;
                    tmax = windowsCenter(tt) + windowLength/2;

                    minTidx = find(originalTaxis >= tmin,1,'first');
                    maxTidx = find(originalTaxis >= tmax,1,'first');
                    
                    dataIn = pseudoPop.pseudoPopulation.(thisSubj).(arrayList{aa})(:,minTidx:maxTidx,:);
                    dataIn = squeeze(sum(dataIn,2));
    
                        [ validationScoreStruct.(arrayList{aa})(tt),SEM.(arrayList{aa})(tt), betaScoreStruct.(arrayList{aa})(tt,:,:), codingMatrix] = SVMClassifierAcrossDays('dataIn', dataIn,...
                        'subject',thisSubj,'trialLabels', trialLabels, 'bhvs2classify',bhvs2classify, 'tmin', tmin, 'tmax', tmax', 'win', windowLength, 'array', arrayList{aa}, ...
                        'threshlowFRFlag', threshlowFRFlag,'threshlowFR', threshlowFR, 'denoiseFlag', denoiseFlag, 'nFolds',nFolds,'plotFlag', false, 'saveFig', false);
                            
                end
            end
        end
    
        if saveResults
            save([outdir '/' subject2analyze{:} '_N' num2str(iter) '_SVMresults_' num2str(windowLength) '_' num2str(windowShiftSize) '.mat'], ...
                'validationScoreStruct','SEM','betaScoreStruct', 'codingMatrix', 'desiredEnd','desiredStart','windowLength','windowShiftSize',...
                'denoiseFlag','nFolds','subject2analyze','popType');
        end
        close all; 
    end
    
    
    %% now run the permutations 
    if permuteFlag
        for ss = 1:length(subject2analyze)
    
            thisSubj = subject2analyze{ss}; 
            trialLabels = pseudoPop.pseudoPopulation.(thisSubj).bhvTrials;
    
            for aa = 1:length(arrayList)
    
                for tt = 1:length(windowsCenter)
                    tmin = windowsCenter(tt) - windowLength/2;
                    tmax = windowsCenter(tt) + windowLength/2;
                    
                    maxTidx = find(originalTaxis >= tmax,1,'first');
                    minTidx = find(originalTaxis >= tmin,1,'first');
                    
                    dataIn = pseudoPop.pseudoPopulation.(thisSubj).(arrayList{aa})(:,minTidx:maxTidx,:);
                    dataIn = squeeze(sum(dataIn,2));
    
                    [ permValidationScoreStruct.(arrayList{aa})(tt,:),permSEM.(arrayList{aa})(tt,:),~] = SVMClassifierPermutationTestingAcrossDays('dataIn', dataIn,...
                        'subject',thisSubj,'trialLabels', trialLabels, 'bhvs2classify',bhvs2classify, 'tmin', tmin, 'tmax', tmax', 'win', windowLength, 'array', arrayList{aa}, ...
                        'threshlowFRFlag', threshlowFRFlag,'threshlowFR', threshlowFR, 'denoiseFlag', denoiseFlag, 'nFolds',nFolds,'nPerms', nPerms);
                end
            end   
        end
        close all;
        
        if saveResults
            save([outdir '/' subject2analyze{:} '_N' num2str(iter) '_permSVMresults_' num2str(windowLength) '_' num2str(windowShiftSize) '.mat'], ...
                'permValidationScoreStruct','permSEM','nPerms','desiredEnd','desiredStart','windowLength','windowShiftSize',...
                'denoiseFlag','nFolds','subject2analyze','popType');
        end
    
    end
end


%% calc 95th%tile of permutations as significance cutoff
prcThres = 95;
P = zeros(length(arrayList),length(windowsCenter)); 

for tt = 1:length(windowsCenter)
    
    for aa = 1:length(arrayList)
        
        P(aa,tt) = prctile(permValidationScoreStruct.(arrayList{aa})(tt,:),prcThres);

    end
    
end

%% Plot results
plotFig = figure; hold on;
for aa = 1:length(arrayList)
    
    errorbar(windowsCenter, validationScoreStruct.(arrayList{aa})*100, SEM.(arrayList{aa})*100,'LineWidth',6); 
end
xlabel('Time (s)')
ylabel('Mean Decoding Accuracy (%)')
title(['Decoding Accuracy per Region -- ' num2str(windowLength) ' ms bins w/ ' num2str(windowShiftSize) ' ms overlap'])
legend(arrayList)

%%
figure;

for aa = 1:length(arrayList)
    subplot(2,2,aa)
    errorbar(windowsCenter, validationScoreStruct.(arrayList{aa})*100, SEM.(arrayList{aa})*100,'LineWidth',6); hold on
    plot(windowsCenter,P(aa,:)*100,'LineWidth',4); hold off
    ylim([40 105])
    title(['Decoding Accuracy -- ' arrayList{aa}])

end

% %% plot feature weights 
% for aa = 1:length(arrayList) 
%     figure; 
% 
%     n = size(pseudoPop.pseudoPopulation.(subject2analyze{:}).(arrayList{aa}),3);
%     for i = 1:3
%         subplot(1,3,i); ax(i) = imagesc(windowsCenter, 1:n,betaScoreStruct.(arrayList{aa})(:,:,i)');colormap(jet)
%         %compare = (((codingMatrix(:,i))) .* bhvs2classify')';
%         %title(['Hyperplane = ' num2str(compare(compare ~=0))])
%     end 
%     ax(1).Parent.YTick = 1:n;
%     %ax(1).Parent.YTickLabel = spNames.(arrayList{aa}){1,1};
%     sgtitle([arrayList{aa} ' SVMDecoder Feature Weights per Timebin'])
% 
% end
% 

