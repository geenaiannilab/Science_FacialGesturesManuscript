
%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%
clear all;

%% intial set-up
subject2analyze = {'combined'};
denoiseFlag = true; % PCA denoising
nFolds = 8;
popType = 'balanced';
nUnitsPerRegion = 50;

workdir = ['/Users/geena/Dropbox/PhD/SUAinfo/PseudoPopulations/' popType '/N_' num2str(nUnitsPerRegion) 'cells/'];
outdir = ['/Users/geena/Dropbox/PhD/SUAinfo/PseudoPopulations/' popType '/N_' num2str(nUnitsPerRegion) 'cells/'];

saveResults = true; 
saveFigsFlag = true; 

iterations = 50;
threshlowFRFlag = true; 
threshlowFR = 0.1;
minRestFlag = 1;

% Timebase parameters for SVM 
desiredStart = -1.0;
desiredEnd = 1.0;
windowLength = 0.2; % 0.1, 0.2, 0.3, 0.4, 0.5
windowShiftSize = 0.01; % 0.05, 0.1
minRest = abs(desiredStart);

for ss = 1:length(subject2analyze)

    thisSubj = subject2analyze{ss};

    % repeat for each pseudopopulation
    for iter = 1:iterations
        
        % load the pseudopopulation
        pseudoPop = load([workdir 'twoSubjPseudopopulation_N' num2str(iter) '.mat']);

        % define some variables
        arrayList = pseudoPop.regionList;
        bhvs2classify = unique(pseudoPop.pseudoPopulation.(subject2analyze{:}).bhvTrials);
        originalTaxis = pseudoPop.tmin:pseudoPop.win:pseudoPop.tmax;
        trialLabels = pseudoPop.pseudoPopulation.(thisSubj).bhvTrials;

        % create timebase
        firstWindowCenter = desiredStart + windowLength/2;
        lastWindowCenter = desiredEnd - windowLength/2;
        windowsCenter = firstWindowCenter:windowShiftSize:lastWindowCenter;

        %% Run the classifier
        % Loop over time windows
        for aa = 1:length(arrayList)

            for tt = 1:length(windowsCenter)
                trainTmin = windowsCenter(tt) - windowLength/2;
                trainTmax = windowsCenter(tt) + windowLength/2;

                minTrainidx = find(originalTaxis >= trainTmin,1,'first');
                maxTrainidx = find(originalTaxis >= trainTmax,1,'first');

                trainNeural = pseudoPop.pseudoPopulation.(thisSubj).(arrayList{aa})(:,minTrainidx:maxTrainidx,:);
                trainNeural = squeeze(sum(trainNeural,2));

                for vv = 1:length(windowsCenter)
                    testTmin = windowsCenter(vv) - windowLength/2;
                    testTmax = windowsCenter(vv) + windowLength/2;

                    minTestidx = find(originalTaxis >= testTmin,1,'first');
                    maxTestidx = find(originalTaxis >= testTmax,1,'first');

                    testNeural = pseudoPop.pseudoPopulation.(thisSubj).(arrayList{aa})(:,minTestidx:maxTestidx,:);
                    testNeural = squeeze(sum(testNeural,2));

                    [ validationScoreStruct.(arrayList{aa})(tt,vv,iter),SEM.(arrayList{aa})(tt,vv,iter), ~] = crossTemporalSVMClassifierAcrossDays('trainNeural', trainNeural,'testNeural',testNeural,...
                        'subject', thisSubj, 'trialLabels', trialLabels, 'bhvs2classify', bhvs2classify,...
                        'denoiseFlag', denoiseFlag, 'minRestFlag', minRestFlag, 'minRest', minRest, 'array', arrayList{aa},...
                        'trainTmin', trainTmin, 'trainTmax', trainTmax, 'testTmin', testTmin, 'testTmax', testTmax,...
                        'plotFlag', false, 'nFolds',nFolds,...
                        'threshlowFRFlag', threshlowFRFlag,'threshlowFR', threshlowFR);

                    
                    disp(['Iteration --- ' num2str(iter)]);
        
                end % end test windows
            end % end train windows
        end % end array list
    end % end iterations
    
    save([outdir '/CTD_' thisSubj '/smallerOverlap/CTD_start' num2str(abs(desiredStart)) '_bin' num2str(windowLength) '.mat'],...
        'nFolds','denoiseFlag','thisSubj','threshlowFRFlag','threshlowFR','minRestFlag','minRest',...
        'desiredStart','desiredEnd','arrayList','bhvs2classify','originalTaxis',...
            'validationScoreStruct','SEM','windowLength','windowsCenter','windowShiftSize','-v7.3');

end % end subjects


