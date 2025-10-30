function [output] = preprocBhvTrialData(varargin)

p = inputParser;

% params
p.addParameter('bhvMeasure', 'bhvMeasure');
p.addParameter('markerPCA', false, @islogical);
p.addParameter('bhvPCs',1:12);
p.addParameter('timeBase','timeBase')
p.addParameter('obhvIn', 'obhvIn');
p.addParameter('bhvs',[1 2]);
p.addParameter('tmin',0.5, @isnumeric);
p.addParameter('tmax',0.5, @isnumeric);
p.addParameter('minRestFlag',true, @islogical);
p.addParameter('minRest',0.5, @isnumeric);
p.addParameter('smoothFlag',false, @islogical)
p.addParameter('smoothSize',20, @isscalar)
p.addParameter('filterType',1,@isscalar)

% parse
p.parse(varargin{:})

%subject = p.Results.subject;
bhvMeasure = p.Results.bhvMeasure;
markerPCA = p.Results.markerPCA;
bhvPCs = p.Results.bhvPCs;
timeBase = p.Results.timeBase;
obhvIn = p.Results.obhvIn;
bhvs =  p.Results.bhvs;
tmin = p.Results.tmin;
tmax = p.Results.tmax;
minRestFlag = p.Results.minRestFlag;
minRest = p.Results.minRest;
smParams.smoothSize = p.Results.smoothSize; % std of gaussian kernel in milisec
smParams.filterType = p.Results.filterType; % 1= gaussian, 2= causal, 3= boxcar
smoothFlag = p.Results.smoothFlag;

% time parameters 
win = mean(diff(timeBase)); % win size for binning is same as bhvSR of camera (or less, if you downsampled)

% indices (in frames) to take per trial 
frMin = floor (tmin / win);
frMax = floor (tmax / win); 
%stepSize = 1; % win represents 1 frame

% take only trials w/ specified amount of rest prior 
obhv = spacedBhvsPSID(obhvIn,minRest,minRestFlag); 
clear obhvIn;

% indices & frame onsets of desired bhv
[evindx, ~] = find(cell2mat(obhv.evScore(:,5)) == bhvs); 
onsets = cell2mat(obhv.evScore(evindx,1));
   
% preallocate
output = cell(length(onsets),1);

% reformat data for PSID 

if smoothFlag
    smWin = genSmWin(smParams);
    gaussKern = dwnsampleKernel(smWin,win);

    for trial = 1:length(onsets)
        thisOnset = (onsets(trial) - frMin); %
        thisOffset =  (onsets(trial) + frMax); %
        data = bhvMeasure(thisOnset:thisOffset,:);
        
        for i = 1:size(bhvMeasure,2)
            smdata(:,i) = conv(data(:,i), gaussKern,'same');
        end
        
        output{trial} = smdata;
    end

elseif markerPCA 
     allTrials = [];
     for trial = 1:length(onsets)
        thisOnset = (onsets(trial) - frMin); %subtract 1 since taking Diff in Frames
        thisOffset =  (onsets(trial) + frMax); %subtract 1 since taking Diff in Frames
        thisTrial = bhvMeasure(thisOnset-1:thisOffset-1,:);
        allTrials = cat(1,allTrials, thisTrial);
     end

    % PCA input is nTimepoints x nMarkers,
    [coeff,score,latent, ~] = pca(allTrials); % centering (subtract column/cell means) done automatically by pca()
    figure; plot(cumsum(latent) ./ sum(latent));
    figure; ax = imagesc(coeff(:,1:10)); colormap(bluewhitered);
    yticks =1:2:size(allTrials,2);
    
    % take top N PCs & reshape them by trial 
    score = score(:,bhvPCs);
    tmp = reshape(score, [], length(onsets),length(bhvPCs));
    for trial = 1:length(onsets)
        output{trial} = squeeze(tmp(:,trial,:));
    end

else

    for trial = 1:length(onsets)
        thisOnset = (onsets(trial) - frMin); %subtract 1 since taking Diff in Frames
        thisOffset =  (onsets(trial) + frMax); %subtract 1 since taking Diff in Frames
        output{trial} = bhvMeasure(thisOnset-1:thisOffset-1,:); % ??? cut off the extra to match neuraldata
     end
end

end