function [output] = preprocBhvTrialDataDLC(varargin)

p = inputParser;

% params
p.addParameter('bhvMeasure', 'bhvMeasure');
p.addParameter('timeBase','timeBase')
p.addParameter('obhvIn', 'obhvIn');
p.addParameter('bhvs',[1 2]);
p.addParameter('tmin',0.5, @isnumeric);
p.addParameter('tmax',0.5, @isnumeric);
p.addParameter('minRestFlag',false, @islogical);
p.addParameter('minRest',0.5, @isnumeric);
p.addParameter('smoothFlag',false, @islogical)
p.addParameter('smoothSize',50, @isscalar)
p.addParameter('filterType',1,@isscalar)

% parse
p.parse(varargin{:})
close all; 

%subject = p.Results.subject;
bhvMeasure = p.Results.bhvMeasure;
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
smWin = genSmWin(smParams);
gaussKern = dwnsampleKernel(smWin,win);

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
    for trial = 1:length(onsets)
        thisOnset = (onsets(trial) - frMin); %subtract 1 since taking Diff in Frames
        thisOffset =  (onsets(trial) + frMax); %subtract 1 since taking Diff in Frames
        data = bhvMeasure(thisOnset:thisOffset,:);
        %data = bhvMeasure(thisOnset-1:thisOffset-1,:); % ??? cut off the extra to match neuraldata
        
        for i = 1:size(bhvMeasure,2)
            smdata(:,i) = conv(data(:,i), gaussKern,'same');
        end
        
        output{trial} = smdata;
    end

else 
     for trial = 1:length(onsets)
        thisOnset = (onsets(trial) - frMin); %subtract 1 since taking Diff in Frames
        thisOffset =  (onsets(trial) + frMax); %subtract 1 since taking Diff in Frames
        output{trial} = bhvMeasure(thisOnset-1:thisOffset-1,:); % ??? cut off the extra to match neuraldata
     end
end

end