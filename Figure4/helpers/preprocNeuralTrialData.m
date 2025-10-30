function [output, dat] = preprocNeuralTrialData(varargin)

p = inputParser;

% Optional params
p.addParameter('subject', 'subject');
p.addParameter('timeBase','timeBase')
p.addParameter('obhvIn', 'obhvIn');
p.addParameter('chls',1:240,@isnumeric);
p.addParameter('bhvs',2);
p.addParameter('sqrtFlag', 0, @isscalar);
p.addParameter('minRestFlag',true, @islogical);
p.addParameter('minRest',0.5, @isnumeric);
p.addParameter('tmin',0.5, @isnumeric);
p.addParameter('tmax',0.5, @isnumeric);
p.addParameter('smoothFlag',false, @islogical)
p.addParameter('smoothSize',20, @isscalar)
p.addParameter('filterType',1,@isscalar)

% parse
p.parse(varargin{:})
%close all; 
%

subject = p.Results.subject;
obhvIn = p.Results.obhvIn;
chls =  p.Results.chls;
bhvs =  p.Results.bhvs;
timeBase = p.Results.timeBase;
tmin = p.Results.tmin;
tmax = p.Results.tmax;
minRestFlag = p.Results.minRestFlag;
minRest = p.Results.minRest;
sqrtFlag = p.Results.sqrtFlag;
smParams.smoothSize = p.Results.smoothSize; % std of gaussian kernel in milisec
smParams.filterType = p.Results.filterType; % 1= gaussian, 2= causal, 3= boxcar
smoothFlag = p.Results.smoothFlag;

% time parameters 
win = mean(diff(timeBase)); % win size for binning is same as bhvSR of camera (or less, if you downsampled)

% count n total possible cells 
n = 0;
for ch = chls
    for unit = 1:length(obhvIn.spikes.el(ch).sp)
        n = n+1;
    end
end

% take only trials w/ specified amount of rest prior 
%obhv = takeSpacedBhvs(obhvIn,minRest,minRestFlag); 
obhv = spacedBhvsPSID(obhvIn,minRest,minRestFlag); 
clear obhvIn;

% indices & frame onsets of desired bhv
[evindx, ~] = find(cell2mat(obhv.evScore(:,5)) == bhvs); 
onsets = cell2mat(obhv.evScore(evindx,1));

% preallocate
dat.binnedSpikes = zeros(length(onsets),length(-tmin:win:tmax)-1,n);
dat.name = cell(n,1);
counter = 1;        % output index 

for ch = chls
    for unit = 1:length(obhv.spikes.el(ch).sp)
        if ~isnan(obhv.spikes.el(ch).sp(unit).times)


            spTimes = obhv.spikes.el(ch).sp(unit).times;
            if strcmpi(subject,'Thor')
                dat.binnedSpikes(:,:, counter) =   binSpikes( obhv.FT2dat.time(onsets), spTimes, win, tmin, tmax); %output is trials x bins
            elseif strcmpi(subject,'Barney') || strcmpi(subject,'Dexter')
                dat.binnedSpikes(:,:, counter) =   binSpikes( obhv.FT1dat.time(onsets), spTimes, win, tmin, tmax); %output is trials x bins
            end
            dat.name{counter} = string([num2str(ch) '-' num2str(unit)]);

            counter = counter +1;

        else % if no spiketimes for this unit, in this subsession, leave it empty

            dat.binnedSpikes(:,:,counter) = NaN;
            dat.name{counter} = NaN;
            counter = counter + 1;

        end
    end
end

% after loop thru all cells, add the immutable behavioral
% data
dat.frameON = cell2mat(obhv.evScore(evindx,1));
dat.trialType = cell2mat(obhv.evScore(evindx,5));
dat.condition = (obhv.evScore(evindx,6));
dat.comments = obhv.evScore(evindx,8);


% remove columns without spikes
idx2keep = cellfun(@isstring,dat.name);
dat.binnedSpikes = dat.binnedSpikes(:,:,idx2keep);
dat.name = dat.name(idx2keep);

% smooth with a 50 ms gaussian kernel  
% (dimensions are ntrials x nbins x ncells) 
if smoothFlag
    smWin = genSmWin(smParams);
    gaussKern = dwnsampleKernel(smWin,win);
    allCellsPSTH = zeros(size(dat.binnedSpikes));

    for trial = 1:size(dat.binnedSpikes,1) %per trial 
        for unit = 1:size(dat.binnedSpikes,3) % per cell 
            allCellsPSTH(trial, :, unit) = conv(dat.binnedSpikes(trial, :, unit), gaussKern, 'same') ./win; %divide by win for FR
        end
    end
    dataOut = allCellsPSTH;

else
    dataOut =  dat.binnedSpikes;
end

% reformat data for PSID
output = cell(length(onsets),1);
for trial = 1:length(onsets)
    if sqrtFlag
        output{trial} = sqrt(squeeze(dataOut(trial,:,:)));
    else
        output{trial} = squeeze(dataOut(trial,:,:));
    end

end


function binnedSpikes = binSpikes(evtime, sptimes, win, tmin, tmax)
%%%%%%%%%
ax = -tmin:win:tmax;
binnedSpikes = zeros(length(evtime),length(ax)-1); %initialize

for i=1:length(evtime)
    curT = evtime(i);
    spt = sptimes(find(sptimes >curT-tmin & sptimes < curT + tmax))-curT;%find spikes times after tmin & before tmax
                                                                         % subtracting curT converts spiketimes to -tmin:tmax scale 
                                                                         % (ie, first possible spktime is -tmin)
                                                                          
    binnedSpikes(i,:) = histcounts(spt, ax); %returns integer of # of spikes occurring in that bin


end
end

end