function [dat] = binSpikesTrialsPSID(subject, obhvIn,chls, bhvs, timeBase,varargin)

minRestFlag = 1;
minRest = 0.5;

% time parameters 
win = mean(diff(timeBase)); % win size for binning is same as bhvSR of camera (or less, if you downsampled)
if length(varargin) == 0    
    tmin = 0.5;
    tmax = 0.5;
else
    tmin = varargin{1};
    tmax = varargin{2};
end

% count n total possible cells 
n = 0;
for ch = chls
    for unit = 1:length(obhvIn.spikes.el(ch).sp)
        n = n+1;
    end
end

% take only trials w/ specified amount of rest prior 
obhv = takeSpacedBhvs(obhvIn,minRest,minRestFlag); 
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