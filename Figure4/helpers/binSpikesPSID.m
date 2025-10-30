function [dat] = binSpikesPSID(varargin)

p = inputParser;
p.addParameter('timeBase','timeBase')
p.addParameter('obhv', 'obhv');
p.addParameter('chls',1:240,@isnumeric);
p.addParameter('sqrtFlag',0,@isscalar);

% parse
p.parse(varargin{:})
close all; 

timeBase = p.Results.timeBase;
obhv = p.Results.obhv;
chls = p.Results.chls;
sqrtFlag = p.Results.sqrtFlag;

% count n total possible cells 
n = 0;
for ch = chls
    for unit = 1:length(obhv.spikes.el(ch).sp)
        n = n+1;
    end
end

% preallocate
counter = 1;
dat.binnedSpikes = zeros(length(timeBase)-1,n);
dat.name = cell(n,1);

% for each possible cell, bin spiketimes & extract channel/unit name
for ch = chls
    for unit = 1:length(obhv.spikes.el(ch).sp)
        if ~isnan(obhv.spikes.el(ch).sp(unit).times)

            spTimes = obhv.spikes.el(ch).sp(unit).times;
            
            dat.binnedSpikes(:,counter) = histcounts(spTimes, timeBase); %returns integer of # of spikes occurring in that bin
            dat.name{counter} = string([num2str(ch) '-' num2str(unit)]);

            counter = counter +1;

        else % if no spiketimes for this unit, in this subsession, leave it empty

            dat.binnedSpikes(:,counter) = NaN;
            dat.name{counter} = NaN;
            counter = counter + 1;

        end
    end
end

% remove columns without spikes
idx2keep = ~any(isnan(dat.binnedSpikes));
dat.binnedSpikes = dat.binnedSpikes(:,idx2keep);
dat.name = dat.name(idx2keep);

if sqrtFlag
    dat.binnedSpikes = sqrt(dat.binnedSpikes);
end

end