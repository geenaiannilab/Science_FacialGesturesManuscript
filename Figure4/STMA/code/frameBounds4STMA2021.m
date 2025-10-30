function [spikeTimesAdj, framesLowerBound, framesUpperBound] = frameBounds4STMA2021(spikeTimes, cameraSync, bhvSR, directionality, windowSize)

%PURPOSE: generate 2 vectors, of frame #s, that correspond to the 1 second
            % periods prior to each spike time, of a single neuron; these
            % frame "windows" will be used to calculate the average
            % stimulus within those windows 
%INPUT: spikeTimes = for a single neuron 
%       cameraSync = FT1 or FT2 sync pulse, depending on subj/face camera
%       bhvSR = sampling rate of bhv; ie., camera frame rate
%       windowSize = period to examine prior to a spike (in this version,
%           it is indeed PRIOR; will need to adjust code for AFTER) in SECONDS 
%       directionality = 'pre' or 'post'; string that determines if you're
%                   looking at 1 sec windows prior to spikes or after spikes 
%       

%OUTPUT: framesLowerBound = 1 x nSpikes vector of frame#s that correspond to
%                           specified window prior to each spike time 
%        framesUpperBound = 1 x nSpikes vector of frame#s that correspond
%                           to specified window prior to each spike time

if strcmpi(directionality, 'pre')
    pre = 1;
    post = 0;
elseif strcmpi(directionality, 'post')
    pre = 0;
    post = 1;
else
    disp('check directionality input; may be either PRE or POST');
end

if pre % examine stimulus window BEFORE spike times
    spikeTimes((spikeTimes <= cameraSync(1) + windowSize)) = 0; %delete spike events that occurred prior to camera pulse & within first window
    spikeTimes((spikeTimes >= cameraSync(end) - (1/bhvSR))) = 0; %delete those after camera pulse turned off
end

if post % examine stimulus winodw AFTER spike times
    spikeTimes((spikeTimes <= cameraSync(1) + (1/bhvSR))) = 0; %delete spike events that occurred prior to camera pulse & within first window
    spikeTimes((spikeTimes >= cameraSync(end) - (windowSize+(1/bhvSR)))) = 0; %delete those after camera pulse turned off
end

spikeIndices = find(spikeTimes)'; % indices of remaining spikes 
nSpikes = numel(spikeIndices); % number of spikes occuring within camera pulse 

% preallocate
wind_lowerBound = zeros(1,nSpikes);
wind_upperBound = zeros(1,nSpikes);
framesLowerBound = zeros(1,nSpikes);
framesUpperBound = zeros(1,nSpikes);

counter = 1;

%for each spike that occurred within the camera pulse, calculate a 1 second
        %window PRIOR OR POST to each spike time, and then find the frame #s that
        %correspond to these one second windows
        
% you only need to do this ONCE per neuron, (since regardless of the stimulus/pixel intensity time series you're looking at, 
%    the neuron's spike times & thus calculated frames to examine, are constant
if pre
    for i = spikeIndices 

        wind_lowerBound(counter) = spikeTimes(i) - windowSize; %in seconds 
        framesLowerBound(counter) = find( cameraSync >= wind_lowerBound(counter), 1,'first'); 
        wind_upperBound(counter) = spikeTimes(i); %maybe this should be moved up closer to spike_time(i)?
        framesUpperBound(counter) = find( cameraSync >= wind_upperBound(counter), 1,'first'); %in seconds
        counter = counter + 1;
    end 
end

if post 
    for i = spikeIndices 

        wind_lowerBound(counter) = spikeTimes(i); %in seconds 
        framesLowerBound(counter) = find( cameraSync >= wind_lowerBound(counter), 1,'first'); 
        wind_upperBound(counter) = spikeTimes(i) + windowSize; 
        framesUpperBound(counter) = find( cameraSync >= wind_upperBound(counter), 1,'first'); %in seconds
        counter = counter + 1;
    end 
end

spikeTimesAdj = spikeTimes; 
