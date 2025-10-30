%%% writen GI 211115 
%%% Calculates gesture preference index 
%%%  GPI = (nExp - [(sum (Ri) / rmax)]) ./ nExp -1 
%%% 
%%% Ri = spike rate for expression, i 
%%% Rmax = max Ri of any expression
%%% nExp = number of expressions compared 
%%%
%%% 0 indicates equal for all expressions, 1 indicates discharge for 1 exp
%%%     only 

close all;
clear all;

date = '171128'; %
subject ='Thor'; %
workdir = (['/Users/geena/Dropbox/PhD/SUAinfo/' subject '_' date '/Data4Analysis']);
chls = 1:192;
subsessions2plot = [ 1 2 4 5];
bhvs2plot = [1 2 4];

minRestFlag = 1;
minRest = 0.5; % in sec; minimal rest prior to move onset (trials to include)

colors = 'rbmgc';
saveFlag = 1;
threshlowFRFlag = 1;
threshlowFR = 0.1;

win = 0.001; % in sec
tmin = 0.5; 
tmax = 0.5; %0.5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate & then downsample gaussKern for convolution w/ spikes binned at "win" resolution
smParams.smoothSize = 50; %std in milisec
smParams.filterType = 1; % 1= gauss, 2= caus half-gauss, 3=box
smWin = genSmWin(smParams);
gaussKern = dwnsampleKernel(smWin,win);

extraBins = length(gaussKern); % extra bins to discard due to conv edges 
tmin = tmin +(extraBins*win); 
tmax = tmax + (extraBins*win);
edges = -tmin:win:tmax;  

% get rid of the edges where convolution results in artifacts & only use
% those bins later 
tmin2take = tmin - (extraBins*win); 
tmax2take = tmax - (extraBins*win);
taxis2take = -tmin2take:win:tmax2take;
bins2take = extraBins+1:length(taxis2take)+extraBins; 

% get relevant files 
flist = natsortfiles(dir([workdir '/bhvspikes_sub*.mat']));
flist = fullfile({flist.folder}, {flist.name}); 

% determine which type of subsession to plot
if strcmpi(subsessions2plot,'face expression')
    subs2show = flist(contains(flist, 'face expression'));
elseif strcmpi(subsessions2plot,'visual expression')
    subs2show = flist(contains(flist, 'visual expression'));
elseif strcmpi(subsessions2plot,'chew expression')
    subs2show = flist(contains(flist, 'chew expression'));
elseif strcmpi(subsessions2plot,'chew')
    subs2show = flist(contains(flist, 'chew'));
elseif strcmpi(subsessions2plot,'rest')
    subs2show = flist(contains(flist, 'rest'));
elseif strcmpi(subsessions2plot, 'all')
    subs2show = flist;
else 
    for i = 1:length(subsessions2plot)
       subs2show(i) = flist(contains(flist, ['sub' num2str(subsessions2plot(i)) '_']));
    end
end


% concatenate subsession spike times 
datAllSubs = [];
    
for session = 1:length(subs2show)   
    
    dat = [];
    obhvIn = load(subs2show{session});
    obhvIn = obhvIn.obhv;
    
    obhv = takeSpacedBhvs(obhvIn,minRest,minRestFlag); % remove trials w/ <500 ms rest prior move onset
    clear obhvIn;
  
    counter = 1;        % output index 

    [evindx, ~] = find(cell2mat(obhv.evScore(:,5)) == bhvs2plot); % indices of desired behavior
    if ~isempty(evindx) 

            for ch = chls
                for unit = 1:length(obhv.spikes.el(ch).sp)
                    if ~isempty(obhv.spikes.el(ch).sp(unit).times)


                        onsets = cell2mat(obhv.evScore(evindx,1));
                        spTimes = obhv.spikes.el(ch).sp(unit).times;
                        if strcmpi(subject,'Thor')
                            dat.binnedSpikes(:,:, counter) =   binSpikes( obhv.FT2dat.time(onsets), spTimes, win, tmin, tmax); %output is trials x bins
                        elseif strcmpi(subject,'Barney')
                            dat.binnedSpikes(:,:, counter) =   binSpikes( obhv.FT1dat.time(onsets), spTimes, win, tmin, tmax); %output is trials x bins
                        end
                        dat.name(counter,:) = string([num2str(ch) '-' num2str(unit)]);

                        counter = counter +1;

                    else % if no spiketimes for this unit, in this subsession, leave it empty

                        dat.binnedSpikes(counter,:,:) =  [];
                        dat.name(counter) = NaN;
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

    else %if there's no relevant behavior in subsession, leave it empty

        dat.frameON = [];
        dat.trialType = [];
        dat.condition = [];
        dat.comments = [];
    end

    datAllSubs = [dat; datAllSubs]; %concatenate across subsessions
    
end

% concatenate all trials across all subsessions 
% dimensions are trials x bins x cells
allSpikesAllSubs.binnedSpikes = (vertcat(datAllSubs(:).binnedSpikes));   
allSpikesAllSubs.name = datAllSubs.name;        

% behavioral / qualitative info (same for all cells) 
allTrials.frameON = vertcat(datAllSubs(:).frameON);
allTrials.trialType = vertcat(datAllSubs(:).trialType);
allTrials.condition = vertcat(datAllSubs(:).condition);
allTrials.comments = vertcat(datAllSubs(:).comments);

% convert spike-labels to strings for plotting later
 spikeLabels = (allSpikesAllSubs(:).name);
for i= 1:length(spikeLabels) 
    spikeLabels2plot(i,:) = str2double(strsplit(spikeLabels(i),'-'));
end 

% dimensions are ntrials x nbins x ncells
allCellsPSTH = zeros(size(allSpikesAllSubs.binnedSpikes));

for trial = 1:size(allSpikesAllSubs.binnedSpikes,1) %per trial 
    for cell = 1:size(allSpikesAllSubs.binnedSpikes,3) % per cell 
        allCellsPSTH(trial, :, cell) = conv(allSpikesAllSubs.binnedSpikes(trial, :, cell), gaussKern, 'same') ./win; %divide by win for FR
    end
end

%  average FRs; this is the trial-averaged response PER cell 
%  allCellsTrialAvgResponse.bhv is nTimepoints x nCells (averaged over
%       nReps)
for i = 1:length(bhvs2plot)
    allCellsTrialAvgResponse(i).bhv = squeeze(mean(allCellsPSTH((allTrials.trialType == bhvs2plot(i)),bins2take,:)));
    allCellsTrialAvgResponse(i).nReps = sum(allTrials.trialType == bhvs2plot(i));
end

%mean FRs per cell, per bhv 
for bhv = 1:length(bhvs2plot)
    meanFRs(:,bhv) = mean(allCellsTrialAvgResponse(bhv).bhv,1);
end

%remove low FR neurons 
if threshlowFRFlag
    thrCells2use = ~sum(meanFRs >= threshlowFR,2) == 0;
    spikeLabels = spikeLabels(thrCells2use,:);
    spikeLabels2plot = spikeLabels2plot(thrCells2use,:);
    
    for i = 1:length(bhvs2plot)
        allCellsTrialAvgResponse(i).bhv = allCellsTrialAvgResponse(i).bhv(:,thrCells2use);
    end
    
end

% assign cortical label 
[regions] = getChannel2CorticalRegionMapping(subject, 1);

for unit = 1:length(spikeLabels2plot)
    el = spikeLabels2plot(unit,1);

    if ismember(el, regions{1,1}.channels)
        cortRegion{unit} = 'S1'; % S1
    elseif ismember(el, regions{1,2}.channels)
        cortRegion{unit} = 'M1'; % M1
    elseif ismember(el, regions{1,3}.channels)
        cortRegion{unit} = 'PMv'; % PMv
    elseif ismember(el, regions{1,4}.channels)
        cortRegion{unit} = 'M3'; % M3
    else
        disp('warning!');
    end

end

%allSpikesAllSubs.binnedSpikes =  allSpikesAllSubs.binnedSpikes(:,bins2take,:);
%save([workdir '/binnedSpikesForMI.mat'],'allTrials', 'allSpikesAllSubs', 'spikeLabels', 'spikeLabels2plot','cortRegion','taxis2take','subject','date');

% finally do the preference index calculation
Ris = meanFRs(thrCells2use,:);
[Rmaxs, RmaxIndx] = max(Ris,[],2);
[Rmins, RminIndx] = min(Ris,[],2);
nExp = length(bhvs2plot);
faceExpPref = nan(1,length(Ris));

% index = (nExp - [(sum (Ri) / rmax)]) ./ nExp -1 
for unit = 1:length(Ris)
    faceExpPref(unit) = (nExp - (sum(Ris(unit,:)) ./ Rmaxs(unit))) ./ (nExp -1);
end

% calculate depth of preference index 
for unit = 1:length(Ris)
    depthPref(unit) = (Rmaxs(unit) - Rmins(unit) ) ./ (Rmaxs(unit) + Rmins(unit));
end

if saveFlag
   save([workdir '/' date '_' subject '_faceExpPref.mat'], 'Ris','cortRegion','depthPref','faceExpPref','spikeLabels', ...
       'spikeLabels2plot')
end

% plot faceExpPref
cmap = colorBeeSwarm(spikeLabels2plot);
figure; subplot 211; 
beeswarm(cmap', faceExpPref');
title([subject '-' date ' : faceExp Preference Index per cell']); %' subs = ' num2str(subsessions2plot)]);
subtitle(['bhvs= ' num2str(bhvs2plot)]);
xticks = 1:6;
xticklabels({'S1','M1m','M1lat','PMv','M3p','M3a'})

subplot 212; 
beeswarm(cmap', depthPref');
title([subject '-' date ' : depth Pref per cell']); %' subs = ' num2str(subsessions2plot)]);
subtitle(['bhvs= ' num2str(bhvs2plot)]);
xticks = 1:6;
xticklabels({'S1','M1m','M1lat','PMv','M3p','M3a'})



function cmap = colorBeeSwarm(spikeLabels2plot)

for k = 1:size(spikeLabels2plot,1)
    if spikeLabels2plot(k,1) <= 32
        cmap(k) = 1;
    elseif spikeLabels2plot(k,1) > 32 & spikeLabels2plot(k,1) <= 64
        cmap(k) = 2 ;
    elseif spikeLabels2plot(k,1) > 64 & spikeLabels2plot(k,1) <= 96
        cmap(k) = 3;
    elseif spikeLabels2plot(k,1) > 96 & spikeLabels2plot(k,1) <= 128
        cmap(k) = 4; 
    elseif spikeLabels2plot(k,1) > 128 & spikeLabels2plot(k,1) <= 160 
        cmap(k) = 5; %5
    elseif spikeLabels2plot(k,1) > 160 & spikeLabels2plot(k,1) <= 192 
        cmap(k) = 6; %6
    elseif spikeLabels2plot(k,1) > 192 & spikeLabels2plot(k,1) <= 225 
        cmap(k) = 7; %7
    elseif spikeLabels2plot(k,1) > 225 & spikeLabels2plot(k,1) <= 240 
        cmap(k) = 8; %8
    end
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
