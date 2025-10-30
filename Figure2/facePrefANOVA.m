%% Calculating for  Figure 2E,  
%% generates '(date_subj_)facePrefAnova,.mat' (needs to be run over all days) 
% this is a two-way ANOVA to look for effect of gesture, time on individual
% cell's FR 
% this data then get loaded by chiSquaredGOF.m to be thresholded, and
% indiviudal Chi-squared tests applied 
% 

%% these data are then plotted by plotBhvTimeANOVAbyRegion.m 
% to generate Figure 2E 

close all;
clear all;

date = '171128';
subject ='Thor'; 
workdir = (['/Users/geena/Dropbox/PhD/SUAinfo/' subject '_' date '/Data4Analysis']);
chls = 1:192;
subsessions2plot = [ 1 2 4 5];
bhvs2plot = [1 2 4];

minRestFlag = 1;
minRest = 0.5; % in sec; minimal rest prior to move onset (trials to include)

colors = 'rbmgc';
threshlowFRFlag = 1;
threshlowFR = 0.1;
saveFlag = 1;

win = 0.5; % in sec
tmin = 0.5;
tmax = 0.5;
edges = -tmin:win:tmax;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

if strcmpi(subject, 'Thor')
    [allSpikesAllSubs, allTrials] = getRidOfChewTrials(allSpikesAllSubs, allTrials);
end

% convert spike-labels to strings for plotting later
 spikeLabels = (allSpikesAllSubs(:).name);
for i= 1:length(spikeLabels) 
    spikeLabels2plot(i,:) = str2double(strsplit(spikeLabels(i),'-'));
end 

nTotalTrials = size(allSpikesAllSubs.binnedSpikes,1);
nTimeBins = size(allSpikesAllSubs.binnedSpikes,2);
totalDataPtsIn = nTotalTrials * nTimeBins;
timeFactors= 1:nTimeBins;
trialTypeIn = repmat(allTrials.trialType',1,nTimeBins);
timeBinIn = sort(repmat(timeFactors,1,nTotalTrials));
  
for unit = 1:size(allSpikesAllSubs.binnedSpikes,3)

    spikes = allSpikesAllSubs.binnedSpikes(:,:,unit);
    anovaIN = reshape(spikes,1, numel(spikes)); 
    [pvalOut,anovaTblOut, statsOut] = anovan(anovaIN, {trialTypeIn timeBinIn},'model',2,'varnames',{'bhv','time'},'display','off');
    allPvals{unit} = pvalOut;
    allResults{unit} = anovaTblOut;
    allStats{unit} = statsOut;
    clear pvalOut; clear anovaTblOut; clear statsOut; clear spikes; clear anovaIN;
   
end

for ii = 1:length(allResults)
    poolPvals(ii,:) = allPvals{1,ii};
end

spikeLabels2plot =spikeLabels2plot(~isnan(poolPvals(:,1)),:);

poolPvals = poolPvals(~isnan(poolPvals(:,1)),:);
%figure; subplot 311; histogram(poolPvals(:,1),histedges); subplot 312; histogram(poolPvals(:,2),histedges); subplot 313; histogram(poolPvals(:,3),histedges);

if saveFlag
    save([workdir '/' date '_' subject '_facePrefAnova.mat'], 'allPvals', 'allResults','allStats','spikeLabels2plot')
end



cmap = colorBeeSwarm(spikeLabels2plot);
figure; beeswarm(cmap', poolPvals'); %ylim([-0.1 0.5])
title([subject '-' date ' subs = ' num2str(subsessions2plot)]);
subtitle(['bhvs= ' num2str(bhvs2plot)]);


function cmap = colorBeeSwarm(spikeLabels2plot)

for k = 1:size(spikeLabels2plot,1)
    if spikeLabels2plot(k,1) <= 32
        cmap(k) = 1;
    elseif spikeLabels2plot(k,1) > 32 & spikeLabels2plot(k,1) <= 64
        cmap(k) = 2;
    elseif spikeLabels2plot(k,1) > 64 & spikeLabels2plot(k,1) <= 96
        cmap(k) = 3;
    elseif spikeLabels2plot(k,1) > 96 & spikeLabels2plot(k,1) <= 128
        cmap(k) = 4;
    elseif spikeLabels2plot(k,1) > 128 & spikeLabels2plot(k,1) <= 160 
        cmap(k) = 5;
    elseif spikeLabels2plot(k,1) > 160 & spikeLabels2plot(k,1) <= 192 
        cmap(k) = 6;
    elseif spikeLabels2plot(k,1) > 192 & spikeLabels2plot(k,1) <= 225 
        cmap(k) = 7;
    elseif spikeLabels2plot(k,1) > 225 & spikeLabels2plot(k,1) <= 240 
        cmap(k) = 8;
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
