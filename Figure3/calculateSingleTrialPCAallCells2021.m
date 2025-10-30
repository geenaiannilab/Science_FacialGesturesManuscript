% written GI 
%
% PUPRPOSE: Pool spikes to calculate neural state space
%   where each point in reduced neural state-space represents a TRIAL 
%   of a particular BEHAVIOR TYPE
%   
% Generates input to PlotSingelTrialsNeuralStateSpace, which is allCellsSpikeCounts (nCells x nTrials) 
%    where each entry is # of spikes that occured on that trial
%

clear all;

date = '210704';
subject ='Barney';
workdir = (['/Users/geena/Dropbox/PhD/SUAinfo/' subject '_' date '/Data4Analysis']);
saveFlag = true; 

chls = [1:240];
subsessions2plot = [2:5 7:10];
bhvs2plot = [1 2 4];
win = 0.02; % in sec
tmin = 0.5; tmax = 1; % 1.5 to -0.5; 1 to 0, 0.5 to 0.5
ax = -tmin:win:tmax;
minRestFlag = 1;
minRest = abs(tmin); % in sec; minimal rest prior to move onset (trials to include)

threshlowFRFlag = 1; % 0/1;  remove low FR neurons 
threshlowFR = 0.01; %  sp/s threshold 

sqrtFlag = 1;
centerOnlyFlag = 0; % ONLY mean center input to PCA (dont divide by std)
centerNormalizeFlag = 1; % normalize cells/variables by their variance in addition to mean-centering (PCA on correlation matrix)    

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
                        elseif strcmpi(subject,'Barney') || strcmpi(subject,'Dexter')
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

        dat.binnedSpikes = [];
        dat.name = [];
        dat.frameON = [];
        dat.trialType = [];
        dat.condition = [];
        dat.comments = [];
    end

    datAllSubs = [dat; datAllSubs]; %concatenate across subsessions
    
end
    
% concatenate all trials across all subsessions 
% dimensions are trials x bins x cells
allSpikesAllSubs.binnedSpikes = (vertcat(datAllSubs(:).binnedSpikes));   % divide by win to convert to FR
allSpikesAllSubs.name = datAllSubs.name;        

% behavioral / qualitative info (same for all cells) 
allTrials.frameON = vertcat(datAllSubs(:).frameON);
allTrials.trialType = vertcat(datAllSubs(:).trialType);
allTrials.condition = vertcat(datAllSubs(:).condition);
allTrials.comments = vertcat(datAllSubs(:).comments);

% sum spike counts over bins (now ntrials x ncells) 
allCellsPerTrialSpikeCount = squeeze(sum(allSpikesAllSubs.binnedSpikes,2));

% calc Fano Factors 
bin.means = mean(allCellsPerTrialSpikeCount,1); % avg # of spikes/trial; NOT RATE
bin.SDs = std(allCellsPerTrialSpikeCount,[],1);
bin.meanFRs = mean(allCellsPerTrialSpikeCount,1)/(tmin+tmax); % now in sp/s

%%%%% do some clean up 
if strcmpi(subject, 'Thor')
    [allSpikesAllSubs, allTrials] = getRidOfChewTrials(allSpikesAllSubs, allTrials);
end

% convert spike-labels to strings 
 spikeLabels = (allSpikesAllSubs(:).name);
for i= 1:length(spikeLabels) 
    spikeLabels2plot(i,:) = str2double(strsplit(spikeLabels(i),'-'));
end 

% remove low FR neurons 
if threshlowFRFlag
    allCellsPerTrialSpikeCount = allCellsPerTrialSpikeCount(:,(bin.meanFRs >= threshlowFR));
    spikeLabels = spikeLabels(bin.meanFRs >= threshlowFR,:);
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


% save results 
if saveFlag 
    save('Fig3A_All.mat','allCellsPerTrialSpikeCount', 'spikeLabels2plot','cortRegion', 'allTrials','chls','bhvs2plot','ax');
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

