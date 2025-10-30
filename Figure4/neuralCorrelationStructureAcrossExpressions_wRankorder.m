% written GI
%%%
%%%%%%%%  analyzing kinematic decoding and non-preserved neural correlations
%%%%%%%%   for Figure 4C, D

% PUPRPOSE: compute the 1) correlation between all neuron pairs over time (trial-averaged)
%%%                            cluster & sort neurons based on
%%%                            similarity of correlational values
%%%                     2) compute the rank order of each matrix (upper
%%%                             triangle) 
%%%                     3) save the sorted and unsorted correlation
%%%                             matrices for plotting 
%%%%%%%%% GRI 09/11/2025 


clear all; close all;
set(0,'defaultAxesFontSize', 24); % bc im blind
set(0,'defaultAxesFontWeight', 'bold'); % bc im blind

date = '171128';
subject ='Thor';
workdir = (['/Users/geena/Dropbox/PhD/SUAinfo/' subject '_' date '/Data4Analysis']);
outdir = (['/Users/geena/Dropbox/PhD/SUAinfo/' subject '_' date '/CorrMatricesV2']);

subsessions2plot = [ 1 2 4 5];
bhvs2plot = [1 2 4];
bhvStrings = {'Threat','Lipsmack','Chew'};
nClust = 5;
defBhv = 2;
nPerm = 10000; % permutations for significant testing; random relabeling of gestures
nBoot = 5000; % bootstrapping for CIs; redraw neurons and recomputed correlation matrices 

saveFlag = false;

win = 0.02; % in sec
tmin = 0.5; tmax = 1;
taxis = -tmin:win:tmax;
minRestFlag = true;
minRest = abs(tmin); % in sec; minimal rest prior to move onset (trials to include)
gaussKernStd = 50; % ms

threshlowFRFlag = 1; % 0/1;  remove low FR neurons
threshlowFR = 0.5; %  sp/s threshold

if ~exist(outdir, 'dir')
    mkdir(outdir)
end

% return longer taxis to remove edge artifacts
[gaussKern, tminLong, tmaxLong, taxisLong, taxis2take] = fixTimeAxisEdgeArtifacts(gaussKernStd, win, tmin, tmax);

% return the channel mapping
[regions] = getChannel2CorticalRegionMapping(subject, 1);

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

% do separately for all regions
for aa = 1:length(regions)

    chls = regions{1,aa}.channels;
    region = regions{1,aa}.label;

    if strcmp(region,'All')
        nClust = 10;
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
                            dat.binnedSpikes(:,:, counter) =   binSpikes( obhv.FT2dat.time(onsets), spTimes, win, tminLong, tmaxLong    ); %output is trials x bins
                        elseif strcmpi(subject,'Barney') || strcmpi(subject,'Dexter')
                            dat.binnedSpikes(:,:, counter) =   binSpikes( obhv.FT1dat.time(onsets), spTimes, win, tminLong, tmaxLong); %output is trials x bins
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

    %%%%%%%%%%%%%

    % clean it up
    allCellsPerTrialSpikeCount = squeeze(sum(allSpikesAllSubs.binnedSpikes,2));
    meanFRs = mean(allCellsPerTrialSpikeCount,1)/(tmax+tmin); % now in sp/s

    % remove NaN channels
    allSpikesAllSubs.binnedSpikes = allSpikesAllSubs.binnedSpikes(:,:,meanFRs ~= 0);
    allSpikesAllSubs.name = allSpikesAllSubs.name(meanFRs ~=0);
    meanFRs = meanFRs(meanFRs ~=0);

    % remove low FR neurons
    if threshlowFRFlag
        allSpikesAllSubs.binnedSpikes = allSpikesAllSubs.binnedSpikes(:,:,(meanFRs >= threshlowFR));
        allSpikesAllSubs.name = allSpikesAllSubs.name(meanFRs>= threshlowFR);
    end

    % dimensions are ntrials x nbins x ncells
    allCellsPSTH = zeros(size(allSpikesAllSubs.binnedSpikes,1), length(taxis), size(allSpikesAllSubs.binnedSpikes,3));

    for trial = 1:size(allSpikesAllSubs.binnedSpikes,1) %per trial
        for cell = 1:size(allSpikesAllSubs.binnedSpikes,3) % per cell
            tmp = conv(allSpikesAllSubs.binnedSpikes(trial, :, cell), gaussKern, 'same') ./win; %divide by win for FR
            allCellsPSTH(trial, :, cell) = tmp(:,taxis2take(1:end-1));
        end
    end

    % compute correlation between all neuron pairs over time (trial-avg exp)
    for bhv = 1:length(bhvs2plot)

        theseTrials = (allTrials.trialType == bhvs2plot(bhv));
        trialAvg(bhv).data = squeeze(mean(allCellsPSTH(theseTrials,:,:))); % per bhv, trial-average (nBins x nCells)
        signalCorrOverTime(bhv).corr = corrcoef(trialAvg(bhv).data);
        signalCorrOverTime(bhv).corr(isnan(signalCorrOverTime(bhv).corr)) = 0;

        % to collect pairwise correlations, take triangular upper matrix & remove
        % diagonal elements, then transform into a row vector with squareform()
        triuBhv = triu(signalCorrOverTime(bhv).corr);
        D = diag(triuBhv);
        prCorr{:,bhv} = [squareform((triuBhv-diag(D)).')];

        tiAlpha(bhv).data = ones(size(signalCorrOverTime(bhv).corr));
        tiAlpha(bhv).data(isnan(signalCorrOverTime(bhv).corr))=0;

    end

    % cluster & sort, according to ONE ordering (unified sorting)
    corrMatDef = signalCorrOverTime(defBhv).corr - eye(size(signalCorrOverTime(defBhv).corr));
    dissimSignalDef = 1 - corrMatDef';

    % compute distances
    ZsigDef = linkage(dissimSignalDef,'complete');

    % cluster the data
    TsigDef = cluster(ZsigDef,'maxclust',nClust);

    % get new sorted indices of defining bhv
    [~, defSortIndxSig] = sort(TsigDef);

    % stash the sorted, self-defining bhv correlation matrix
    signalCorrOverTime(defBhv).unifiedSort.data = dissimSignalDef(defSortIndxSig,defSortIndxSig);
    signalCorrOverTime(defBhv).selfSort.data = dissimSignalDef(defSortIndxSig,defSortIndxSig);

    % Build aligned correlation matrices (same neuron order for all behaviors)
    R1 = signalCorrOverTime(1).corr(defSortIndxSig, defSortIndxSig);
    R2 = signalCorrOverTime(2).corr(defSortIndxSig, defSortIndxSig);
    R3 = signalCorrOverTime(3).corr(defSortIndxSig, defSortIndxSig);

    % (Optional but recommended) zero the diagonal to avoid any stray NaNs
    R1(1:size(R1,1)+1:end) = 0;
    R2(1:size(R2,1)+1:end) = 0;
    R3(1:size(R3,1)+1:end) = 0;

    % plot for sanity
    figure; imagesc(R1'); colormap(bluewhitered)
    figure; imagesc(R2'); colormap(bluewhitered)
    figure; imagesc(R3'); colormap(bluewhitered)

    % Run rank-order similarity (Spearman) with Mantel permutations + bootstrap CIs
    out = rankOrderSpearmanKendall3(R1, R2, R3, ...
        'NPerm', nPerm, 'NBoot', nBoot, 'Tail', 'both', 'Seed', 123, 'UseDiag', false);
    disp('=== Rank-order similarity across the three behaviors (signal) ===');
    disp(out.summary);

    % save the cluster IDs and spikeLabels
    signalCorrOverTime(defBhv).unifiedSort.clusterID = TsigDef;
    signalCorrOverTime(defBhv).unifiedSort.sortIndx = defSortIndxSig;
    signalCorrOverTime(defBhv).unifiedSort.spikeLabels = allSpikesAllSubs.name(defSortIndxSig);
    signalCorrOverTime(defBhv).selfSort.clusterID = TsigDef;
    signalCorrOverTime(defBhv).selfSort.sortIndx = defSortIndxSig;
    signalCorrOverTime(defBhv).selfSort.spikeLabels = allSpikesAllSubs.name(defSortIndxSig);

    % keep track of original indices
    spikeIdxOriginal{:,defBhv} = allSpikesAllSubs.name;

    % Now order remaining bhv correlation matrices by hierarchical clusters
    % defined by defBhv
    nonDefBhvs = setxor(defBhv, 1:length(bhvs2plot));

    for bhv = nonDefBhvs

        % remove diagonal elements & convert to dissimilarity matrix
        corrMat = signalCorrOverTime(bhv).corr - eye(size(signalCorrOverTime(bhv).corr));
        dissimSignal = 1 - corrMat';

        % stash defining bhv sorted correlation matrices
        signalCorrOverTime(bhv).unifiedSort.data = dissimSignal(defSortIndxSig,defSortIndxSig);

        % compute distances
        Zsig = linkage(dissimSignal,'complete');

        % cluster the data
        Tsig = cluster(Zsig,'maxclust',nClust);

        % get self-defined sorted indices
        [~, sortIndxSig] = sort(Tsig);

        % stash the self-sorted correlation matrices
        signalCorrOverTime(bhv).selfSort.data = dissimSignal(sortIndxSig,sortIndxSig);

        % save the signalCorr cluster IDs and spike Labels
        signalCorrOverTime(bhv).unifiedSort.clusterID = TsigDef; % def-bhv sorted clusters
        signalCorrOverTime(bhv).unifiedSort.sortIndx = defSortIndxSig;
        signalCorrOverTime(bhv).unifiedSort.spikeLabels = allSpikesAllSubs.name(defSortIndxSig);

        signalCorrOverTime(bhv).selfSort.clusterID = Tsig; % self-sorted clusters
        signalCorrOverTime(bhv).selfSort.sortIndx = sortIndxSig;
        signalCorrOverTime(bhv).selfSort.spikeLabels = allSpikesAllSubs.name(sortIndxSig);

        % keep track of original indices
        spikeIdxOriginal{:,bhv} = allSpikesAllSubs.name;
    end

    if saveFlag
        if ~exist(outdir, 'dir')
            mkdir(outdir)
        end

        save([outdir '/' regions{aa}.label], 'signalCorrOverTime','noiseCorr','trialAvg',...
            'region','allTrials','corrMatDef','defBhv','date','subject',...
            'meanFRs','nClust','prCorr','prNoise','subsessions2plot','taxis',...
            'threshlowFR','threshlowFRFlag','allCellsPSTH','spikeIdxOriginal','win');

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







