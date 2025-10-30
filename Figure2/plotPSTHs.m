
%%%%%%%% PLOTTING Figure 2C
%%%%%%%%  Single-Cell Activity and Selectivity in Cortical Face-Motor Regions
%%%%%%%%
%%%%%%%% Peri-event time histograms of four cells in the cortical face motor system 
%%%%%%%%  aligned to onset of three facial gestures (color-coded). 
%%%%%%%%  Heavy traces indicate means, shaded bars are 2 SEMs of trial-averaged response
%%%%%%%%  
%%%%%%%%  Input data includes 'psthCell', a cell array, each of nTrials x nTimepoints 
%%%%%%%%  Input data indludes 'trialTypescell', a cell array of nTrials (categorical gesture) 
%%%%%%%%% GRI 09/11/2025 

close all;
clear all;
set(0,'defaultAxesFontSize',36);
set(0,'defaultAxesFontWeight','bold');

% --- Load combined spikes file
S = load('matfiles/Fig2C.mat');
nFiles = numel(S.bins2takePerFile);  

win = 0.001; % in sec
frThresh = 0.1; % in hz

tmin = 1; % window to take 
tmax = 1;

minRestFlag = 1;
minRest = abs(tmin); % in sec; minimal rest prior to move onset (trials to include)

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

% loop each cell's PSTH 
for fileIdx = 1:nFiles

    % Indices of combined slices that came from this file
    cellIdx = find(S.sourceFileIdx == fileIdx);

    % nTrials that cell saw 
    nTrials_i = max(S.trialsPerCell(cellIdx));
    nBins_i   = numel(S.bins2takePerFile{fileIdx});

    % pool PSTH 
    allCellsPSTH = S.combinedPSTH(1:nTrials_i, 1:nBins_i, cellIdx);

    % Trial types for this cell 
    allTrials = struct();
    allTrials.trialType = S.trialTypesCell{cellIdx(1)};  % aligned to rows of allCellsPSTH

    % Time bins used in this file (indices into the original bin grid)
    bins2take = S.bins2takePerFile{fileIdx};

    % Spike labels 
    labels_this = S.combinedLabels(cellIdx, :);  % [nCells x 2] = [channel unit]
    spikeLabels = arrayfun(@(r) sprintf('%d-%d', labels_this(r,1), labels_this(r,2)), ...
        (1:size(labels_this,1)).', 'uni', false);
    allSpikesAllSubs = struct('name', spikeLabels(:));

    %  average FRs; this is the trial-averaged response PER cell
    %  allCellsTrialAvgResponse.bhv is nTimepoints x nCells (averaged over
    %       nReps)
    bhvs2plot = [1 2 4];
    for i = 1:length(bhvs2plot)
        allCellsTrialAvgResponse(i).bhv = squeeze(mean(allCellsPSTH((allTrials.trialType == bhvs2plot(i)),:,:)));
        allCellsTrialAvgResponse(i).nReps = sum(allTrials.trialType == bhvs2plot(i));
    end

    % mean FRs per cell, per bhv
    for bhv = 1:length(bhvs2plot)
        meanFRs(:,bhv) = mean(allCellsTrialAvgResponse(bhv).bhv); %,1
    end

    taxis2take = edges(bins2take);
    colors = [1 0 0; 0 1 0; 0 0 1];

    for unit =1:size(allCellsPSTH,3) % per cell
        if sum(meanFRs(unit,:),2) <= frThresh
            continue
        end

        % get GPI 
        [faceExpPref,depthPref] = calcFacePrefIndex(meanFRs(unit,:), length(bhvs2plot));
        
        % plot 
        fig1 = figure('units','normalized','outerposition',[0.3783    0.0835    0.6217    0.7974]);

        for bhv = 1: length(bhvs2plot) % per bhv

            thisBhvTrials = allTrials.trialType == bhvs2plot(bhv);

            fig1 = shadedErrorBar(taxis2take,(allCellsTrialAvgResponse(bhv).bhv(:,unit)),...
                {@(m) (mean(allCellsPSTH(thisBhvTrials,:,unit))),...
                @(x) 2*(std(allCellsPSTH(thisBhvTrials,:,unit))) ./ sqrt(sum(thisBhvTrials))},'lineprops',{'Color',colors(bhv,:),'linew',8});

            hold on
            set(gca,'box','off','tickdir','out');
            grid off
        end
        xline(0,'--','k','linew',5)
        hold off

        axes=gca;
        axes.FontWeight = 'bold'; Taxes.FontSize = 34;
        xlabel('Time, s');
        ylabel('Firing Rate, Hz');
        %legend({'Thr','LS','Chew'},'FontSize',17);

        title(['Unit-' spikeLabels{unit}]);
        subtitle(['GPI = ' num2str(round(faceExpPref,2))],'FontSize',26);


    end
end

function [faceExpPref,depthPref] = calcFacePrefIndex(Ris, nExp)

[Rmaxs, ~] = max(Ris,[],2);
[Rmins, ~] = min(Ris,[],2);

faceExpPref = (nExp - (sum(Ris) ./ Rmaxs)) ./ (nExp -1);
depthPref = (Rmaxs - Rmins ) ./ (Rmaxs + Rmins);


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

