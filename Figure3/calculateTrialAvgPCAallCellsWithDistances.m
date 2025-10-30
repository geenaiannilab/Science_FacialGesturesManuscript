% written GI updated 250817
%
% PUPRPOSE: Calculate & plot reduced neural "state space"
%   where each point in reduced neural state-space represents a BIN in TIME
%   
% Input to PCA is allCellsTrialAvgResponse, which is nBins x nCells, where
% each entry is the smoothed FR at that timepoint for that cell 
%
% This script will calculate the 1) loadings of each neuron onto each PC 
%                                2) the projections of the data onto the
%                                PCs
%                                3) variance explained by each PC 
%                                4) distance between neural trajectories
%                                over time
% 


clear all; close all;
set(0,'defaultAxesFontSize',36); set(0, 'DefaultLineLineWidth', 4);
set(0,'defaultAxesFontWeight','bold')

date = '210704';
subject ='Barney';
workdir = (['/Users/geena/Dropbox/PhD/SUAinfo/' subject '_' date '/Data4Analysis']);
saveFlag = true; 

regions = getChannel2CorticalRegionMapping(subject, 1);
chls = 1:240;

subsessions2plot = [2:5 7:10];
bhvs2plot = [1 2 4];

% DO NOT CHANGE!
win = 0.02;  % in sec
tmin = 1;
tmax = 1;
minRestFlag = 1;
minRest = abs(tmin); % in sec; minimal rest prior to move onset (trials to include)

centerOnlyFlag = 0; % ONLY mean center input to PCA (dont divide by std)
centerNormalizeFlag = 1; % normalize cells/variables by their variance in addition to mean-centering (PCA on correlation matrix)    
threshlowFRFlag = 1; % 0/1;  remove low FR neurons 
threshlowFR = 0.1; %  sp/s threshold 
%%% NOTE RE: PREPROCESSING INPUT
%%% pca() in matlab does mean-centering *by default*, but does not divide by stdev by default 

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

restSubs = subs2show(contains(subs2show,'rest'));

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

if strcmpi(subject, 'Thor')
    [allSpikesAllSubs, allTrials] = getRidOfChewTrials(allSpikesAllSubs, allTrials);
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

% mean FRs per cell, per bhv
for bhv = 1:length(bhvs2plot)
    meanFRs(:,bhv) = mean(allCellsTrialAvgResponse(bhv).bhv,1);
end

% convert spike-labels to strings for plotting later
spikeLabels = (allSpikesAllSubs(:).name);
for i= 1:length(spikeLabels) 
    spikeLabels2plot(i,:) = str2double(strsplit(spikeLabels(i),'-'));
end 

% concatente all averaged responses for PCA input 
neuralInput = vertcat(allCellsTrialAvgResponse(:).bhv);

% remove low FR neurons 
if threshlowFRFlag
    thrCells2use = ~sum(meanFRs >= threshlowFR,2) == 0;
    neuralInput = neuralInput(:,thrCells2use);
    spikeLabels = spikeLabels(thrCells2use,:);
    spikeLabels2plot = spikeLabels2plot(thrCells2use,:);
end
   
% double check you have enough cells
 if size(neuralInput,2) <= 5
    disp(['fewer than 5 cells, (' num2str(sum(cells2use)) ') quitting'])
    return
 end


%% ============================
%% Trial-Averaged PCA 
% Input is nBins x nCells,
%   where each entry is the trial-averaged smoothed FR within a time-bin
% Each point in state-space is a TIMEPOINT
%% ============================

for rr = 1:length(regions)

    % chls in this array
    chls2take = ismember(spikeLabels2plot(:,1), regions{rr}.channels);

    pcaInput = neuralInput(:,chls2take);

    if centerOnlyFlag
        [coeff,score,latent, ~] = pca(pcaInput); % centering (subtract column/cell means) done automatically by pca()
        input2plot = pcaInput - mean(pcaInput); %for plotting only, later
    elseif centerNormalizeFlag
        [coeff,score,latent, ~] = pca(pcaInput,'VariableWeights', 'variance');
        input2plot = (pcaInput - mean(pcaInput)) ./ std(pcaInput,[],1);
    else
        [coeff,score,latent, ~] = pca(pcaInput,'Centered','off');
        input2plot = pcaInput;
        disp('performing PCA on uncentered, non-normalized data!')
    end

    % store in structure

    allDataOut.(regions{rr}.label).coeff = coeff;
    allDataOut.(regions{rr}.label).score = score;
    allDataOut.(regions{rr}.label).latent = latent;
    allDataOut.(regions{rr}.label).explained = latent/sum(latent);
    allDataOut.(regions{rr}.label).nDimGreater90 = find(round(cumsum(allDataOut.(regions{rr}.label).explained),2) >= 0.9, 1);

    % pull out per bhv neural Traj
    start = 1;
    for bhv = 1:length(bhvs2plot)
        stop = length(taxis2take)*bhv;
        allDataOut.(regions{rr}.label).bhv(bhv).neuralTraj = score(start:stop,:);
        start = start + length(taxis2take);
    end

end


%% ============================
%% calculate euclidean neural space starting positions of each trajectory
%% ============================

for rr = 1:length(regions)

    maxDim = allDataOut.(regions{rr}.label).nDimGreater90;

    % position coordinates
    ThrCoords = allDataOut.(regions{rr}.label).bhv(1).neuralTraj(:,1:maxDim);
    LSCoords = allDataOut.(regions{rr}.label).bhv(2).neuralTraj(:,1:maxDim);
    ChewCoords = allDataOut.(regions{rr}.label).bhv(3).neuralTraj(:,1:maxDim);

    % euclidean distance as a function of time, from dim =1 to max dimensions
    for dd = 1:maxDim

        % nNeuralDimensional distance vector, per bhv
        sosThrChDistance = 0; sosThrLSDistance = 0; sosLSChDistance = 0 ;

        for dim = 1:dd
            % nNeuralDimensional distance vector, per bhv
            sosThrChDistance = sosThrChDistance + ((ThrCoords(:,dim) - ChewCoords(:,dim)) .^2);
            sosThrLSDistance = sosThrLSDistance + ((ThrCoords(:,dim) - LSCoords(:,dim)) .^2);
            sosLSChDistance = sosLSChDistance + ((LSCoords(:,dim) - ChewCoords(:,dim)) .^2);
        end

        % nNeuralDimensional distance vector, per bhv
        allDataOut.(regions{rr}.label).distanceThrVCh(:,dd) = sqrt(sosThrChDistance);
        allDataOut.(regions{rr}.label).distanceThrVLS(:,dd) = sqrt(sosThrLSDistance);
        allDataOut.(regions{rr}.label).distanceLSVCh(:,dd) = sqrt(sosLSChDistance);
        allDataOut.(regions{rr}.label).distanceAverage(:,dd) = mean([allDataOut.(regions{rr}.label).distanceThrVCh(:,dd), allDataOut.(regions{rr}.label).distanceThrVLS(:,dd),allDataOut.(regions{rr}.label).distanceLSVCh(:,dd)],2);

    end % end max dimensionality

end

%% ============================
%% SHUFFLED BASED SIGNIFICANCE TESTING 
%% ============================

% Example: significance mask for M1, Threat vs Chew
% shuffleStats.M1.ThrVCh.sig_fdr   -> [T x 1] logical
% shuffleStats.M1.ThrVCh.p_fdr     -> FDR p-values
% shuffleStats.M1.ThrVCh.z         -> z-score vs null

cfg = struct();
cfg.nShuf       = 1000;      % number of label shuffles
cfg.alphaFDR    = 0.05;      % FDR target per region×pair
cfg.useMaxDim   = true;      
cfg.dimIdx      = [];        
cfg.rngSeed     = 13;        % for reproducibility

[shuffleStats, shuffleNull] = test_distance_shuffle( ...
    allDataOut, allCellsPSTH, allTrials, ...
    regions, bhvs2plot, spikeLabels2plot, bins2take, taxis2take, ...
    cfg);

if saveFlag
    save([ workdir '/neuralTraj_' regions{rr}.label '_0.2.mat'], 'allDataOut','chls',...
        'date','subject','subsessions2plot','bhvs2plot','win','tmin','tmax','taxis2take',...
        'shuffleStats','shuffleNull','cfg');
end


%% ============================
%% subfunctions for significance testing  
%% ============================

function [stats, nullStore] = test_distance_shuffle( ...
    allDataOut_obs, allPSTH, allTrials, ...
    regions, bhvs2plot, spikeLabels2plot, bins2take, taxis2take, cfg)
% Label-shuffle permutation test for trajectory distances over time.
% Uses precomputed PSTHs and your exact PCA+distance pipeline.
%
% Inputs:
%   allDataOut_obs : observed struct from your code (contains distance* fields)
%   allPSTH        : [nTrials x nBins x nCells] smoothed FRs (your allCellsPSTH)
%   allTrials      : struct with .trialType (vector over trials)
%   regions        : your regions cell array (each has .label and .channels)
%   bhvs2plot      : 1x3 vector of behavior codes (order must match your code)
%   spikeLabels2plot : [nCells x 2] numeric [chan, unit] per cell
%   bins2take      : vector of time bin indices used for trajectories
%   taxis2take     : time axis (for length T)
%   cfg            : struct with fields:
%       .nShuf, .alphaFDR, .useMaxDim, .dimIdx, .rngSeed
%
% Outputs:
%   stats.(regionLabel).(pairName) with fields:
%       .obs   [T x 1] observed distance curve (chosen dim)
%       .p     [T x 1] one-sided p (P(null >= obs))
%       .p_fdr [T x 1] BH-FDR corrected p
%       .sig_fdr [T x 1] logical mask (p_fdr < alpha)
%       .z     [T x 1] z-score vs null
%   nullStore.(regionLabel).(pairName) : [T x nShuf] null distances
%
% Notes:
% - Pair names: 'ThrVCh', 'ThrVLS', 'LSVCh' (and 'Average' kept separately)
% - One-sided test (greater-than): conservative finite-sample estimate.
% - FDR applied per region×pair across timepoints.

if ~isfield(cfg,'nShuf'),     cfg.nShuf = 1000; end
if ~isfield(cfg,'alphaFDR'),  cfg.alphaFDR = 0.05; end
if ~isfield(cfg,'useMaxDim'), cfg.useMaxDim = true; end
if ~isfield(cfg,'dimIdx'),    cfg.dimIdx = []; end
if ~isfield(cfg,'rngSeed'),   cfg.rngSeed = 13; end

rng(cfg.rngSeed);

% Decide which PCA dimension to test for each region
dimPerRegion = containers.Map();
for rr = 1:numel(regions)
    rlab = regions{rr}.label;
    if cfg.useMaxDim
        dimPerRegion(rlab) = allDataOut_obs.(rlab).nDimGreater90;
    else
        if isempty(cfg.dimIdx)
            error('cfg.dimIdx must be provided when cfg.useMaxDim == false.');
        end
        dimPerRegion(rlab) = cfg.dimIdx;
    end
end

% Determine T from observed curves
T = numel(taxis2take);

pairs = {'ThrVCh','ThrVLS','LSVCh'};  % keep Average separately
stats = struct();
nullStore = struct();

% ---------- OBSERVED: pick the chosen dimension column ----------
for rr = 1:numel(regions)
    rlab = regions{rr}.label;
    dd   = dimPerRegion(rlab);

    obs.ThrVCh = allDataOut_obs.(rlab).distanceThrVCh(:,dd);
    obs.ThrVLS = allDataOut_obs.(rlab).distanceThrVLS(:,dd);
    obs.LSVCh  = allDataOut_obs.(rlab).distanceLSVCh(:,dd);
    obs.Average= allDataOut_obs.(rlab).distanceAverage(:,dd); % optional

    stats.(rlab) = struct();
    for p = 1:numel(pairs)
        stats.(rlab).(pairs{p}).obs = obs.(pairs{p});
    end
    stats.(rlab).Average.obs = obs.Average; % keep for plotting/reference
end

% ---------- NULL: build via label shuffles ----------
% Preallocate
for rr = 1:numel(regions)
    rlab = regions{rr}.label;
    for p = 1:numel(pairs)
        nullStore.(rlab).(pairs{p}) = nan(T, cfg.nShuf);
    end
    nullStore.(rlab).Average = nan(T, cfg.nShuf);
end

trialType_obs = allTrials.trialType(:);
nTrials = numel(trialType_obs);

for s = 1:cfg.nShuf
    % permute labels across trials (preserves per-trial FRs & trial count)
    permIdx = randperm(nTrials);
    trialType_shuf = trialType_obs(permIdx);

    % recompute trajectories, PCA, and distances for this shuffle
    allDataOut_shuf = compute_distances_for_trialtype( ...
        allPSTH, trialType_shuf, bhvs2plot, spikeLabels2plot, ...
        regions, bins2take);

    % store chosen dim curves
    for rr = 1:numel(regions)
        rlab = regions{rr}.label;
        dd   = min(dimPerRegion(rlab), size(allDataOut_shuf.(rlab).distanceAverage,2));

        nullStore.(rlab).ThrVCh(:,s) = allDataOut_shuf.(rlab).distanceThrVCh(:,dd);
        nullStore.(rlab).ThrVLS(:,s) = allDataOut_shuf.(rlab).distanceThrVLS(:,dd);
        nullStore.(rlab).LSVCh(:,s)  = allDataOut_shuf.(rlab).distanceLSVCh(:,dd);
        nullStore.(rlab).Average(:,s)= allDataOut_shuf.(rlab).distanceAverage(:,dd);
    end
end

% ---------- Stats: one-sided p, FDR over time, z-score ----------
for rr = 1:numel(regions)
    rlab = regions{rr}.label;

    for p = 1:numel(pairs)
        obsCurve  = stats.(rlab).(pairs{p}).obs;
        nullMat   = nullStore.(rlab).(pairs{p}); % [T x nShuf]

        % one-sided p = P(null >= obs)
        geCounts = sum(bsxfun(@ge, nullMat, obsCurve), 2);
        p_one = (geCounts + 1) / (size(nullMat,2) + 1);

        % BH-FDR across T time bins
        p_fdr = bh_fdr(p_one, cfg.alphaFDR);

        % z-score vs null
        mu  = mean(nullMat, 2, 'omitnan');
        sig = std(nullMat, 0, 2, 'omitnan');
        z   = (obsCurve - mu) ./ sig;

        stats.(rlab).(pairs{p}).p      = p_one;
        stats.(rlab).(pairs{p}).p_fdr  = p_fdr;
        stats.(rlab).(pairs{p}).sig_fdr= p_fdr < cfg.alphaFDR;
        stats.(rlab).(pairs{p}).z      = z;
    end

    % Keep Average stats too (optional; not typically hypothesis-tested)
    obsCurve  = stats.(rlab).Average.obs;
    nullMat   = nullStore.(rlab).Average;
    geCounts  = sum(bsxfun(@ge, nullMat, obsCurve), 2);
    p_one     = (geCounts + 1) / (size(nullMat,2) + 1);
    p_fdr     = bh_fdr(p_one, cfg.alphaFDR);
    mu        = mean(nullMat, 2, 'omitnan');
    sig       = std(nullMat, 0, 2, 'omitnan');
    z         = (obsCurve - mu) ./ sig;
    stats.(rlab).Average.p       = p_one;
    stats.(rlab).Average.p_fdr   = p_fdr;
    stats.(rlab).Average.sig_fdr = p_fdr < cfg.alphaFDR;
    stats.(rlab).Average.z       = z;
end
end

function allDataOut = compute_distances_for_trialtype( ...
    allPSTH, trialType_vec, bhvs2plot, spikeLabels2plot, regions, bins2take)

% Build trial-averaged responses PER behavior (nTime x nCells)
allCellsTrialAvgResponse = struct([]);
for i = 1:numel(bhvs2plot)
    sel = (trialType_vec == bhvs2plot(i));
    % average across trials, keep only bins2take
    allCellsTrialAvgResponse(i).bhv = squeeze(mean(allPSTH(sel, bins2take, :), 1));
    allCellsTrialAvgResponse(i).nReps = sum(sel);
end

% Concatenate behaviors -> PCA input (time x cells)
neuralInput = vertcat(allCellsTrialAvgResponse(:).bhv);

% ---- Per region: PCA on its channels, then distances over time ----
for rr = 1:numel(regions)
    rlab = regions{rr}.label;
    % pick cells in this region by channel match
    chls2take = ismember(spikeLabels2plot(:,1), regions{rr}.channels);
    pcaInput = neuralInput(:, chls2take);

    % compute PCA 
    [coeff, score, latent] = pca(pcaInput); 

    explained = latent / sum(latent);
    nDimGreater90 = find(round(cumsum(explained), 2) >= 0.9, 1);
    if isempty(nDimGreater90), nDimGreater90 = min(3, size(score,2)); end

    % pull out trajectories per behavior
    T  = size(allCellsTrialAvgResponse(1).bhv,1);
    s0 = 0;
    for b = 1:numel(bhvs2plot)
        s1 = s0 + T;
        allDataOut.(rlab).bhv(b).neuralTraj = score(s0+1:s1, :);
        s0 = s1;
    end

    allDataOut.(rlab).coeff = coeff;
    allDataOut.(rlab).score = score;
    allDataOut.(rlab).latent = latent;
    allDataOut.(rlab).explained = explained;
    allDataOut.(rlab).nDimGreater90 = nDimGreater90;

    % Euclidean distances over time up to dd dims
    maxDim = nDimGreater90;
    ThrCoords  = allDataOut.(rlab).bhv(1).neuralTraj(:, 1:maxDim);
    LSCoords   = allDataOut.(rlab).bhv(2).neuralTraj(:, 1:maxDim);
    ChewCoords = allDataOut.(rlab).bhv(3).neuralTraj(:, 1:maxDim);

    % prep
    Tpts = size(ThrCoords,1);
    allDataOut.(rlab).distanceThrVCh  = nan(Tpts, maxDim);
    allDataOut.(rlab).distanceThrVLS  = nan(Tpts, maxDim);
    allDataOut.(rlab).distanceLSVCh   = nan(Tpts, maxDim);
    allDataOut.(rlab).distanceAverage = nan(Tpts, maxDim);

    for dd = 1:maxDim
        d1 = sqrt(sum((ThrCoords(:,1:dd)  - ChewCoords(:,1:dd)).^2, 2));
        d2 = sqrt(sum((ThrCoords(:,1:dd)  - LSCoords(:,1:dd)).^2, 2));
        d3 = sqrt(sum((LSCoords(:,1:dd)   - ChewCoords(:,1:dd)).^2, 2));

        allDataOut.(rlab).distanceThrVCh(:,dd)  = d1;
        allDataOut.(rlab).distanceThrVLS(:,dd)  = d2;
        allDataOut.(rlab).distanceLSVCh(:,dd)   = d3;
        allDataOut.(rlab).distanceAverage(:,dd) = mean([d1 d2 d3], 2);
    end
end
end


% --------- helpers ---------

function fill_between(x, y1, y2, faceColor)
% Shade the area between y1 and y2 over x
x = x(:);
y1 = y1(:);
y2 = y2(:);
X = [x; flipud(x)];
Y = [y1; flipud(y2)];
patch('XData',X,'YData',Y, 'FaceColor',faceColor, 'EdgeColor','none', 'FaceAlpha',0.6);
end

function shade_sig_runs(tAxis, sigMask, ylims, baseColor)
% Shade contiguous significant time runs as a light background patch.
% sigMask: logical [T x 1]
if ~any(sigMask), return; end

% lighter version of baseColor
tint = 0.2 + 0.8*baseColor;  % pull toward white
tint(tint>1)=1;

% find contiguous runs
dSig = diff([false; sigMask(:); false]);
runStarts = find(dSig==1);
runEnds   = find(dSig==-1) - 1;

for k = 1:numel(runStarts)
    i0 = runStarts(k); i1 = runEnds(k);
    x0 = tAxis(i0);    x1 = tAxis(i1);
    patch('XData',[x0 x1 x1 x0], 'YData',[ylims(1) ylims(1) ylims(2) ylims(2)], ...
          'FaceColor', tint, 'FaceAlpha', 0.15, 'EdgeColor','none');
end

uistack(findobj(gca,'Type','patch'),'bottom'); % keep shading behind curves
end

function p_fdr = bh_fdr(p, alpha)
% Benjamini-Hochberg FDR across a vector p (length m).
% Returns the *adjusted* p-values; threshold with p_fdr < alpha.
m = numel(p);
[ps, idx] = sort(p(:));
th = (1:m)'/m * alpha;
rej = ps <= th;
k = find(rej, 1, 'last');
p_fdr = nan(size(p));
if isempty(k)
    % No rejections; conservative adjusted p-values
    adj = min(1, ps .* m ./ (1:m)');
else
    % Compute adjusted p-values monotonized
    adj = zeros(m,1);
    for i = 1:m
        adj(i) = min(1, min(ps(i:end) .* m ./ (i:m)'));
    end
end
p_fdr(idx) = adj;
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
