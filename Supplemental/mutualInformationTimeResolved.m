
%%%% PLOTTING Fig S10
%%%% Median bias-corrected mutual information time courses are shown per region, aligned to movement onset (0 sec). 
%%%% Shaded colored curves show median MIc across neurons (400 ms sliding window, 10 ms steps). 
%%%% Shaded areas indicate time intervals of significant information above zero, identified via cluster-based permutation testing 
%%%% (10k sign flips, cluster-forming threshold one sided p < 0.01, family-error controlled). 

%%%% Vertical dashed lines mark the earliest significant latency, defined as the first occurrence of at least two contiguous supra-threshold bins. 

% written GI 250831


%% ================== USER SETTINGS ==================
dates   = {'Thor_171005','Thor_171010','Thor_171027','Thor_171128',...
    'Barney_210704','Barney_210706','Barney_210805'};
baseDir = '~/Dropbox/PhD/SUAInfo/';
runSim = false; 
set(0,'defaultAxesFontSize',20)

% Time-resolved MI parameters
win_ms    = 400;           % window length (ms)
step_ms   = 10;            % step between windows (ms)
nShuffles = 200;           % gesture label shuffles per window
alpha_bin = 0.05;          % per-bin (cluster-forming) alpha (one-sided)
nPerm     = 10000;          % cluster permutations (sign-flips across neurons)
minRun    = 2;             % minimum consecutive bins to declare latency (e.g., 1)

baseColors = [0.494 0.184 0.556; 0.635 0.078 0.184; 0.85 0.325 0.098; 0.929 0.694 0.125];

%% ================== POOLS ACROSS DAYS ==================
MIcorr_days   = {};   % [nWin x nCells_day]
MIraw_days    = {};   % [nWin x nCells_day]
pvals_days    = {};   % [nWin x nCells_day]
regions_days  = {};   % 1 x nCells_day (cellstr)
centers_days  = {};   % 1 x nWin (s)

% store best example neuron per region (across all days)
examples = struct('region',{},'day',{},'cellIdx',{},'spikes',{},'labels',{},'binSize',{});

%% ================== MAIN DAY LOOP ==================
for dd = 1:numel(dates)
    dayID  = dates{dd};
   
    if runSim 
        dat = simulate_timeResolved_MI_dataset();
    else
        workDir = [baseDir dayID '/Data4Analysis/'];
        dat = load([workDir 'binnedSpikesForMI.mat']); 
    end
    
    % Quick sanity check:
    s = squeeze(sum(dat.allSpikesAllSubs.binnedSpikes,2));  % trial x cells total counts

    spikes = dat.allSpikesAllSubs.binnedSpikes;    % [trials x time x cells]
    labels = dat.allTrials.trialType;
    regions_full = dat.cortRegion;

    % Bin size (s)
    if isfield(dat,'taxis2take') && numel(dat.taxis2take) > 1
        binSize = mode(diff(dat.taxis2take));
        t0      = dat.taxis2take(1);
    else
        binSize = 0.001; t0 = 0;
    end

    % Remove near-silent cells
    cellTotals = squeeze(sum(sum(spikes,1),2));
    keepCells  = cellTotals > 1;
    spikes     = spikes(:,:,keepCells);
    regions    = regions_full(keepCells);
    nCells_day = size(spikes,3);
    nTime      = size(spikes,2);

    % Sliding windows
    winBins  = max(1, round((win_ms/1000)/binSize));
    stepBins = max(1, round((step_ms/1000)/binSize));
    starts   = 1:stepBins:(nTime - winBins + 1);
    nWin     = numel(starts);
    centersS = t0 + ((starts - 1) + winBins/2) * binSize;  % seconds

    MI_raw_win  = nan(nWin, nCells_day);
    MI_corr_win = nan(nWin, nCells_day);
    pvals_win   = nan(nWin, nCells_day);

    % Time-resolved MI per window
    for w = 1:nWin
        idx = starts(w):(starts(w)+winBins-1);
        counts2D = squeeze(sum(spikes(:, idx, :), 2));  % [trials x cells]
        [MI_corr_w, MI_raw_w, pvals_w] = computeMI_perCell_fast(counts2D, labels, nShuffles);
        MI_raw_win(w,:)  = MI_raw_w;
        MI_corr_win(w,:) = MI_corr_w;
        pvals_win(w,:)   = pvals_w;
    end

    % Save per day
    MIcorr_days{dd}  = MI_corr_win;
    MIraw_days{dd}   = MI_raw_win;
    pvals_days{dd}   = pvals_win;
    regions_days{dd} = regions(:)';
    centers_days{dd} = centersS(:)';
end

%% Align time across days (trim to shortest)
nWin_each = cellfun(@numel, centers_days);
nWin_min  = min(nWin_each);
for dd = 1:numel(dates)
    MIcorr_days{dd}  = MIcorr_days{dd}(1:nWin_min, :);
    MIraw_days{dd}   = MIraw_days{dd}(1:nWin_min, :);
    pvals_days{dd}   = pvals_days{dd}(1:nWin_min, :);
    centers_days{dd} = centers_days{dd}(1:nWin_min);
end
timeVec = centers_days{1};
nWin    = nWin_min;

%% Pool per window across days
regions_all = [regions_days{:}];
uRegs = unique(regions_all, 'stable');
nRegs = numel(uRegs);
colors = (nRegs <= size(baseColors,1)) * baseColors(1:nRegs,:) + (nRegs > size(baseColors,1)) * lines(nRegs);

% build region-wise matrices: for each region, MIcorr_region: [nWin x nCellsInRegion]
MIcorr_reg = cell(nRegs,1);
for r = 1:nRegs
    MIc = [];
    for dd = 1:numel(dates)
        idx = strcmp(regions_days{dd}, uRegs{r});
        MIc = [MIc, MIcorr_days{dd}(:, idx)];
    end
    MIcorr_reg{r} = MIc;
end

%% -------- Cluster-based permutation per region (one-sided > 0) --------
% where does mean MI(t) across neurons exceed chance more consistently than would happen in sign-flipped null data?

alpha_cluster = alpha_bin;         % cluster-forming per-bin alpha
[minMask, latencies] = deal(cell(nRegs,1), nan(nRegs,1));
sigMask_reg = cell(nRegs,1);       % logical mask over time for each region

% do the permutation testing and determine latency by first two contiguous
% bins above threshold 
for r = 1:nRegs
    M = MIcorr_reg{r};                   % [time x cells]
    [sigMask, ~, ~] = cluster_perm_time(M, nPerm, alpha_cluster); % returns 1 x nWin logical
    sigMask_reg{r} = sigMask(:)';
    
    % Latency: first index with >= minRun consecutive true
    latIdx = first_run(sigMask, minRun);
    if ~isnan(latIdx), latencies(r) = timeVec(latIdx); end
end

%% -------- Plot per-region MI(t) with shaded significance + latency and example PSTH --------
fig = figure; 
for r = 1:nRegs
    M = MIcorr_reg{r};                       % [time x cells]
    medMI = median(M,2,'omitnan');
    
    subplot(2,2,r)
    plot(timeVec, medMI, 'LineWidth', 2, 'Color', colors(r,:));

    % shade sig regions
    mask = sigMask_reg{r};
    shade_sig_regions(timeVec, medMI, mask, colors(r,:));
    
    % latency line
    if ~isnan(latencies(r))
        xline(latencies(r), '--k', 'LineWidth', 1.5);
        text(latencies(r), max(medMI)*1.02, sprintf('Latency = %.0f ms', latencies(r)*1000), ...
            'HorizontalAlignment','left','VerticalAlignment','bottom','FontSize',16);
    end
    xlabel('Time (s)'); ylabel('Median MI_{corr} (bits)');
    xlim([timeVec(1) timeVec(end)])
    grid on;
end
sgtitle(['All Regions: ' num2str(win_ms) 'ms window, 10ms step: cluster-perm p < ', num2str(alpha_cluster)]);

%% ================== HELPERS ==================

function [MI_corr, MI_raw, pvals] = computeMI_perCell_fast(spikeCounts, labels, nShuffles)
    labels = labels(:);
    [~,~,labBins] = unique(labels); % 1..K
    nCells = size(spikeCounts,2);
    MI_raw = zeros(1,nCells);
    pvals  = zeros(1,nCells);

    for c = 1:nCells
        MI_raw(c) = mi_discrete(spikeCounts(:,c), labBins);
    end

    shufMI = zeros(nShuffles,nCells);
    for s = 1:nShuffles
        shufLab = labBins(randperm(numel(labBins)));
        for c = 1:nCells
            shufMI(s,c) = mi_discrete(spikeCounts(:,c), shufLab);
        end
    end
    MI_corr = MI_raw - mean(shufMI,1);

    for c = 1:nCells
        pvals(c) = (1 + sum(shufMI(:,c) >= MI_raw(c))) / (1 + nShuffles);
    end
end

function mi = mi_discrete(resp, labBins)
    [respVals,~,rIdx] = unique(resp);
    K = max(labBins);
    R = numel(respVals);
    joint = accumarray([rIdx, labBins], 1, [R, K]);
    joint = joint / sum(joint(:));
    p_r = sum(joint,2);
    p_l = sum(joint,1);
    mi = 0;
    for i = 1:R
        for j = 1:K
            pij = joint(i,j);
            if pij > 0
                mi = mi + pij * log2( pij / (p_r(i) * p_l(j)) );
            end
        end
    end
end

function [sigMask, clusterStats, thr] = cluster_perm_time(M, nPerm, alpha_bin)
% M: [time x cells] bias-corrected MI per neuron. Null ~ 0.
% Sign-flip permutation across neurons preserves temporal structure.
% Returns sigMask over time (logical), cluster stats, and per-bin Z threshold used.

    [T,N] = size(M);
    if N < 2, sigMask = false(1,T); clusterStats=[]; thr=NaN; return; end

    % observed summary (mean across neurons)
    obs = mean(M,2,'omitnan');

    % permutations: sign-flip each neuron's entire timecourse
    permMeans = zeros(T, nPerm);
    signs = 2*(rand(N,nPerm)>0.5)-1;  % ±1 per neuron per perm
    for p = 1:nPerm
        Mp = M .* signs(:,p)';        % broadcast across time
        permMeans(:,p) = mean(Mp,2,'omitnan');
    end

    % Z-score using permutation distribution
    mu  = mean(permMeans,2);
    sig = std(permMeans,0,2) + eps;
    zObs = (obs - mu) ./ sig;
    zPerm = (permMeans - mu) ./ sig;  % T x nPerm

    % cluster-forming threshold (one-sided)
    thr = norminv(1 - alpha_bin, 0, 1);   %
    binSig = zObs > thr;

    % observed cluster masses (sum of z within contiguous supra-threshold bins)
    obsClusters = contiguous_clusters(binSig);
    obsMasses   = cellfun(@(idx) sum(zObs(idx)), obsClusters);

    % permutation max-cluster-mass distribution
    maxMass = zeros(nPerm,1);
    for p = 1:nPerm
        b = zPerm(:,p) > thr;
        cl = contiguous_clusters(b);
        if isempty(cl), maxMass(p) = 0;
        else, maxMass(p) = max(cellfun(@(idx) sum(zPerm(idx,p)), cl));
        end
    end

    % determine which observed clusters are significant (family-wise control)
    sigMask = false(1,T);
    for k = 1:numel(obsClusters)
        p_cluster = (1 + sum(maxMass >= obsMasses(k))) / (1 + nPerm);
        if p_cluster < 0.05
            sigMask(obsClusters{k}) = true;
        end
    end
    clusterStats = struct('obsMasses',obsMasses,'maxMass',maxMass);
end

function clusters = contiguous_clusters(mask)
% Returns cell array of index vectors for contiguous true segments
    mask = mask(:)';
    d = diff([false, mask, false]);
    starts = find(d==1);
    ends   = find(d==-1)-1;
    clusters = arrayfun(@(a,b) a:b, starts, ends, 'UniformOutput', false);
end

%% ================== PLOTTING ==================

function idx = first_run(mask, minRun)
% first index where mask has >= minRun consecutive trues; NaN if none
    clusters = contiguous_clusters(mask);
    idx = NaN;
    for i = 1:numel(clusters)
        if numel(clusters{i}) >= minRun
            idx = clusters{i}(1);
            return;
        end
    end
end

function shade_sig_regions(t, y, mask, col)
% Shade regions where mask==true under the curve
    yl = ylim; y0 = min(yl(1), min(y)); yTop = max(yl(2), max(y));
    cl = contiguous_clusters(mask);
    for i = 1:numel(cl)
        ii = cl{i};
        patch([t(ii(1)) t(ii) t(ii(end))], ...
              [y0 y(ii)' y0], ...
              col, 'FaceAlpha', 0.15, 'EdgeColor','none');
    end
    uistack(findobj(gca,'Type','line'),'top');
end


