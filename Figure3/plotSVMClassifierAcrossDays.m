%%%%%%%%
%%%%%%%% PLOTTING Figure 3B, S3A, S3B
%%%%%%%%  Decoding gesture type from neural population activity
%%%%%%%%
%%%%%%%% Categorical facial gesture decoded from region-specific population spiking activity on a trial-by-trial basis
%%%%%%%% Mean decoding accuracies ± 2 SEMs , across pseudopopulations (N=50 iterations, 50 cells per iteration, drawn from entire database)
%%%%%%%% from -1000 ms to +1000 ms (400 ms bins in 50 ms steps)
%%%%%%%%
%%%%%%%% Grey curve indicates significance threshold based on 98th percentile of one hundred permuted decoding accuracies 
%%%%%%%%
%%%%%%%%  Statistical significance determined as (Maris & Oostenveld, 2007), by comparing per region decoding at each timepoint to the 98th percentile of permuted decoding accuracies 
%%%%%%%   Δ(t)=Accreal (t)−q98null (t) per pseudopopulation)
%%%%%%%%  one-sample (right-tailed) t-test within a cluster-based permutation framework 
%%%%%%%%  At each time point, across-iteration t-values computed for real decoding accuracy relative to its permutation-derived baseline (random sign-flipping of decoding time course), contiguous supra-threshold bins grouped into clusters. Cluster masses were compared against the distribution of maximum cluster masses obtained by permutations, providing family-wise error control over time. 

%%%%%%%%  Will also plot decoding accuracies, compared pairwise by region
%%%%%%%%  Black bars indicate significant difference in accuracy; two-sample Welch s t-test, cluster-based permutation correction for multiple comparisons
%%%%%%%%
%%%%%%%   Input data contains data.results (true decoding results of N= 50
%%%%%%%   pseudopopulations) and data.perm (null results from 100
%%%%%%%   permutations per population) 
%%%%%%%%% GRI 09/11/2025 

clear all; close all; 
set(0,'defaultAxesFontSize',20)

%% intial set-up 
nIterations = 50;
windowLength = 0.4; windowShiftSize = 0.05;
colorArray = [0.4940 0.1840 0.5560;0.6350 0.0780 0.1840;0.8500 0.3250 0.0980;0.9290 0.6940 0.1250];

%% load data 
data = load('matfiles/Fig3B.mat');
results = data.results;

%% -------- Cluster-based permutation parameters --------
alpha_cf = 0.01;     % cluster-forming p-threshold (per-time)
nPerm    = 10000;     % number of permutations for max cluster mass null
useDirectionalBetweenShading = false;  % true: shade by winner color; false: blended color
errorType = 1; % 2SEMs across iter 

%% define time axis 
firstWindowCenter = results(1).desiredStart + results(1).windowLength/2;
lastWindowCenter = results(1).desiredEnd - results(1).windowLength/2;
windowsCenter = firstWindowCenter:results(1).windowShiftSize:lastWindowCenter;

arrayList = fieldnames(results(1).validationScoreStruct);
validationScores = [results(:).validationScoreStruct];
SEMScores = [results(:).SEM];
permScores = [data.perm(:).permValidationScoreStruct];

%% get averages across N pseudopopulations and define errorBars
for aa = 1:length(arrayList)

    allIterResults = [];
    allIterSEMs = [];

    for iter = 1:nIterations

        thisIterResults = validationScores(iter).(arrayList{aa});
        allIterResults = [thisIterResults; allIterResults];

        thisIterSEM = SEMScores(iter).(arrayList{aa});
        allIterSEMs = [thisIterSEM; allIterSEMs];
        
    end

    meanIterValidationScores(aa,:) = mean(allIterResults,1);
    semErrBars(aa,:) = ((std(allIterResults,1) ./ sqrt(nIterations))) ;
    semErrBars_2(aa,:) = (mean(allIterSEMs,1));
    stdErrBars(aa,:) = (std(allIterResults,[],1));

end

if errorType == 1
    errorBars = semErrBars;
elseif errorType == 2
    errorBars = semErrBars_2;
elseif errorType == 3
    errorBars = stdErrBars;

end

%% permutations 98th% as cutoff for significance of decoder
prcThres = 98;
P = zeros(length(arrayList),length(windowsCenter),nIterations);

for aa = 1:length(arrayList)

    for iter = 1:nIterations

        for tt = 1:length(windowsCenter)

            P(aa,tt,iter) = prctile(permScores(iter).(arrayList{aa})(tt,:),prcThres);

        end

    end

end


%% plot decoding accuracy by region over time 
fig = figure; fig.Position =[476 92 1037 774];

for aa = 1:length(arrayList)
  
    chance = mean(P(aa,:,:),3)*100;
    chanceSEM = (std(P(aa,:,:),[],3)./ sqrt(nIterations) ) *100;
    
    errorbar(windowsCenter, meanIterValidationScores(aa,:)*100, errorBars(aa,:)*100,...
        'Color', colorArray(aa,:), 'LineWidth',8); hold on
    errorbar(windowsCenter,chance, chanceSEM,'Color', [0.5 0.5 0.5],'LineWidth',3); hold on
  
end
hold off
xlim([- 0.8 0.8]);
axes = gca;
axes.FontWeight = 'bold';
axes.FontSize = 34;
axes.XTick = windowsCenter(1:4:end);
xlabel('Time (s)')
ylabel('Mean Accuracy (%)')
ylim([50 103])
title('Decoding Performance By Region','FontSize',45);


%% Collect iteration matrices per region (real) and per-iteration null means
nTime = numel(windowsCenter);
nReg  = numel(arrayList);

realAcc = cell(nReg,1);        % {aa}: nIter x nTime (real)

%% -------- Build real matrices 

for r = 1:nReg
    A = zeros(nIterations, nTime);
    B = zeros(nIterations, nTime);
    for iter = 1:nIterations
        
        A(iter,:) = validationScores(iter).(arrayList{r})(:).';
        
        % Null: nTime x nPerms -> mean over perms (per-iteration null)
        P = permScores(iter).(arrayList{r});   % nTime x nPerms
        B(iter,:) = mean(P, 2, 'omitnan').';
    end
    realAcc{r}  = A;
end

%% Per-iteration 98th percentile of the permuted (null) accuracies
nullQ98  = cell(nReg,1);          % {r}: nIter x nTime
permCell = cell(nReg,1);          % {r}: 1 x nIter cell of [nTime x nPerms] 

for r = 1:nReg
    Q = zeros(nIterations, nTime);
    C = cell(1, nIterations);
    for iter = 1:nIterations
        P = permScores(iter).(arrayList{r});   % [nTime x nPerms]
        Q(iter,:) = prctile(P, 98, 2).';       % upper 98th percentile per time (row)
        C{iter}   = P;                         % keep all perms for the cluster null
    end
    nullQ98{r}  = Q;
    permCell{r} = C;
end

%% -------- Within-region, right-tailed: real > q98(null), cluster-based permutation
within = cell(nReg,1);
for r = 1:nReg
    within{r} = cluster_perm_within_percentile(realAcc{r}, permCell{r}, 0.98, ...
                    'alpha_cf', alpha_cf, 'nPerm', nPerm);
    within{r}.label = arrayList{r};
end


%% -------- Between-region tests (two-group; two-sided) --------
pairs = nchoosek(1:nReg, 2);
between = cell(size(pairs,1),1);
for p = 1:size(pairs,1)
    i = pairs(p,1); j = pairs(p,2);
    between{p} = cluster_perm_2group(realAcc{i}, realAcc{j}, 'tail','both', 'alpha_cf',alpha_cf, 'nPerm',nPerm);
    between{p}.labelA = arrayList{i};
    between{p}.labelB = arrayList{j};
end

clear allIterResults allIterSEMs meanIterValidationScores permScores perm validationScores; 

%% -------- Plot: regions vs null (one figure, subplots) --------
plot_regions_vs_null_grid(windowsCenter, realAcc, nullQ98, within, arrayList, colorArray);

%% -------- Plot: between-region contrasts (one figure, subplots) --------
plot_between_regions_grid(windowsCenter, realAcc, between, pairs, arrayList, colorArray, useDirectionalBetweenShading);

%% ----------- Diagnostics (optional) -----------
% Quick textual summary of between-region clusters 

for p = 1:numel(between)
    i = pairs(p,1); j = pairs(p,2);
    S = between{p};
    fprintf('\n=== %s vs %s ===\n', arrayList{i}, arrayList{j});
    fprintf('Any cluster-corrected bins? %s\n', string(any(S.sig_mask)));
    if ~isempty(S.clusters)
        for c = 1:numel(S.clusters)
            seg = S.clusters(c).start:S.clusters(c).end;
            fprintf('  Cluster %d: %d bins, [%g..%g], mass=%.3f, p=%.4g\n', ...
                c, numel(seg), windowsCenter(seg(1)), windowsCenter(seg(end)), S.clusters(c).mass, S.clusters(c).p);
        end
    else
        fprintf('  (No significant clusters)\n');
    end
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  PERMUTATION BUILDERS                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function out = cluster_perm_within_percentile(realMat, permCell, q, varargin)
% Cluster permutation vs per-iteration null q-quantile baseline.
% realMat : [nIter x T] real decoding accuracies
% permCell: 1 x nIter cell, each [T x nPerms] permuted accuracies for that iteration
% q       : quantile for cut off 

p = inputParser;
addParameter(p,'alpha_cf',0.05,@(x)isnumeric(x)&&isscalar(x)&&x>0&&x<1);
addParameter(p,'nPerm',5000,@(x)isnumeric(x)&&isscalar(x)&&x>=100);
parse(p,varargin{:});
alpha_cf = p.Results.alpha_cf;
nPerm    = p.Results.nPerm;

[nIter, T] = size(realMat);

% --- Build observed baseline (per-iter q-quantile) and difference D ---
Q = zeros(nIter, T);
for i = 1:nIter
    P = permCell{i};                 % [T x nPerms]
    Q(i,:) = prctile(P, q*100, 2).';
end
Dobs = realMat - Q;                  % positive means real > q-quantile
[m, s, n] = msd(Dobs);
tObs = m ./ (s ./ sqrt(n));          % one-sample t across iterations
df   = max(n - 1, 1);

% Cluster-forming mask from per-time right-tailed p
pRight = 1 - tcdf(tObs, df);
supra  = (pRight < alpha_cf);
obsClusters = label_clusters(supra);
for c = 1:numel(obsClusters)
    seg = obsClusters(c).start:obsClusters(c).end;
    obsClusters(c).mass = sum(tObs(seg));     % positive mass (right-tailed)
end

% --- Permutation null: choose one permuted trajectory per iteration ---
maxMass0 = zeros(nPerm,1);
for pidx = 1:nPerm
    Dp = zeros(nIter, T);
    for i = 1:nIter
        P = permCell{i};                        % [T x nPerms]
        k = randi(size(P,2));                  
        permTraj = P(:,k).';
        % Use the *same* Q baseline (built from all perms) for speed; 
        Dp(i,:) = permTraj - Q(i,:);
    end
    [m0, s0, n0] = msd(Dp);
    t0  = m0 ./ (s0 ./ sqrt(n0));
    df0 = max(n0 - 1, 1);
    pR0 = 1 - tcdf(t0, df0);
    supra0 = (pR0 < alpha_cf);
    maxMass0(pidx) = max_cluster_mass(t0, supra0, @(u) sum(u));  % right-tailed mass
end

% --- Cluster p-values & sig mask ---
sig_mask = false(1,T);
for c = 1:numel(obsClusters)
    cmass = obsClusters(c).mass;
    pCl = (sum(maxMass0 >= cmass) + 1) / (nPerm + 1);
    obsClusters(c).p = pCl;
    if pCl < 0.05
        sig_mask(obsClusters(c).start:obsClusters(c).end) = true;
    end
end

out.t        = tObs;
out.sig_mask = sig_mask;
out.clusters = obsClusters;
out.maxMass0 = maxMass0;
out.info     = struct('q',q,'alpha_cf',alpha_cf,'nPerm',nPerm,'nIter',nIter,'T',T);
end

function out = cluster_perm_2group(A, B, varargin)

% Two-group cluster permutation across time (comparing region-region
% decoding) 
p = inputParser;
addParameter(p,'tail','both',@(s)ischar(s)||isstring(s));
addParameter(p,'alpha_cf',0.05,@(x)isnumeric(x)&&isscalar(x)&&x>0&&x<1);
addParameter(p,'nPerm',5000,@(x)isnumeric(x)&&isscalar(x)&&x>=100);
parse(p,varargin{:});
tail     = lower(string(p.Results.tail));
alpha_cf = p.Results.alpha_cf;
nPerm    = p.Results.nPerm;

% bookeeping
nA_t = sum(isfinite(A),1);
nB_t = sum(isfinite(B),1);
keep = (nA_t >= 2) & (nB_t >= 2);
A = A(:,keep);  B = B(:,keep);
[nA,T] = size(A);  nB = size(B,1);

% Observed t, df, and two-sided p (per time), the stat 
[tObs, p2Obs, dfObs] = welch_t(A,B);

% Cluster-forming mask from per-time p-values 
switch tail
    case "both"
        supra = (p2Obs < alpha_cf);
        massfun = @(t) sum(abs(t));
        tForMass = abs(tObs);
    case "right"
        pRight = 1 - tcdf(tObs, dfObs);
        supra  = (pRight < alpha_cf);
        massfun = @(t) sum(t);
        tForMass = tObs;
    case "left"
        pLeft = tcdf(tObs, dfObs);
        supra = (pLeft < alpha_cf);
        massfun = @(t) sum(t);
        tForMass = -tObs;   % flip so mass is positive when B>A
    otherwise
        error('tail must be both/right/left');
end

obsClusters = label_clusters(supra);
for c = 1:numel(obsClusters)
    seg = obsClusters(c).start:obsClusters(c).end;
    obsClusters(c).mass = massfun(tForMass(seg));
end

% Permutation: shuffle whole decoding trajecotries; recompute our
% statistics each time 
X = [A; B]; nTot = size(X,1);
maxMass0 = zeros(nPerm,1);
for pidx = 1:nPerm
    idx = randperm(nTot);
    Ai = X(idx(1:nA),:);
    Bi = X(idx(nA+1:end),:);
    [t0, p20, df0] = welch_t(Ai,Bi);
    switch tail
        case "both"
            supra0 = (p20 < alpha_cf);
            mass0  = max_cluster_mass(abs(t0), supra0, @(u) sum(u));
        case "right"
            pR0    = 1 - tcdf(t0, df0);
            supra0 = (pR0 < alpha_cf);
            mass0  = max_cluster_mass(t0, supra0, @(u) sum(u));
        case "left"
            pL0    = tcdf(t0, df0);
            supra0 = (pL0 < alpha_cf);
            mass0  = max_cluster_mass(-t0, supra0, @(u) sum(u));
    end
    maxMass0(pidx) = mass0;
end

sig_mask = false(1,T);
for c = 1:numel(obsClusters)
    cmass = obsClusters(c).mass;
    pcl = (sum(maxMass0 >= cmass) + 1) / (nPerm + 1);
    obsClusters(c).p = pcl;
    if pcl < 0.05
        sig_mask(obsClusters(c).start:obsClusters(c).end) = true;
    end
end

% Pack back into output 
fullT = false(1, numel(keep)); fullT(keep) = true;
out.t        = nan(1, numel(keep)); out.t(keep) = tObs;
out.p_unc    = nan(1, numel(keep)); out.p_unc(keep) = p2Obs;
out.sig_mask = false(1, numel(keep)); out.sig_mask(keep) = sig_mask;
out.clusters = obsClusters;
out.maxMass0 = maxMass0;
out.info     = struct('tail',tail,'alpha_cf',alpha_cf,'nPerm',nPerm,'nA',nA,'nB',nB,'T',T,'keptMask',fullT);

end

function [t, p2, df] = welch_t(A,B)
% Vectorized Welch t across time with per-time df and two-sided p
mA = mean(A,1,'omitnan'); vA = var(A,0,1,'omitnan'); nA = sum(isfinite(A),1);
mB = mean(B,1,'omitnan'); vB = var(B,0,1,'omitnan'); nB = sum(isfinite(B),1);
se = sqrt(vA./nA + vB./nB);
t  = (mA - mB) ./ se;

% Welch–Satterthwaite df (per time); guard small denominators
num = (vA./nA + vB./nB).^2;
den = ((vA./nA).^2) ./ max(nA-1,1) + ((vB./nB).^2) ./ max(nB-1,1);
df  = num ./ max(den, eps);
p2  = 2 * (1 - tcdf(abs(t), df));
end

function [m,s,n] = msd(X)
m = mean(X,1,'omitnan');
s = std(X,0,1,'omitnan');
n = sum(isfinite(X),1);
s(s==0) = inf; % t=0 where variance is 0
end

function clusters = label_clusters(mask)
idx = find(mask(:)');
clusters = struct('start',{},'end',{},'mass',{});
if isempty(idx), return; end
splits = [1, find(diff(idx) > 1) + 1, numel(idx)+1];
for s = 1:numel(splits)-1
    seg = idx(splits(s):splits(s+1)-1);
    clusters(end+1).start = seg(1); 
    clusters(end).end     = seg(end);
end
end

function m = max_cluster_mass(t, mask, massfun)
cls = label_clusters(mask);
if isempty(cls), m = 0; return; end
ms = zeros(numel(cls),1);
for c = 1:numel(cls)
    seg = cls(c).start:cls(c).end;
    ms(c) = massfun(t(seg));
end
m = max(ms);
end


%% -------------------- Plotting helpers --------------------

function plot_regions_vs_null_grid(tCenters, realAcc, nullMean, within, arrayList, colorArray)
nReg = numel(arrayList);
[rows, cols] = bestSubplotGrid(nReg);
fig = figure('Color','w','Name','Regions vs Null (cluster-corrected)');
fig.Position =[476 92 1037 774];

for r = 1:nReg
    ax = subplot(rows, cols, r, 'Parent', fig); hold(ax,'on');
    c  = colorArray(min(r,size(colorArray,1)),:);
    cDark  = darkenColor(c, 0.45);   % shading
    cLight = lightenColor(c, 0.55);  % null line

    mR = mean(realAcc{r},1,'omitnan'); sR = std(realAcc{r},0,1,'omitnan') ./ sqrt(size(realAcc{r},1));
    mN = mean(nullMean{r},1,'omitnan'); sN = std(nullMean{r},0,1,'omitnan') ./ sqrt(size(nullMean{r},1));

    hR = errorbar(ax, tCenters, 100*mR, 100*sR, '-', 'LineWidth',1.8, 'Color', c,      'CapSize',0);
    hN = errorbar(ax, tCenters, 100*mN, 100*sN, '-', 'LineWidth',2, 'Color', 'k', 'CapSize',0, 'LineStyle','-.');

    ylim([50 100]);
    % Shade significant bins (after lines; uses bin edges, handles 1-bin clusters)
    shade_sig_ax(ax, tCenters, within{r}.sig_mask, cDark, 0.30);

    if r==1
        legend(ax, [hR hN], {'real','null mean'}, 'Location','best');
    end
    xlabel(ax,'Time'); ylabel(ax,'Decoding accuracy');
    title(ax, sprintf('%s: real > null (cluster-corr)', arrayList{r}));
    grid(ax,'on'); box(ax,'off'); 
end
end

function plot_between_regions_grid(tCenters, realAcc, between, pairs, arrayList, colorArray, ~)
% One figure, subplots = between-region contrasts.
% Significance is shown as a horizontal black bar along the bottom (x-axis).
% All subplots share the same y-limits.

nPairs = size(pairs,1);
[rows, cols] = bestSubplotGrid(nPairs);
fig = figure('Color','w','Name','Between-Region Contrasts (cluster-corrected)');
fig = figure; fig.Position =[476 92 1037 774];

% ---------- compute global y-lims across all pairs ----------
yMin =  inf; 
yMax = -inf;
for p = 1:nPairs
    i = pairs(p,1); j = pairs(p,2);
    mA = mean(realAcc{i},1,'omitnan'); sA = std(realAcc{i},0,1,'omitnan')./sqrt(size(realAcc{i},1));
    mB = mean(realAcc{j},1,'omitnan'); sB = std(realAcc{j},0,1,'omitnan')./sqrt(size(realAcc{j},1));
    localMin = min([mA - sA, mB - sB], [], 'all', 'omitnan');
    localMax = max([mA + sA, mB + sB], [], 'all', 'omitnan');
    yMin = min(yMin, localMin);
    yMax = max(yMax, localMax);
end
if ~isfinite(yMin) || ~isfinite(yMax) || yMax <= yMin
    yMin = 0; yMax = 1;   % fallback
end
yRange  = yMax - yMin;
pad     = 0.08 * yRange;       % 8% padding
yMinAll = yMin - pad;
yMaxAll = yMax + pad;
barH    = 0.03 * (yMaxAll - yMinAll);  % bar height = 3% of axis range

% ---------- plot each pair ----------
for p = 1:nPairs
    i = pairs(p,1); j = pairs(p,2);
    ax = subplot(rows, cols, p, 'Parent', fig); hold(ax,'on');

    cA = colorArray(min(i,size(colorArray,1)),:);
    cB = colorArray(min(j,size(colorArray,1)),:);

    mA = mean(realAcc{i},1,'omitnan'); sA = std(realAcc{i},0,1,'omitnan') ./ sqrt(size(realAcc{i},1));
    mB = mean(realAcc{j},1,'omitnan'); sB = std(realAcc{j},0,1,'omitnan') ./ sqrt(size(realAcc{j},1));

    hA = errorbar(ax, tCenters, mA, sA, '-', 'LineWidth',1.8, 'Color', cA, 'CapSize',0);
    hB = errorbar(ax, tCenters, mB, sB, '-', 'LineWidth',1.8, 'Color', cB, 'CapSize',0, 'LineStyle','--');

    % set unified y-lims before drawing significance bar
    ylim(ax, [yMinAll, yMaxAll]);

    % draw horizontal significance bar along bottom (solid black)
    bar_sig_ax(ax, tCenters, between{p}.sig_mask, [0 0 0], yMinAll, barH);

    %if p == 1
        legend(ax, [hA hB], {arrayList{i}, arrayList{j}}, 'Location','best');
   % end
    xlabel(ax,'Time'); ylabel(ax,'Decoding accuracy');
    title(ax, sprintf('%s vs %s (cluster-corr)', arrayList{i}, arrayList{j}));
    grid(ax,'on'); box(ax,'off');
end
end

function bar_sig_ax(ax, tCenters, mask, colorRGB, y0, h)
% Draw solid horizontal bars (rectangles) along bottom for each significant segment.
    if isempty(mask) || ~any(mask), return; end
    edges = bin_edges_from_centers(tCenters);
    segs  = logical_mask_to_segments(mask);
    for s = 1:size(segs,1)
        i0 = segs(s,1);
        i1 = segs(s,2);
        x0 = edges(i0);
        x1 = edges(i1+1);      % close at right edge
        patch('Parent',ax, ...
              'XData',[x0 x1 x1 x0], 'YData',[y0 y0 y0+h y0+h], ...
              'FaceColor', colorRGB, 'EdgeColor','none', 'FaceAlpha', 1.0);
    end
    % keep bars behind data lines but visible
    uistack(findobj(ax,'Type','patch'),'bottom');
end

function shade_sig_ax(ax, tCenters, mask, colorRGB, alphaVal)

% Fill full-height vertical bands for significant clusters
    if isempty(mask) || ~any(mask), return; end
    edges = bin_edges_from_centers(tCenters);
    segs  = logical_mask_to_segments(mask);

    % use current y-limits so it spans the whole plot
    yl = ylim(ax);
    y0 = yl(1); y1 = yl(2);

    for s = 1:size(segs,1)
        i0 = segs(s,1); i1 = segs(s,2);
        x0 = edges(i0); x1 = edges(i1+1);
        patch('Parent',ax, ...
              'XData',[x0 x1 x1 x0], 'YData',[y0 y0 y1 y1], ...
              'FaceColor', colorRGB, 'FaceAlpha', alphaVal, 'EdgeColor','none');
    end
    uistack(findobj(ax,'Type','patch'),'bottom');
end


function edges = bin_edges_from_centers(c)
% Convert center times to edges (works for non-uniform spacing)
c = c(:)'; 
if numel(c) == 1
    d = 1; edges = [c(1)-d/2, c(1)+d/2]; return;
end
mid   = (c(1:end-1) + c(2:end))/2;
left  = c(1)  - (mid(1) - c(1));
right = c(end) + (c(end) - mid(end));
edges = [left, mid, right];
end

function segs = logical_mask_to_segments(mask)
% Return Nx2 [startIdx endIdx] for contiguous true segments (1-based)
mask = mask(:)'; idx = find(mask);
segs = zeros(0,2);
if isempty(idx), return; end
brks = [1, find(diff(idx) > 1)+1, numel(idx)+1];
for b = 1:numel(brks)-1
    seg = idx(brks(b):brks(b+1)-1);
    segs(end+1,:) = [seg(1) seg(end)]; 
end
end

function [rows, cols] = bestSubplotGrid(N)
rows = floor(sqrt(N));
cols = ceil(N / max(rows,1));
if rows*cols < N, rows = rows + 1; end
end

function c2 = darkenColor(c, amt)  % blend toward black
c2 = max(0, (1-amt).*c);
end

function c2 = lightenColor(c, amt) % blend toward white
c2 = min(1, c + (1-c).*amt);
end

