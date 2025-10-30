% 
% %% This script will calculate an omnibus statistic (chi-sq) to tell you if 
% %% 1) fraction of significantly modulated cells is same v. different
% %%     between regions
% %%
% %% 2) posthoc tests will tell you which region pairs had different
% %%     fractions of modulated cells 
% 
% %% loads '(date_subj_)facePrefAnova,.mat', which is the output of facePrefAnova.m
% %% (which performed a 2-way anova (bhv, time, interaction) on FRs per cell
% 
% %%%%%%%% GRI 09/11/2025 

clear all; 
subjects = {'barney','thor'};

pathToMatFiles = '/Users/geena/Dropbox/PhD/SUAinfo/';
pThres = 0.05;

%% get barney's data 
barneyRawData1 = load(fullfile(pathToMatFiles, 'Barney_210704/Data4Analysis/210704_barney_facePrefAnova.mat'));
barneyRawData2 = load(fullfile(pathToMatFiles, 'Barney_210706/Data4Analysis/210706_barney_facePrefAnova.mat'));
barneyRawData3 = load(fullfile(pathToMatFiles, 'Barney_210805/Data4Analysis/210805_barney_facePrefAnova.mat'));

barneyRawData.allResults = [barneyRawData1.allResults, barneyRawData2.allResults, barneyRawData3.allResults];
barneyRawData.allPvals = [barneyRawData1.allPvals, barneyRawData2.allPvals, barneyRawData3.allPvals];
barneyRawData.spikeLabels2plot = [barneyRawData1.spikeLabels2plot; barneyRawData2.spikeLabels2plot; barneyRawData3.spikeLabels2plot];

%% get thor's data 
thorRawData1 = load(fullfile(pathToMatFiles, 'Thor_171010/Data4Analysis/171010_Thor_facePrefAnova.mat'));
thorRawData2 = load(fullfile(pathToMatFiles, 'Thor_171027/Data4Analysis/171027_Thor_facePrefAnova.mat'));
thorRawData3 = load(fullfile(pathToMatFiles, 'Thor_171005/Data4Analysis/171005_Thor_facePrefAnova.mat'));
thorRawData4 = load(fullfile(pathToMatFiles, 'Thor_171128/Data4Analysis/171128_Thor_facePrefAnova.mat'));

thorRawData.allResults = [thorRawData1.allResults, thorRawData2.allResults, thorRawData3.allResults, thorRawData4.allResults];
thorRawData.allPvals = [thorRawData1.allPvals, thorRawData2.allPvals, thorRawData3.allPvals, thorRawData4.allPvals];
thorRawData.spikeLabels2plot = [thorRawData1.spikeLabels2plot; thorRawData2.spikeLabels2plot; thorRawData3.spikeLabels2plot; thorRawData4.spikeLabels2plot];

%% get the cortical regions per cell 
cortRegionOut = [];

for ss = 1:length(subjects)
    if strcmpi(subjects{ss},'Barney')
        data = barneyRawData;
    elseif strcmpi(subjects{ss},'Thor')
        data = thorRawData;
    end

    % pool all the p-values 
    for ii = 1:length(data.allResults)
        
        poolPvals(ii,:) = data.allPvals{1,ii};
        
    end
    poolPvalsOut{ss}(:,:) = poolPvals(~isnan(poolPvals(:,1)),:);
    spikeLabels2plot = data.spikeLabels2plot;
    
    % get region mapping
    [regions] = getChannel2CorticalRegionMapping(subjects{ss}, 1);
    
    % code the cortical region for each cell 
    for unit = 1:length(spikeLabels2plot)
        el = spikeLabels2plot(unit,1);
       
        if ismember(el, regions{1,1}.channels)
            cortRegion(unit) = 1; % S1
        elseif ismember(el, regions{1,2}.channels)
            cortRegion(unit) = 2; % M1
        elseif ismember(el, regions{1,3}.channels)
            cortRegion(unit) = 3; % PMv
        elseif ismember(el, regions{1,4}.channels)
            cortRegion(unit) = 4;
        else
            disp('warning!');
        end
       
    end

    cortRegionOut = [cortRegionOut; cortRegion'];

    clear cortRegion;
    clear poolPvals;
end

% threshold significance for two ANOVA factors + their interactions
% 
sigBhv = [poolPvalsOut{1}(:,1) ;poolPvalsOut{2}(:,1)] <= pThres;
sigTime = [poolPvalsOut{1}(:,2) ;poolPvalsOut{2}(:,2)] <= pThres;
sigInteract = [poolPvalsOut{1}(:,3) ;poolPvalsOut{2}(:,3)]  <= pThres;

% run the chi-square tests:
[tblBhv,chi2statBhv,pvalBhv] = crosstab(sigBhv,cortRegionOut);

[tblTime,chi2statTime,pvalTime] = crosstab(sigTime,cortRegionOut);

[tblInt,chi2statInt,pvalInt] = crosstab(sigInteract,cortRegionOut);

% now set-up the post-hoc tests 

S1Time = sigTime(cortRegionOut == 1)';
M1Time = sigTime(cortRegionOut == 2)';
PMvTime = sigTime(cortRegionOut == 3)';
M3Time = sigTime(cortRegionOut == 4)';

S1Int = sigInteract(cortRegionOut == 1)';
M1Int = sigInteract(cortRegionOut == 2)';
PMvInt = sigInteract(cortRegionOut == 3)';
M3Int = sigInteract(cortRegionOut == 4)';

%% time posthocs 
% S1 vs M1 time
[tbl_S1M1time,chi2stat_S1M1time,pval_S1M1time,labels] = crosstab([S1Time M1Time],sort(cortRegionOut(cortRegionOut == 1 | cortRegionOut == 2)));
% S1 v PMv time
[tbl_S1PMvtime,chi2stat_S1PMvtime,pval_S1PMvtime] = crosstab([S1Time PMvTime],sort(cortRegionOut(cortRegionOut == 1 | cortRegionOut == 3)));
% S1 v M3 time
[tbl_S1M3time,chi2stat_S1M3time,pval_S1M3time] = crosstab([S1Time M3Time],sort(cortRegionOut(cortRegionOut == 1 | cortRegionOut == 4)));
% M1 v PMv time
[tbl_M1PMvtime,chi2stat_M1PMvtime,pval_M1PMvtime] = crosstab([M1Time PMvTime],sort(cortRegionOut(cortRegionOut == 2 | cortRegionOut == 3)));
% M1 v M3 time
[tbl_M1M3time,chi2stat_M1M3time,pval_M1M3time] = crosstab([M1Time M3Time],sort(cortRegionOut(cortRegionOut == 2 | cortRegionOut == 4)));
% PMv v M3 time
[tbl_PMvM3time,chi2stat_PMvM3time,pval_PMvM3time] = crosstab([PMvTime M3Time],sort(cortRegionOut(cortRegionOut == 3 | cortRegionOut == 4)));

%% interaction posthocs
% S1 vs M1 interact
[tbl_S1M1int,chi2stat_S1M1int,pval_S1M1int] = crosstab([S1Int M1Int],sort(cortRegionOut(cortRegionOut == 1 | cortRegionOut == 2)));
% S1 v PMv interact
[tbl_S1PMvint,chi2stat_S1PMvint,pval_S1PMvint] = crosstab([S1Int PMvInt],sort(cortRegionOut(cortRegionOut == 1 | cortRegionOut == 3)));
% S1 v M3 interact
[tbl_S1M3int,chi2stat_S1M3int,pval_S1M3int] = crosstab([S1Int M3Int],sort(cortRegionOut(cortRegionOut == 1 | cortRegionOut == 4)));
% M1 v PMv interact
[tbl_M1PMvint,chi2stat_M1PMvint,pval_M1PMvint] = crosstab([M1Int PMvInt],sort(cortRegionOut(cortRegionOut == 2 | cortRegionOut == 3)));
% M1 v M3 interact
[tbl_M1M3int,chi2stat_M1M3int,pval_M1M3int] = crosstab([M1Int M3Int],sort(cortRegionOut(cortRegionOut == 2 | cortRegionOut == 4)));
% PMv v M3 interact
[tbl_PMvM3int,chi2stat_PMvM3int,pval_PMvM3int] = crosstab([PMvInt M3Int],sort(cortRegionOut(cortRegionOut == 3 | cortRegionOut == 4)));

%% FDR 
pairwiseP = {
    pval_S1M1time, pval_S1PMvtime, pval_S1M3time, ...
    pval_M1PMvtime, pval_M1M3time, pval_PMvM3time, ...
    pval_S1M1int,  pval_S1PMvint,  pval_S1M3int,  ...
    pval_M1PMvint, pval_M1M3int,   pval_PMvM3int};

[qBH, sigBH, critBH] = bh_fdr_from_cell(pairwiseP, 0.05);

%% ---------- Build a single tidy table: Stats ----------
regionLabels = ["S1","M1","PMv","M3"];
alpha = 0.05;

% Region index vectors
idxR = arrayfun(@(r) find(cortRegionOut == r), 1:4, 'uni', 0);

% Helpers: per-region N and fraction significant
propPerRegion = @(sigX) cellfun(@(ix) mean(sigX(ix)), idxR);   % fractions 0..1
nPerRegion    = @(sigX) cellfun(@(ix) numel(ix),      idxR);

% ---------- Overall rows ----------
Overall_Beh = makeOverallRow("Behavior", tblBhv,  chi2statBhv,  pvalBhv,  sigBhv,  alpha, regionLabels, nPerRegion, propPerRegion);
Overall_Time= makeOverallRow("Time",     tblTime, chi2statTime, pvalTime, sigTime, alpha, regionLabels, nPerRegion, propPerRegion);
Overall_Int = makeOverallRow("Interact", tblInt,  chi2statInt,  pvalInt,  sigInteract, alpha, regionLabels, nPerRegion, propPerRegion);

% ---------- Pairwise rows (Time & Interact) ----------
pairs     = ["S1 vs M1","S1 vs PMv","S1 vs M3","M1 vs PMv","M1 vs M3","PMv vs M3"];
timeChi2  = [chi2stat_S1M1time, chi2stat_S1PMvtime, chi2stat_S1M3time, chi2stat_M1PMvtime, chi2stat_M1M3time, chi2stat_PMvM3time];
timeP_raw = [pval_S1M1time,     pval_S1PMvtime,     pval_S1M3time,     pval_M1PMvtime,     pval_M1M3time,     pval_PMvM3time];
intChi2   = [chi2stat_S1M1int,  chi2stat_S1PMvint,  chi2stat_S1M3int,  chi2stat_M1PMvint,  chi2stat_M1M3int,  chi2stat_PMvM3int];
intP_raw  = [pval_S1M1int,      pval_S1PMvint,      pval_S1M3int,      pval_M1PMvint,      pval_M1M3int,      pval_PMvM3int];

% Flatten BH outputs
qBH_vec   = qBH(:);      % force column vector
sigBH_vec = sigBH(:);    % force column vector

q_time    = qBH_vec(1:6);   sig_timeBH = sigBH_vec(1:6);
q_int     = qBH_vec(7:12);  sig_intBH  = sigBH_vec(7:12);

Pair_Time = makePairRows("Time",     pairs, timeChi2, timeP_raw, q_time, sig_timeBH, alpha);
Pair_Int  = makePairRows("Interact", pairs, intChi2,  intP_raw,  q_int,  sig_intBH,  alpha);

% ---------- Assemble ----------
Stats = [Overall_Beh; Overall_Time; Overall_Int; Pair_Time; Pair_Int];

% (Optional) save
 save('matfiles/chi2_summary_full.mat','Stats');
 writetable(Stats,'matfiles/chi2_summary_full.csv');

%% ---------- Local functions ----------
function T = makeOverallRow(factorName, chi2tbl, chi2stat, pval, sigVec, alpha, regionLabels, nPerRegion, propPerRegion)
    nVec   = nPerRegion(sigVec);
    pVec   = propPerRegion(sigVec);
    df     = (size(chi2tbl,1)-1)*(size(chi2tbl,2)-1);
    T = table( ...
        "Overall", string(factorName), string(strjoin(regionLabels,",")), ...
        chi2stat, df, pval, NaN, alpha, double(pval < alpha), ...
        nVec(1), pVec(1), nVec(2), pVec(2), nVec(3), pVec(3), nVec(4), pVec(4), ...
        'VariableNames', {'Level','Factor','Comparison', ...
                          'Chi2','df','p_raw','p_adj','alpha','Sig', ...
                          'N_S1','PropSig_S1','N_M1','PropSig_M1', ...
                          'N_PMv','PropSig_PMv','N_M3','PropSig_M3'});
end

function T = makePairRows(factorName, pairNames, chi2v, p_raw, q_adj, sigBH, alpha)
    n = numel(pairNames);
    T = table( ...
        repmat("Pairwise", n,1), repmat(string(factorName), n,1), string(pairNames(:)), ...
        chi2v(:), repmat(1,n,1), p_raw(:), q_adj(:), repmat(alpha,n,1), double(sigBH(:)), ...
        NaN(n,1), NaN(n,1), NaN(n,1), NaN(n,1), NaN(n,1), NaN(n,1), NaN(n,1), NaN(n,1), ...
        'VariableNames', {'Level','Factor','Comparison', ...
                          'Chi2','df','p_raw','p_adj','alpha','Sig', ...
                          'N_S1','PropSig_S1','N_M1','PropSig_M1', ...
                          'N_PMv','PropSig_PMv','N_M3','PropSig_M3'});
end



function [qvals, sigMask, crit_p, pvec, idx_orig] = bh_fdr_from_cell(p_in, alpha)
% bh_fdr_from_cell  Benjamini–Hochberg FDR from cell array or numeric vector.
% USAGE:
%   [qvals, sigMask, crit_p] = bh_fdr_from_cell(p_in, alpha)
% INPUTS:
%   p_in  : cell array of p-values (any shape) OR numeric vector/matrix
%   alpha : desired FDR level (e.g., 0.05). Default = 0.05
% OUTPUTS:
%   qvals    : FDR-adjusted p-values (same size as p_in)
%   sigMask  : logical mask of discoveries under BH at 'alpha' (same size)
%   crit_p   : BH critical p-value (largest p that passes)
%   pvec     : sorted p-values used internally (column vector, NaNs removed)
%   idx_orig : indices mapping sorted p-values back to linear indices of p_in

    if nargin < 2 || isempty(alpha), alpha = 0.05; end

    if iscell(p_in)
        p_all = cellfun(@(x) double(x), p_in, 'UniformOutput', true);
    else
        p_all = double(p_in);
    end
    sz = size(p_all);
    p_lin = p_all(:);

    % Track which entries are valid numbers
    valid = isfinite(p_lin);
    p_vec = p_lin(valid);

    % Edge cases
    qvals = nan(sz);
    sigMask = false(sz);
    crit_p = NaN;
    if isempty(p_vec)
        return
    end

    % Sort ascending
    [p_sorted, sort_idx] = sort(p_vec, 'ascend');
    m = numel(p_sorted);
    ranks = (1:m)';               % 1..m

    % Compute BH q-values (adjusted p)
    q_sorted = (m ./ ranks) .* p_sorted;

    q_sorted = min(1, cummin(q_sorted(end:-1:1)));
    q_sorted = q_sorted(end:-1:1);

    q_vec = nan(size(p_vec));
    q_vec(sort_idx) = q_sorted;

    % Significance under BH at alpha
    thresh = (ranks / m) * alpha;
    pass = p_sorted <= thresh;
    if any(pass)
        k = find(pass, 1, 'last');     % largest k meeting the criterion
        crit_p = p_sorted(k);
        sig_sorted = false(m,1);
        sig_sorted(1:k) = true;
        % map back to original order among valid entries
        sig_vec = false(size(p_vec));
        sig_vec(sort_idx) = sig_sorted;
    else
        sig_vec = false(size(p_vec));
        crit_p = NaN;  % nothing passes
    end

    % Write back into full-size outputs
    q_out = nan(numel(p_lin),1);  q_out(valid) = q_vec;
    s_out = false(numel(p_lin),1); s_out(valid) = sig_vec;

    qvals(:)   = reshape(q_out, sz);
    sigMask(:) = reshape(s_out, sz);

    % Optional: also return the sorted valid p-values and their original indices
    % (linear indices into p_in)
    idx_valid_lin = find(valid);
    idx_orig = idx_valid_lin(sort_idx);
end
