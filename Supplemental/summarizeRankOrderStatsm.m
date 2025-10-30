
%%%%%%%% PLOTTING Figure S6A, S6B 
%%%%%%%%  Kinematic decoding and non-preserved neural correlations
%%%%%%%%
%%%%%%%%  Plots all Spearman's rho values (large for preserved rank-order, small if not) 
%%%%%%% (1 value per day, between pairwise correlations:
%%%%%%        threat-lipsmack, threat-chew, lipsmack-chew); 
%%%%%%% GRI 09/11/2025 
%%%%%%%
%%%%%%% Rank-order similarity of gesture matrix pairs across all neurons, one value per recording and per gesture-matrix comparison (threat–lipsmack, threat–chew, lipsmack–chew). 
%%%%%%%  Per day we computed Spearman%s ρ. 
%%%%%%% Significance for each day's ρ was obtained with a Mantel permutation test on the matrix labels, yielding a permutation p-value. 
%%%%%%% Jittered points show per-day ρ. Circles indicate median ρ across days, error bars are bootstrapped 95% CIs. 
%%%%%%% Asterisks above each pair indicate combined significance across days using Fisher%s method on the per-day Mantel p-values. 
%%%%%%% *** = p < 0.001, ** = p < 0.01, at α=0.05. 

set(0,'defaultAxesFontSize', 24); % bc im blind
set(0,'defaultAxesFontWeight', 'bold'); 
metric = 'Kendall'; % or 'Spearman'; will print these results to console 

combined = summarize_rankorder_results_dirs({''}, 'DayLabel','dir', ...
    'FilePattern','*rankSim*.mat', ... 
    'Recurse',false, ...
    'SaveCSV',true, ...
    'CSVPath','matfiles/combined_rankorder_summary.csv');     

% plot spearman by region summary 
%%
plot_rankorder_by_region(combined, ...
    'Metric', metric, ...
    'RegionSrc','File', ...
    'RegionPatterns', {'All','All'; 'S1','S1'; 'M1','M1'; 'PMv','PMv'; 'M3','M3'}, ...
    'NBoot',5000, 'ShowPoints',true);


% Spearman summary
%plot_rankorder_summary(combined, 'Metric','Spearman', 'Alpha',0.05, 'NBoot',5000, 'ShowPoints',true);

% Kendall summary (optional)
%plot_rankorder_summary(combined, 'Metric','Kendall', 'Alpha',0.05, 'NBoot',5000, 'ShowPoints',true);


function [combined] = summarize_rankorder_results_dirs(dirList, varargin)
% summarize_rankorder_results_dirs
% Walk a list of directories, load rank-order result tables, and print/aggregate.
%
% INPUT
%   dirList : string array or cellstr of directories (absolute or relative)
%
% Name-Value options (all optional):
%   'FilePattern' : '*.mat'             % which files to read in each dir
%   'Recurse'     : false               % search subfolders too (requires '**' support)
%   'SaveCSV'     : true                % write combined CSV
%   'CSVPath'     : 'combined_rankorder_summary.csv'
%   'DayLabel'    : 'dir'               % 'dir' (folder name) or 'file' (file base name)
%
% OUTPUT
%   combined : table with Day + all summary columns stacked across dirs/files
%
% Each result file must contain either:
%   - variable 'out' with field 'summary' (table), OR
%   - variable 'summary' directly (table)

% ------------- options -------------
p = inputParser;
addParameter(p,'FilePattern','*.mat',@(s)ischar(s)||isstring(s));
addParameter(p,'Recurse',false,@islogical);
addParameter(p,'SaveCSV',true,@islogical);
addParameter(p,'CSVPath','combined_rankorder_summary.csv',@(s)ischar(s)||isstring(s));
addParameter(p,'DayLabel','dir',@(s)any(strcmpi(string(s),["dir","file"])));
parse(p,varargin{:});
filePattern = string(p.Results.FilePattern);
recurse     = p.Results.Recurse;
saveCSV     = p.Results.SaveCSV;
csvPath     = string(p.Results.CSVPath);
daySource   = lower(string(p.Results.DayLabel));

% --- fast path: if CSV already exists, just load and return ---
if exist(csvPath, 'file') == 2
    combined = readtable(csvPath);
    fprintf('Existing CSV found. Loaded: %s\n', csvPath);
    return;
end

% required columns
reqCols = {'Pair','Spearman_r','Spearman_p','Spearman_CI_low','Spearman_CI_high', ...
           'Kendall_tau','Kendall_p','Kendall_CI_low','Kendall_CI_high'};

% normalize dirList
if ischar(dirList) || isstring(dirList), dirList = cellstr(dirList); end
dirList = string(dirList(:));

allChunks = {};
fprintf('Scanning %d directories...\n\n', numel(dirList));

for d = 1:numel(dirList)
    % Build the expected folder path
    dpath = fullfile('/Users/geena/Dropbox/PhD/SUAInfo', dirList(d), 'CorrMatricesV2');

    if ~isfolder(dpath)
        warning('Skipping: %s (not a folder)', dpath);
        continue;
    end

    % discover files
    if recurse
        S = dir(fullfile(dpath, '**', filePattern));
    else
        S = dir(fullfile(dpath, filePattern));
    end
    if isempty(S)
        warning('No files in %s matching %s', dpath, filePattern);
        continue;
    end

    % pretty day label (folder name)
    [~, dayLabelFromDir] = fileparts(char(dpath));
    fprintf('=== %s ===\n', dayLabelFromDir);

    for k = 1:numel(S)
        fn = fullfile(S(k).folder, S(k).name);
        data = load(fn);

        % pull summary table
        if isfield(data, 'out') && isfield(data.out, 'summary') && istable(data.out.summary)
            T = data.out.summary;
        elseif isfield(data, 'summary') && istable(data.summary)
            T = data.summary;
        else
            warning('  Skipping %s (no out.summary / summary table).', fn);
            continue;
        end

        % ensure required columns exist
        missing = setdiff(reqCols, T.Properties.VariableNames);
        if ~isempty(missing)
            warning('  Skipping %s (missing columns: %s)', fn, strjoin(missing,', '));
            continue;
        end

        % pick day label
        dirLabel = string(dirList(d));  % e.g., 'Barney_210704'
        if daySource == "file"
            [~, baseName, ~] = fileparts(S(k).name);
            dayLabel = string(baseName);            % label by file basename
        else
            dayLabel = dirLabel;                    % label by dirList entry
        end

        % ---- print human-readable summary for this file/day ----
        fprintf('   Date: %s\n', dayLabel);
        fprintf('  File: %s\n', S(k).name);
        for i = 1:height(T)
            pr  = T.Spearman_r(i);
            psp = T.Spearman_p(i);
            cil = T.Spearman_CI_low(i);
            cih = T.Spearman_CI_high(i);

            kt  = T.Kendall_tau(i);
            pkp = T.Kendall_p(i);
            kil = T.Kendall_CI_low(i);
            kih = T.Kendall_CI_high(i);

            fprintf('    %-7s  ρ=% .3f (95%% CI [% .3f, % .3f]);  τ=% .3f (95%% CI [% .3f, % .3f])\n', ...
                T.Pair{i}, pr, cil, cih, kt, kil, kih);
            fprintf('             Mantel p(ρ)=%s   Mantel p(τ)=%s\n', pformat(psp), pformat(pkp));
        end
        fprintf('\n');

        % augment and stash
        Day   = repmat(dayLabel, height(T), 1);
        File  = repmat(string(S(k).name), height(T), 1);
        Chunk = addvars(T, Day, File, 'Before', 1);
        allChunks{end+1} = Chunk; %#ok<AGROW>
    end
end

if isempty(allChunks)
    warning('No valid result tables found. Returning empty table.');
    combined = table();
    return;
end

combined = vertcat(allChunks{:});
combined = movevars(combined, 'Day', 'Before', 1);

% save CSV
if saveCSV
    writetable(combined, csvPath);
    fprintf('Combined table saved to: %s\n', csvPath);
end
end


function plot_rankorder_by_region(combined, varargin)
% plot_rankorder_by_region  Make one summary plot per cortical region.
%
% Usage:
%   plot_rankorder_by_region(combined, ...
%       'Metric','Spearman', ...          % or 'Kendall'
%       'RegionSrc','File', ...           % which column to parse: 'File' or 'Day' or a custom col
%       'RegionPatterns', {'All','All'; 'S1','S1'; 'M1','M1'; 'PMv','PMv'; 'M3','M3'}, ...
%       'Alpha',0.05, 'NBoot',5000, 'ShowPoints',true)
%
% Notes:
%   - Collapses across days within each region.
%   - Pairs are displayed as: R1–R2 => threat–lipsmack, R1–R3 => threat–chew, R2–R3 => lipsmack–chew.
%   - Significance per pair is combined across days via Fisher’s method from per-day Mantel p-values.
%
p = inputParser;
addParameter(p,'Metric','Spearman',@(s)ischar(s)||isstring(s));
addParameter(p,'RegionSrc','File',@(s)ischar(s)||isstring(s));
addParameter(p,'RegionPatterns', {'All','All'; 'S1','S1'; 'M1','M1'; 'PMv','PMv'; 'M3','M3'});
addParameter(p,'Alpha',0.05,@(x)isnumeric(x)&&isscalar(x)&&x>0&&x<1);
addParameter(p,'NBoot',5000,@(x)isnumeric(x)&&isscalar(x)&&x>=100);
addParameter(p,'ShowPoints',true,@islogical);
parse(p,varargin{:});
metric   = lower(string(p.Results.Metric));
regionSrc= char(p.Results.RegionSrc);
patMap   = p.Results.RegionPatterns;
alpha    = p.Results.Alpha;
NBoot    = p.Results.NBoot;
showPts  = p.Results.ShowPoints;

% ---- columns & labels ----
switch metric
    case "spearman", effCol='Spearman_r'; pCol='Spearman_p'; ylab='Rank-order similarity (Spearman \rho)';
    case "kendall",  effCol='Kendall_tau'; pCol='Kendall_p';  ylab='Rank-order similarity (Kendall \tau)';
    otherwise, error('Metric must be ''Spearman'' or ''Kendall''.');
end
needCols = {'Pair',effCol,pCol,regionSrc};
missing = setdiff(needCols, combined.Properties.VariableNames);
if ~isempty(missing), error('Combined table missing: %s', strjoin(missing,', ')); end

pairOrder   = {'R1-R2','R1-R3','R2-R3'};
pairDisplay = {'Threat–LS','Threat–Chew','LS–Chew'};
M = numel(pairOrder);

% ---- assign Region from pattern map (first match wins) ----
region = strings(height(combined),1);
srcVals = string(combined.(regionSrc));
for i = 1:size(patMap,1)
    pat   = string(patMap{i,1});          % pattern to search
    label = string(patMap{i,2});          % label to assign
    hit = contains(srcVals, pat, 'IgnoreCase', true);
    region(hit) = label;
end
% keep only rows that matched a known region
keep = region ~= "";
combined_region = combined(keep,:);
combined_region.Region = region(keep);

% ---- unique regions in the requested order ----
regionList = unique(string(patMap(:,2)),'stable');

% region colormap
regionCmap = [0.5 0.5 0.5
    0.4940 0.1840 0.5560	
    0.6350 0.0780 0.1840
    0.8500 0.3250 0.0980	
    0.9290 0.6940 0.1250
    ];

%  
for r = 1:numel(regionList)
    reg = regionList(r);
    T   = combined_region(strcmp(combined_region.Region, reg), :);
    if isempty(T)
        warning('No rows found for region "%s". Skipping.', reg);
        continue;
    end

    % --- collapse across days for this region (same logic as single-plot version) ---
    mu    = nan(M,1);
    CI    = nan(2,M);
    pComb = nan(M,1);
    nPer  = nan(M,1);

    for i = 1:M
        sel   = strcmp(T.Pair, pairOrder{i});
        vals  = T.(effCol)(sel);
        pvals = T.(pCol)(sel);

        vals  = vals(isfinite(vals));
        pvals = pvals(isfinite(pvals));

        nPer(i) = numel(vals);
        if isempty(vals), continue; end

        mu(i) = median(vals,'omitnan');

        if numel(vals)==1
            CI(:,i) = [vals; vals];
        else
            bs = nan(NBoot,1);
            for b = 1:NBoot
                ii = randi(numel(vals), numel(vals), 1);
                bs(b) = median(vals(ii),'omitnan');
            end
            CI(:,i) = prctile(bs,[2.5 97.5])';
        end

        if ~isempty(pvals)
            pvals = max(pvals, 1e-12);    % floor to avoid -Inf
            X  = -2*sum(log(pvals));
            df = 2*numel(pvals);
            pComb(i) = 1 - chi2cdf(X, df);
        end
    end

    % --- plotting for this region ---
    figure('Color','w','Name',sprintf('Region %s (%s)', reg, char(metric)));
    hold on;
    x = (1:M).';

    if showPts
        gcol = regionCmap(r,:);
        for i = 1:M
            sel  = strcmp(T.Pair, pairOrder{i});
            vals = T.(effCol)(sel);
            vals = vals(isfinite(vals));
            if isempty(vals), continue; end
            jitter = (rand(size(vals))-0.5)*0.15;
            plot(i + jitter, vals, '.', 'Color', gcol, 'MarkerSize', 34);
        end
    end

    mu_col  = mu(:);
    ci_low  = CI(1,:).';
    ci_high = CI(2,:).';
    yneg    = mu_col - ci_low;
    ypos    = ci_high - mu_col;

    errorbar(x, mu_col, yneg, ypos, 'o', 'LineWidth',2, 'CapSize',10, 'MarkerSize',26,'Color', gcol);

    finiteCI = isfinite(ci_low) & isfinite(ci_high);
    rngY = diff([min(ci_low(finiteCI)), max(ci_high(finiteCI))]);
    if isempty(rngY) || ~isfinite(rngY) || rngY==0, rngY = 0.1; end

    for i = 1:M
        if ~isfinite(pComb(i)), continue; end
        if     pComb(i) < 0.001, stars='***';
        elseif pComb(i) < 0.01,  stars='**';
        elseif pComb(i) < alpha, stars='*';
        else,  stars='n.s.';
        end
        yTop = ci_high(i); if ~isfinite(yTop), yTop = mu_col(i); end
        text(x(i), yTop + 0.3*rngY, stars, 'HorizontalAlignment','center','FontSize', 18, 'FontWeight','bold');
    end

    xlim([0.5 M+0.5]); ylim([-1 1])
    xticks(x); xticklabels(pairDisplay);
    ylabel(ylab);
    title({sprintf('Region: %s — Median Spearman \\rho, bootstrap 95%% CI', reg), ...
       'p via Fisher''s method'});    grid on; box off;

    annotation('textbox',[0.13 0.01 0.8 0.05], 'String', ...
        sprintf('Points = per-day estimates. Combined p via Fisher. n per pair: %s', ...
        strjoin(arrayfun(@(n)sprintf('%d',n), nPer,'UniformOutput',false), ', ')), ...
        'EdgeColor','none','HorizontalAlignment','left','FontSize',14);
    
    print_rankorder_APA(sprintf('Region %s — %s', reg, char(metric)), ...
    char(metric), pairDisplay, mu, CI, pComb);

end % end regions 
end

function plot_rankorder_summary(combined, varargin)
% plot_rankorder_summary  Collapse across days and plot summary with CIs & significance.
% Labels: R1='threat', R2='lipsmack', R3='chew'. Collapses over Day/File.

p = inputParser;
addParameter(p,'Metric','Spearman',@(s)ischar(s)||isstring(s));
addParameter(p,'Alpha',0.05,@(x)isnumeric(x)&&isscalar(x)&&x>0&&x<1);
addParameter(p,'NBoot',5000,@(x)isnumeric(x)&&isscalar(x)&&x>=100);
addParameter(p,'ShowPoints',true,@islogical);
parse(p,varargin{:});
metric = lower(string(p.Results.Metric));
alpha  = p.Results.Alpha;
NBoot  = p.Results.NBoot;
showPts= p.Results.ShowPoints;

switch metric
    case "spearman"
        effCol = 'Spearman_r'; pCol = 'Spearman_p';
        ylab   = 'Rank-order similarity (Spearman \rho)';
    case "kendall"
        effCol = 'Kendall_tau'; pCol = 'Kendall_p';
        ylab   = 'Rank-order similarity (Kendall \tau)';
    otherwise
        error('Metric must be ''Spearman'' or ''Kendall''.');
end

needCols = {'Pair',effCol,pCol};
missing = setdiff(needCols, combined.Properties.VariableNames);
if ~isempty(missing), error('Combined table missing: %s', strjoin(missing,', ')); end

pairOrder   = {'R1-R2','R1-R3','R2-R3'};
pairDisplay = {'threat–lipsmack','threat–chew','lipsmack–chew'};
M = numel(pairOrder);

mu    = nan(M,1);
CI    = nan(2,M);
pComb = nan(M,1);
nPer  = nan(M,1);

for i = 1:M
    sel   = strcmp(combined.Pair, pairOrder{i});
    vals  = combined.(effCol)(sel);
    pvals = combined.(pCol)(sel);

    vals  = vals(isfinite(vals));
    pvals = pvals(isfinite(pvals));

    nPer(i) = numel(vals);
    if isempty(vals), continue; end

    mu(i) = median(vals, 'omitnan');

    if numel(vals)==1
        CI(:,i) = [vals; vals];
    else
        bs = nan(NBoot,1);
        for b = 1:NBoot
            ii = randi(numel(vals), numel(vals), 1);
            bs(b) = median(vals(ii), 'omitnan');
        end
        CI(:,i) = prctile(bs, [2.5 97.5])';
    end

    if ~isempty(pvals)
        pvals = max(pvals, 1e-12); % floor to avoid -Inf
        X = -2 * sum(log(pvals));
        df = 2 * numel(pvals);
        pComb(i) = 1 - chi2cdf(X, df);
    end
end

figure('Color','w'); hold on;
x = (1:M).';

% (optional) per-day points
if showPts
    gcol = 0.7*[1 1 1];
    for i = 1:M
        sel  = strcmp(combined.Pair, pairOrder{i});
        vals = combined.(effCol)(sel);
        vals = vals(isfinite(vals));
        if isempty(vals), continue; end
        jitter = (rand(size(vals))-0.5)*0.15;
        plot(i + jitter, vals, '.', 'Color', gcol, 'MarkerSize', 24);
    end
end

% --- FIX: ensure column vectors for errorbar inputs ---
mu_col   = mu(:);                % Mx1
ci_low   = CI(1,:).';            % Mx1
ci_high  = CI(2,:).';            % Mx1
yneg     = mu_col - ci_low;      % Mx1
ypos     = ci_high - mu_col;     % Mx1

% Guard against negatives from NaNs/order issues (optional)
% yneg(yneg<0) = NaN; ypos(ypos<0) = NaN;

eh = errorbar(x, mu_col, yneg, ypos, 'o', 'LineWidth',3, 'CapSize',18, 'MarkerSize',24);

% significance stars
finiteCI = isfinite(ci_low) & isfinite(ci_high);
rngY = diff([min(ci_low(finiteCI)), max(ci_high(finiteCI))]);
if isempty(rngY) || ~isfinite(rngY) || rngY==0, rngY = 0.1; end
for i = 1:M
    if ~isfinite(pComb(i)), continue; end
    if pComb(i) < 0.001, stars='***';
    elseif pComb(i) < 0.01, stars='**';
    elseif pComb(i) < alpha, stars='*';
    else, stars='n.s.';
    end
    yTop = ci_high(i);
    if ~isfinite(yTop), yTop = mu_col(i); end
    text(x(i), yTop + 0.3*rngY, stars, 'HorizontalAlignment','center', 'FontSize', 16,'FontWeight','bold');
end

xlim([0.5 M+0.5]); ylim([-1 1])
xticks(x); xticklabels(pairDisplay);
ylabel(ylab);
title(sprintf('Across-day summary (%s) — Median with bootstrap 95%% CI; p via Fisher''s method', char(metric)));
grid on; box off;

annotation('textbox',[0.13 0.01 0.8 0.05], 'String', ...
    sprintf('Stars reflect combined Mantel p (Fisher). Points = per-day estimates (n per pair: %s).', ...
    strjoin(arrayfun(@(n)sprintf('%d',nPer(i)), (1:M), 'UniformOutput',false), ', ')), ...
    'EdgeColor','none','HorizontalAlignment','left','FontSize',14);

print_rankorder_APA('Across-day summary', char(metric), pairDisplay, mu, CI, pComb);


end

% ---------- helper ----------
% ========== Helper: APA-style console print for each figure ==========
function print_rankorder_APA(contextLabel, metric, pairLabels, mu, CI, pComb)
% Prints APA-style summary lines for each figure/panel:
% "<context> | pair: Median ρ = .42, 95% CI [.31, .53], p = .004"
%
% metric: 'Spearman' or 'Kendall'
% pairLabels: 1xM cellstr of display names (e.g., {'Threat–LS','Threat–Chew','LS–Chew'})
% mu:   Mx1 vector of medians (per pair)
% CI:   2xM matrix of 95% CIs (rows: [low; high])
% pComb: Mx1 vector of (combined) p-values

    if isstring(pairLabels); pairLabels = cellstr(pairLabels); end
    if isrow(mu),   mu   = mu(:);   end
    if isrow(pComb),pComb= pComb(:);end

    switch lower(metric)
        case 'spearman', sym = '\rho';
        case 'kendall',  sym = '\tau';
        otherwise,       sym = '\rho';
    end

    fprintf('\n--- APA summary: %s ---\n', contextLabel);
    for i = 1:numel(mu)
        if ~isfinite(mu(i)) || ~all(isfinite(CI(:,i)))
            fprintf('  %s: (insufficient data)\n', pairLabels{i});
            continue;
        end
        m  = mu(i);
        lo = CI(1,i);
        hi = CI(2,i);
        p  = pComb(i);

        fprintf('  %s: Median %s = %s, 95%% CI [%s, %s], %s\n', ...
            pairLabels{i}, sym, fmtNum(m), fmtNum(lo), fmtNum(hi), fmtP(p));
    end
end

% --- tiny formatters (APA-ish) ---
function s = fmtNum(x)
    if ~isfinite(x), s = 'NaN'; return; end
    s = sprintf('%.3f', x);
    s = strrep(s,'-0.','-.');     % drop leading zero for negatives
    if startsWith(s,'0.'),  s = s(2:end); end
end

function s = fmtP(p)
    if ~isfinite(p) || isnan(p)
        s = 'p = NaN'; return;
    end
    if p < 0.001
        s = 'p < .001';
    else
        s = sprintf('p = %.3f', p);
        s = strrep(s,'0.','.');
    end
end


function s = pformat(p)
    if ~isfinite(p) || isnan(p), s = 'NaN'; return; end
    if p < 1e-3
        s = sprintf('%.2e', p);
    else
        s = sprintf('%.3f', p);
    end
end
