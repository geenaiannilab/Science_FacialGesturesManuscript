
%%%%%%%% PLOTTING Figure S4A, S4B
%%%%%%%%  Euclidean distances between gesture-specific neural trajectories  
%%%%%%%%
%%%%%%%  Inter-gesture neural trajectory distances were first averaged across all pairs (Threat-Lipsmack, Threat-Chew, Lipsmack-Chew) for each day. 
%%%%%%%   This global neural trajectory distance is plotted for each region. 
%%%%%%%   Thick lines, across-day mean distance for each region, shaded bands indicate 2 SEMs across days. 
%%%%%%%   Faint shaded bands at bottom represent shuffle-null envelopes (2.5–97.5% across all shuffles). 
%%%%%%%%  
%%%%%%%%   Each subpanel shows region-specific inter-trajectory distances. 
%%%%%%%%   Within each panel, subplots show the unique gesture trajectory pairwise comparisons (from left to right, Threat vs Chew, Threat vs. Lipsmack, Lipsmack v. Chew). 
%%%%%%%%  Thick lines, across-day mean distance, shaded bands indicate 1 SEMs across days. 
%%%%%%%%  Thin lines, per day distances. Faint shaded bands at bottom represent shuffle-null envelopes (2.5–97.5% across all shuffles). 
%%%%%%%%  Time bins at which average inter-trajectory distance was significantly greater than null...
%%%%%%%%    combining per-day p-values across days (Fisher s method) and applying FDR across time, are background shaded. 
%%%%%%%%
% written GRI 250911

clear all; close all; 

data = load('matfiles/FigS4.mat');
regions = {'S1','M1','PMv','M3','All'};

allDataOut_days = data.allDataOut_days;
shuffleStats_days = data.shuffleStats_days;
shuffleNull_days = data.shuffleNull_days;
taxis = data.taxis;

% allDataOut_days = [];
% shuffleStats_days = [];
% shuffleNull_days = [];
% 

% for dd = 1:length(dates2analyze)
%     workdir = ['~/Dropbox/PhD/SUAinfo/' dates2analyze{dd} '/Data4Analysis/'];
%     data(dd) = load([workdir '/neuralTraj_All_0.2.mat']);
%     allDataOut_days{dd} = data(dd).allDataOut;
%     shuffleStats_days{dd} = data(dd).shuffleStats;
%     shuffleNull_days{dd} = data(dd).shuffleNull;
% end

colorMap = [0.4940 0.1840 0.5560;...
    0.6350 0.0780 0.1840;...
    0.8500 0.3250 0.0980;...
    0.9290 0.6940 0.1250;
    0.75 0.75 0.75];

%% summarize distances and stats across days 
summary = summarize_distance_across_days( ...
    allDataOut_days, shuffleStats_days, shuffleNull_days, regions, ...
    'Pairs', {'ThrVCh','ThrVLS','LSVCh','Average'}, ...
    'UseMaxDim', true, ...
    'AlphaFDR', 0.01);


%% do the plotting per region 
plot_distance_summary(summary, regions, taxis, ...
    'Pairs', {'ThrVCh','ThrVLS','LSVCh'}, ...
    'RegionColors', colorMap, ...
    'ShowPerDay', true, ...
    'PerDayAlpha', 0.15);   % how much to lighten per-day lines

%% plot regions against one another 
plot_average_distance_overlay(summary, regions, taxis, ...
    'PairName','Average', ...
    'RegionColors', colorMap, ...
    'MeanLineWidth', 5, ...
    'SEMAlpha', 0.28, ...
    'NullAlpha', 0.22, ...
    'LegendLocation','eastoutside', ...
    'Title','Average Distance Across Regions: mean ± 2 SEM with null bands');
%% ====================== Subfunctions ======================


function summary = summarize_distance_across_days(allDataOut_days, shuffleStats_days, shuffleNull_days, regions, varargin)
% Summarize pairwise trajectory distances across days (mean/variability)
% and combine per-day significance (Fisher/Stouffer) per time bin.
%
% Inputs
%   allDataOut_days   : 1xD cell, day-wise allDataOut (observed) structs
%   shuffleStats_days : 1xD cell, day-wise shuffleStats structs
%   shuffleNull_days  : 1xD cell, day-wise shuffleNull structs
%   regions           : cell array; each has .label (and .channels)
%
% Name/Value options (InputParser)
%   'Pairs'     : cellstr/string/char, default {'ThrVCh','ThrVLS','LSVCh','Average'}
%   'UseMaxDim' : logical, default true (uses min across-day nDimGreater90)
%   'DimIdx'    : positive scalar, default [] (required if UseMaxDim=false)
%   'AlphaFDR'  : scalar in (0,1], default 0.05
%
% Output: summary.(regionLabel).(pair) with fields:
%   .curvesPerDay   [T x D] observed curves at aligned dim
%   .meanCurve      [T x 1]
%   .sdCurve        [T x 1]
%   .semCurve       [T x 1] 
%   .ci95           [T x 2] (row-wise percentiles, NaN-safe)
%   .p_fisher       [T x 1] Fisher-combined p across days
%   .p_fisher_fdr   [T x 1] BH-FDR adjusted across time
%   .sig_fisher     [T x 1] logical (p_fisher_fdr < AlphaFDR)
%   .p_stouffer_fdr [T x 1] BH-FDR adjusted
%   .sig_stouffer   [T x 1] logical (p_stouffer_fdr < AlphaFDR)
%   .prop_sig_day   [T x 1] fraction of days significant (per-day FDR)
%   .any_sig_day    [T x 1] any day significant (per-day FDR)
%   .dd_used        scalar   PCA dim used across days for this region
%   .null_prctiles  [T x 2] aggregated null 2.5/97.5% across days (context)

% ------------- Parse inputs (InputParser) -------------
p = inputParser;
p.FunctionName = mfilename;

defaultPairs    = {'ThrVCh','ThrVLS','LSVCh','Average'}; % trajectory comparisons 
defaultUseMax   = true;
defaultDimIdx   = [];
defaultAlphaFDR = 0.05;

addParameter(p, 'Pairs',     defaultPairs,  @(x) validatePairs(x));
addParameter(p, 'UseMaxDim', defaultUseMax, @(x) islogical(x) && isscalar(x));
addParameter(p, 'DimIdx',    defaultDimIdx, @(x) isempty(x) || (isscalar(x) && isnumeric(x) && x>=1));
addParameter(p, 'AlphaFDR',  defaultAlphaFDR, @(x) isscalar(x) && isnumeric(x) && x>0 && x<=1);

parse(p, varargin{:});
Pairs     = normalizePairs(p.Results.Pairs);
UseMaxDim = p.Results.UseMaxDim;
DimIdx    = p.Results.DimIdx;
AlphaFDR  = p.Results.AlphaFDR;

D = numel(allDataOut_days);
if D == 0, error('No day-wise inputs provided.'); end

% Infer T from first available region
T = [];
for rr = 1:numel(regions)
    rlab = regions{rr};
    if isfield(allDataOut_days{1}, rlab) && isfield(allDataOut_days{1}.(rlab),'distanceAverage')
        distAvg = allDataOut_days{1}.(rlab).distanceAverage;
        if ~isempty(distAvg)
            T = size(distAvg,1);
            break;
        end
    end
end
if isempty(T), error('Could not infer T (timepoints) from inputs.'); end

summary = struct();

% ------------- Per region aggregation -------------
for rr = 1:numel(regions)
    rlab = regions{rr};

    % Choose aligned PCA dim across days (fair comparisons) 
    if UseMaxDim
        dd_days = nan(D,1);
        for d = 1:D
            if isfield(allDataOut_days{d}, rlab) && isfield(allDataOut_days{d}.(rlab),'nDimGreater90')
                dd_days(d) = allDataOut_days{d}.(rlab).nDimGreater90;
            end
        end
        dd_used = min(dd_days(~isnan(dd_days)));
        if isempty(dd_used) || dd_used < 1, dd_used = 1; end
    else
    end

    % --------- Per pair aggregation ----------
    for pi = 1:numel(Pairs)
        pairName = Pairs{pi};                 % e.g., 'ThrVCh'
        fieldName = ['distance' pairName];    

        curves       = nan(T, D);
        perDay_p     = nan(T, D);
        perDay_sig   = false(T, D);
        allNullCols  = []; % concatenated nulls across days 

        for d = 1:D
            % Observed curve at aligned dim
            if isfield(allDataOut_days{d}, rlab) && isfield(allDataOut_days{d}.(rlab), fieldName)
                distMat = allDataOut_days{d}.(rlab).(fieldName); % [T x Nd]
                if ~isempty(distMat)
                    dd_eff = min(dd_used, size(distMat,2));
                    curves(:,d) = distMat(:, dd_eff);
                end
            end

            % Per-day p and FDR sig mask 
            if isfield(shuffleStats_days{d}, rlab) && isfield(shuffleStats_days{d}.(rlab), pairName)
                st = shuffleStats_days{d}.(rlab).(pairName);
                if isfield(st,'p'),       perDay_p(:,d)   = st.p(:); end
                if isfield(st,'sig_fdr'), perDay_sig(:,d) = logical(st.sig_fdr(:)); end
            end

            % Null columns for aggregated context band
            if isfield(shuffleNull_days{d}, rlab) && isfield(shuffleNull_days{d}.(rlab), pairName)
                nullMat_d = shuffleNull_days{d}.(rlab).(pairName); % [T x nShuf_d]
                if ~isempty(nullMat_d)
                    allNullCols = [allNullCols, nullMat_d]; 
                end
            end
        end

        % --------- Across-day variability ---------
        meanCurve = mean(curves, 2, 'omitnan');
        sdCurve   = std(curves, 0, 2, 'omitnan');
        N_t       = sum(isfinite(curves), 2);
        semCurve  = sdCurve ./ max(N_t,1).^.5;

        % Row-wise 2.5/97.5 percentiles (null envelope) 
        ci95 = rowwise_prctile_nan(curves, [2.5 97.5]);

        % --------- Combine p-values across days per time bin ---------
        p_fisher   = nan(T,1);
        p_stouffer = nan(T,1);
        z_stouffer = nan(T,1);

        for t = 1:T
            pvec = perDay_p(t, :);
            pvec = pvec(isfinite(pvec) & pvec > 0 & pvec <= 1);
            if isempty(pvec), continue; end
            p_fisher(t) = fisher_combine_p(pvec);
            [z_stouffer(t), p_stouffer(t)] = stouffer_combine_p(pvec);
        end

        % FDR across time (per region × pair)
        p_fisher_fdr   = bh_fdr_vec(p_fisher,   AlphaFDR);
        p_stouffer_fdr = bh_fdr_vec(p_stouffer, AlphaFDR);
        sig_fisher     = p_fisher_fdr   < AlphaFDR;
        sig_stouffer   = p_stouffer_fdr < AlphaFDR;

        % Day-level summaries
        prop_sig_day = mean(perDay_sig, 2, 'omitnan'); % fraction of days
        any_sig_day  = any(perDay_sig, 2);

        % Aggregated null percentiles for context band
        if ~isempty(allNullCols)
            null_prct = rowwise_prctile_nan(allNullCols, [2.5 97.5]);
        else
            null_prct = nan(T,2);
        end

        % --------- Write out ---------
        S = struct();
        S.curvesPerDay   = curves;
        S.meanCurve      = meanCurve;
        S.sdCurve        = sdCurve;
        S.semCurve       = semCurve;
        S.ci95           = ci95;
        S.p_fisher       = p_fisher;
        S.p_fisher_fdr   = p_fisher_fdr;
        S.sig_fisher     = sig_fisher;
        S.z_stouffer     = z_stouffer;
        S.p_stouffer     = p_stouffer;
        S.p_stouffer_fdr = p_stouffer_fdr;
        S.sig_stouffer   = sig_stouffer;
        S.prop_sig_day   = prop_sig_day;
        S.any_sig_day    = any_sig_day;
        S.null_prctiles  = null_prct;
        S.dd_used        = dd_used;

        summary.(rlab).(pairName) = S;
    end
end
end

function plot_distance_summary(summary, regionLabels, tAxis, varargin)
% Plot across-day distance summaries per region and pair, with region-specific colors.
%
% Inputs
%   summary      : struct from summarize_distance_across_days(...)
%   regionLabels : cellstr/char/string of region names, e.g., {'M1','S1',...}
%   tAxis        : [T x 1] time vector
%
% Name/Value options
%   'Pairs'         : which pairs to show (default {'ThrVCh','ThrVLS','LSVCh'})
%   'Use'           : 'fisher' (default) or 'stouffer' for significance shading
%   'ShowNullBand'  : true/false (default true)  — aggregated null 2.5–97.5%
%   'ShowSEM'       : true/false (default true)  — mean ± SEM band
%   'ShowPerDay'    : true/false (default false) — per-day spaghetti lines
%   'TitlePrefix'   : text (default '')
%   'RegionColors'  : [R x 3] colormap for gestures 

% ---------------- options ----------------
ip = inputParser; ip.FunctionName = mfilename;
addParameter(ip,'Pairs', {'ThrVCh','ThrVLS','LSVCh'});
addParameter(ip,'Use','fisher', @(s) any(strcmpi(s,{'fisher','stouffer'})));
addParameter(ip,'ShowNullBand',true,@islogical);
addParameter(ip,'ShowSEM',true,@islogical);
addParameter(ip,'ShowPerDay',false,@islogical);
addParameter(ip,'PerDayAlpha',0.35,@(x) isscalar(x) && x>=0 && x<=1);
addParameter(ip,'LinkYAcross',false,@islogical);
addParameter(ip,'YLabel','Distance',@(s) ischar(s) || isstring(s));
addParameter(ip,'XLabel','Time',@(s) ischar(s) || isstring(s));
addParameter(ip,'TitlePrefix','',@(s) ischar(s) || isstring(s));
addParameter(ip,'RegionColors', lines(7), @(C) isnumeric(C) && size(C,2)==3);
parse(ip,varargin{:});
Pairs         = normalizeCellStr(ip.Results.Pairs);
Use           = lower(ip.Results.Use);
ShowNullBand  = ip.Results.ShowNullBand;
ShowSEM       = ip.Results.ShowSEM;
ShowPerDay    = ip.Results.ShowPerDay;
PerDayAlpha   = ip.Results.PerDayAlpha;
LinkYAcross   = ip.Results.LinkYAcross;
ylab          = char(ip.Results.YLabel);
xlab          = char(ip.Results.XLabel);
titlePrefix   = char(ip.Results.TitlePrefix);
RegionColors  = ip.Results.RegionColors;

regionLabels = normalizeCellStr(regionLabels);
T = numel(tAxis);
nP = numel(Pairs);

% -------------- per-region figures --------------
for rr = 1:numel(regionLabels)
    rlab = regionLabels{rr};
    baseCol = get_region_color(rr, RegionColors);          % region color
    semCol  = lighten_color(baseCol, 0.55);                % SEM band tint
    dayCol  = lighten_color(baseCol, min(max(PerDayAlpha,0),1)); % per-day lines tint

    fig = figure('Name',sprintf('Summary %s', rlab),'Color','w');
    fig.Position =  [626 47 1298 774];
    tl = tiledlayout(1, nP, 'Padding','compact','TileSpacing','compact');

    axList = gobjects(1,nP);
    for pi = 1:nP
        pairName = Pairs{pi};
        axList(pi) = nexttile; hold(axList(pi), 'on');

        if ~isfield(summary, rlab) || ~isfield(summary.(rlab), pairName)
            title(axList(pi), sprintf('%s %s (missing)', rlab, pairName), 'Interpreter','none');
            continue;
        end
        S = summary.(rlab).(pairName);

        % Aggregated null band (gray)
        if ShowNullBand && isfield(S,'null_prctiles') && isequal(size(S.null_prctiles), [T 2])
            fill_between(tAxis, S.null_prctiles(:,1), S.null_prctiles(:,2), [0.90 0.90 0.90], 0.8);
        end

        % Per-day spaghetti (lightened region color)
        if ShowPerDay && isfield(S,'curvesPerDay') && ~isempty(S.curvesPerDay)
            C = S.curvesPerDay;   % [T x D]
            for d = 1:size(C,2)
                plot(tAxis, C(:,d), 'LineWidth', 0.75, 'Color', dayCol);
            end
        end

        % Mean ± SEM band (tinted region color)
        if ShowSEM && isfield(S,'semCurve') && isfield(S,'meanCurve')
            fill_between(tAxis, S.meanCurve - S.semCurve, S.meanCurve + S.semCurve, semCol, 0.35);
        end

        % Mean curve (region color)
        hMean = plot(tAxis, S.meanCurve, 'LineWidth', 2, 'Color', baseCol);

        % Significance shading (based on region color)
        sigMask = false(T,1);
        if strcmp(Use,'fisher')
            if isfield(S,'sig_fisher') && ~isempty(S.sig_fisher), sigMask = logical(S.sig_fisher(:)); end
        else
            if isfield(S,'sig_stouffer') && ~isempty(S.sig_stouffer), sigMask = logical(S.sig_stouffer(:)); end
        end
        if numel(sigMask) ~= T, sigMask = pad_or_trim(sigMask, T); end

        yl = get(axList(pi),'YLim');
        shade_sig_runs(tAxis, sigMask, yl, baseCol);

        grid on;
        xlabel(axList(pi), xlab);
        ylabel(axList(pi), ylab);
        title(axList(pi), strtrim(sprintf('%s %s %s', titlePrefix, rlab, pairName)), 'Interpreter','none');
        %legend(axList(pi), hMean, {'Across-day mean'}, 'Location','best');
        uistack(hMean,'top');
    end

    % unify y-limits if requested
    if LinkYAcross
        ymins = arrayfun(@(a) a.YLim(1), axList(isgraphics(axList)));
        ymaxs = arrayfun(@(a) a.YLim(2), axList(isgraphics(axList)));
        if ~isempty(ymins)
            yl = [min(ymins), max(ymaxs)];
            set(axList(isgraphics(axList)), 'YLim', yl);
        end
    end

    title(tl, sprintf('Across-day Distance Summary — %s', rlab), 'FontWeight','bold');
end
end

function plot_average_distance_overlay(summary, regionLabels, tAxis, varargin)
% Overlay the "Average" Euclidean distance for each region:
%   - thick mean line (region color)
%   - shaded ±2 SEM (same color, translucent)
%   - shaded null band (2.5–97.5% envelope) in a greyed region color, behind
%
% Inputs
%   summary      : struct from summarize_distance_across_days(...)
%   regionLabels : cellstr/char/string, e.g., {'M1','S1','PMv','M3'}
%
% Name/Value options
%   'PairName'        : which pair to plot (default 'Average')
%   'RegionColors'    : [R x 3] colormap 
%   'MeanLineWidth'   : width for mean lines 
%   'SEMAlpha'        : face alpha for ±SEM fills 
%   'NullAlpha'       : face alpha for null band fills 

% ---- options ----
ip = inputParser; ip.FunctionName = mfilename;
addParameter(ip,'PairName','Average', @(s) ischar(s) || isstring(s));
addParameter(ip,'RegionColors', lines(7), @(C) isnumeric(C) && size(C,2)==3);
addParameter(ip,'MeanLineWidth', 2.5, @(x) isscalar(x) && x>0);
addParameter(ip,'SEMAlpha',  0.28, @(x) isscalar(x) && x>=0 && x<=1);
addParameter(ip,'NullAlpha', 0.22, @(x) isscalar(x) && x>=0 && x<=1);
addParameter(ip,'LegendLocation','eastoutside', @(s) ischar(s) || isstring(s));
addParameter(ip,'YLabel','Distance', @(s) ischar(s) || isstring(s));
addParameter(ip,'XLabel','Time',     @(s) ischar(s) || isstring(s));
addParameter(ip,'Title','Average Distance (mean ± 2 SEM) with Null Bands', @(s) ischar(s) || isstring(s));
parse(ip, varargin{:});

pairName       = char(ip.Results.PairName);
RegionColors   = ip.Results.RegionColors;
LW             = ip.Results.MeanLineWidth;
alphaSEM       = ip.Results.SEMAlpha;
alphaNull      = ip.Results.NullAlpha;
LegendLocation = char(ip.Results.LegendLocation);
ylab           = char(ip.Results.YLabel);
xlab           = char(ip.Results.XLabel);
figTitle       = char(ip.Results.Title);

regionLabels = normalizeCellStr(regionLabels);
T = numel(tAxis);

% ---- plot it ----
fig=figure('Name','Average Distance Overlay','Color','w');
fig.Position =  [626 47 1298 774];
ax = axes; hold(ax,'on');

% 1) Draw region null percentile bands FIRST (behind lines)
for rr = 1:numel(regionLabels)
    rlab = regionLabels{rr};
    if ~isfield(summary, rlab) || ~isfield(summary.(rlab), pairName), continue; end
    S = summary.(rlab).(pairName);
    if isfield(S,'null_prctiles') && ~isempty(S.null_prctiles)
        band = S.null_prctiles;
        baseCol = get_region_color(rr, RegionColors);
        bandCol = grey_tint(baseCol, 0.65);         % greyed region color
        fill_between(tAxis, band(:,1), band(:,2), bandCol, alphaNull);
    end
end

% 2) For each region: shaded ±2*SEM and thick mean line
hForLegend = gobjects(0); legNames = {};
for rr = 1:numel(regionLabels)
    rlab = regionLabels{rr};
    if ~isfield(summary, rlab) || ~isfield(summary.(rlab), pairName), continue; end
    S = summary.(rlab).(pairName);
    if ~isfield(S,'meanCurve') || ~isfield(S,'semCurve') || isempty(S.meanCurve), continue; end

    mu  = S.meanCurve(:);
    se  = S.semCurve(:);
    mu  = pad_or_trim(mu, T);
    se  = pad_or_trim(se, T);
    lo  = mu - 2*se;
    hi  = mu + 2*se;

    baseCol = get_region_color(rr, RegionColors);
    semCol  = lighten_color(baseCol, 0.55);

    % shaded ±2 SEM
    fill_between(tAxis, lo, hi, semCol, alphaSEM);

    % thick mean line (legend)
    h = plot(tAxis, mu, 'LineWidth', LW, 'Color', baseCol, 'DisplayName', rlab);
    hForLegend(end+1) = h; %#ok<AGROW>
    legNames{end+1}   = rlab; %#ok<AGROW>
end

% aesthetics 
grid(ax,'on');
xlabel(ax, xlab);
ylabel(ax, ylab);
title(ax, figTitle, 'Interpreter','none');

if ~isempty(hForLegend)
    legend(ax, hForLegend, legNames, 'Location', LegendLocation);
end

% keep null bands behind everything
uistack(findobj(ax,'Type','patch'), 'bottom');
end



% ==================== helpers (for plotting) ====================

function tf = validatePairs(P)
% Accept cellstr, string array, char row/matrix
tf = iscellstr(P) || isstring(P) || ischar(P);
end

function Pout = normalizePairs(P)
% Normalize to 1xN cell array of char
if isstring(P)
    Pout = cellstr(P);
elseif ischar(P)
    if size(P,1) == 1
        Pout = {P};
    else
        Pout = cellstr(P); % char matrix rows -> cellstr
    end
else
    Pout = P;
end
Pout = Pout(:).';
Pout = cellfun(@char, Pout, 'UniformOutput', false);
end

function p = fisher_combine_p(pvals)
% Fisher's method (right tail of chi^2 with 2k df)
k  = numel(pvals);
X2 = -2 * sum(log(pvals));
p  = 1 - chi2cdf(X2, 2*k);
end

function [z, p] = stouffer_combine_p(pvals, w)
% Stouffer's Z (one-sided). Optional weights w.
if nargin < 2 || isempty(w), w = ones(size(pvals)); end
z_i = norminv(1 - pvals);           % one-sided conversion
z   = sum(w .* z_i) / sqrt(sum(w.^2));
p   = 1 - normcdf(z);               % one-sided p
end

function p_fdr = bh_fdr_vec(p, alpha)
% Benjamini–Hochberg adjusted p-values (vector). Threshold with < alpha.
if nargin < 2 || isempty(alpha), alpha = 0.05; end
m = numel(p);
p = p(:);
[ps, idx] = sort(p);
adj = nan(m,1);
minv = 1;
for i = m:-1:1
    val = ps(i) * m / i;
    minv = min(minv, val);
    adj(i) = minv;
end
p_fdr = nan(size(p));
p_fdr(idx) = min(adj, 1);
p_fdr = reshape(p_fdr, size(p));
end

function P = rowwise_prctile_nan(X, pr)
% Row-wise percentiles ignoring NaNs. Returns [size(X,1) x numel(pr)].
T = size(X,1);
P = nan(T, numel(pr));
for t = 1:T
    x = X(t, :);
    x = x(isfinite(x));
    if ~isempty(x)
        P(t, :) = prctile(x, pr);
    end
end
end

function C = normalizeCellStr(C)
if isstring(C), C = cellstr(C); end
if ischar(C)
    if size(C,1)==1, C = {C}; else, C = cellstr(C); end
end
C = C(:).';
C = cellfun(@char, C, 'UniformOutput', false);
end

function fill_between(x, y1, y2, faceColor, alphaVal)
x = x(:); y1 = y1(:); y2 = y2(:);
X = [x; flipud(x)];
Y = [y1; flipud(y2)];
patch('XData',X,'YData',Y, ...
      'FaceColor',faceColor, 'FaceAlpha',alphaVal, ...
      'EdgeColor','none');
end
function shade_sig_runs(tAxis, sigMask, ylims, baseColor)
if ~any(sigMask), return; end
tint = 0.2 + 0.8*baseColor; tint(tint>1)=1;
dSig = diff([false; sigMask(:); false]);
runStarts = find(dSig==1); runEnds = find(dSig==-1)-1;
for k = 1:numel(runStarts)
    i0 = runStarts(k); i1 = runEnds(k);
    patch('XData',[tAxis(i0) tAxis(i1) tAxis(i1) tAxis(i0)], ...
          'YData',[ylims(1) ylims(1) ylims(2) ylims(2)], ...
          'FaceColor',tint, 'FaceAlpha',0.15, 'EdgeColor','none');
end
uistack(findobj(gca,'Type','patch'),'bottom');
end

function v = pad_or_trim(v, T)
v = v(:);
if numel(v) < T
    v = [v; false(T - numel(v), 1)];
elseif numel(v) > T
    v = v(1:T);
end
end

function col = get_region_color(rr, Cmap)
% Map region index rr -> row in colormap, cycling if needed
nC = size(Cmap,1);
idx = mod(rr-1, nC) + 1;
col = Cmap(idx, :);
col = max(min(col,1),0);
end

function out = lighten_color(col, amount)
% Mix color with white by 'amount' (0=no change, 1=white)
amount = min(max(amount,0),1);
out = (1-amount)*col + amount*[1 1 1];
end

function out = grey_tint(col, greyAmt)
% Blend toward neutral grey by greyAmt (0=no change, 1=grey)
% "Grey" here is mid-grey [0.5 0.5 0.5] for a subdued band color.
greyAmt = min(max(greyAmt,0),1);
out = (1-greyAmt)*col + greyAmt*[0.5 0.5 0.5];
end
