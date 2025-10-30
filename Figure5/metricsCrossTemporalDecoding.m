
%%%% PLOTTING Fig S6 (also will plot Fig 5, left and middle panels) 
%%%% Generates summary metrics regarding stability/dynamnicity of codes 
%%%%
%%%% DI = Fraction of decoding "mass" on/near the diagonal vs off-diagonal.
%%%% Dynamic ⇢ high DI (narrow diagonal ribbon). Stable ⇢ lower DI (broad off-diagonal generalization).
%%%%
%%%% Temporal generalization width (TGW): For each train time, find the full-width at half-max (FWHM) of accuracy across test times; 
%%%% report the median width. Dynamic ⇢ narrow TGW.

% written GI 250831

close all; clear all ;
set(0,'defaultAxesFontSize',24)
set(0,'defaultAxesFontWeight','bold')

%% intial set-up
subject = 'combined';
nIterations = 50; nUnitsPerRegion = 50;
plotIndividualIterations = false; 
binSize = 0.4;
colorArray = [0.4940 0.1840 0.5560;0.6350 0.0780 0.1840;0.8500 0.3250 0.0980;0.9290 0.6940 0.1250];

%% load data 
f2load = dir('matfiles/Fig5_left.mat');
ctdat = load([f2load.folder '/' f2load.name]);

%% Create timebase, channels 
arrayList = ctdat.arrayList;
windowsCenter = ctdat.windowsCenter; 
time2test = -0.8:0.4:0.8; 
skipFactor = 10;


%% ===============================
%% calculate metrics : DI and TGW
%% ===============================

chance = 0.6;
band   = 10;                  % half-width (in bins) of the diagonal strip you count as "on/near the diagonal3 when computing DI;                               % pick band roughly matching your temporal resolution or expected generalization (e.g., 1–2 bins for sharp, dynamic codes
                              % larger band → more matrix “mass” counted as diagonal ⇒ higher DI; smaller band makes the metric stricter (more “dynamic” if mass is tightly on the diagonal).
                   
out = summarize_ctd_dynamics(ctdat, arrayList, ...
     'chance',chance, 'band',band, 'tAxis',windowsCenter,'Colors',colorArray);

%% DI comparison 
% Build per-iteration DI vectors (set chance and diagonal band as you use in figures)
diCell = build_di_vectors(ctdat, arrayList, 'Chance', chance, 'Band', band);

% Run all pairwise, unpaired permutation comparisons
resDI = compare_DI_unpaired_all(diCell, arrayList, 'Seed', 123);

%% TGW comparison
% Build per-iteration TGW vectors (choose units via tAxis)
tgwCell = build_tgw_vectors(ctdat, arrayList, 'Chance', chance, 'TAxis', windowsCenter);

% Run all pairwise comparisons (unpaired permutation)
resTGW = compare_TGW_unpaired_all(tgwCell, arrayList, ...
           'Seed', 123);


%% Plot pairwise comparisons (positive = A more dynamic)
% Plot with FDR correction
outDI = plot_pairwise_metric_matrix(resDI, ...
        'MetricLabel','DI (higher = more dynamic)', ...
        'Alpha',0.05, 'Adjust','fdr', 'ShowValues',true);

outTGW = plot_pairwise_TGW_matrix(resTGW, 'Alpha',0.05, 'Adjust','fdr', ...
                                  'FlipForDynamic', true, 'ShowValues', true);
%% Plot all results per region
% Summary bars with overlaid brackets (FDR across pairs; use two-sided p)
summarize_ctd_dynamics(ctdat, arrayList, ...
   'chance', chance, 'band', band, 'tAxis', windowsCenter, 'Colors', colorArray, ...
   'resDI', resDI, 'resTGW', resTGW, ...
   'Alpha', 0.05, 'Adjust', 'fdr', 'UsePColumn', 'p_two', ...
   'ShowStars', true);

%% Plot the actual CTDs
for aa = 1:length(arrayList)

    data = mean(ctdat.validationScoreStruct.(arrayList{aa}),3);
    sem = std(ctdat.validationScoreStruct.(arrayList{aa}),[],3) ./ sqrt(nIterations);
    cmax = max(sem,[],'all');

    diagData = diag(data);
    diagSEM = diag(sem);

    plotFig = figure('Position',[ 476    86   964   780]);
    zeroPos = find(windowsCenter == 0);
    imagesc(data); hold on; 
    line([zeroPos zeroPos], [1 length(windowsCenter)],'Color', 'w', 'LineStyle',':','LineWidth',2.5); 
    line([1 length(windowsCenter)],[zeroPos zeroPos], 'Color', 'w', 'LineStyle',':','LineWidth',2.5); 

    hold off
    set ( gca, 'ydir', 'normal' );
    xlabel('Generalization Time')
    ylabel('Training Time')
    xticks(1:skipFactor:(length(windowsCenter)))
    xticklabels(windowsCenter(xticks))
    yticks(1:skipFactor:(length(windowsCenter)))
    yticklabels(windowsCenter(xticks))
    
    colorbar; colorbarpzn(0.5,1,'full',0.6); %colormap(jet); caxis([0.5 1]);%
    set ( gca, 'FontSize', 24);
    title([arrayList{aa} ' - Generalization Across Time'], 'FontSize',28)
 
    plotFig = figure('units','normalized','outerposition',[  0    0.0886    0.2817    1]);
    nSubplots = length(time2test);
    ha = tight_subplot(nSubplots,1, 0.02, [0.05 0.03], [0.12 0.1]);

    for tt = 1:length(time2test)
        t2testIdx = find(windowsCenter >= time2test(tt),1,'first');
        t2testPerform = data(t2testIdx,:);
        t2testErr = sem(t2testIdx,:);
        
        axes(ha((length(time2test)-tt)+1))
        errorbar(windowsCenter, diagData, diagSEM,'k','linew',2); hold on; 
        errorbar(windowsCenter,t2testPerform, t2testErr,'Color', [0 0.4470 0.7410],'linew',2); hold off
        xlim([-0.8 0.8])
        str = {['Train = ' num2str(time2test(tt)) ' s']}; 
        yyaxis right ;ylabel(str); tmp = gca; tmp.YTick = []; tmp.YLabel.Color = 'k'; tmp.FontSize = 18; tmp.FontWeight = 'bold';
       
        if tt > 1
            tmp.XTickLabel = '';
        end
        box off
    end
    
    han=axes(plotFig,'visible','off'); 
    han.Title.Visible='on'; 
    han.XLabel.Visible='on';
    han.YLabel.Visible='on';
    yyaxis right; 
    yyaxis left; ylabel(han,'Accuracy');han.YLabel.Color = 'k';han.YLabel.FontSize = 20;
    box off; 
    han.YLabel.Position = [-0.1 0.5 0 ];

    xlabel(han,'Test Time Bin','FontSize',20);
    han.XLabel.Color = 'k'; 
    han.XLabel.Position = [0.5 -0.1 0];

    han.Title.Position = [0.5 1.06 0.5];    
    title(han,'Decoder Performance','FontSize',20);
    
    plotFig = gcf;
    plotFig.Color = [1 1 1 ];

end


function out = summarize_ctd_dynamics(ctdat, arrayList, varargin)
% Compute DI, TGW for each region and plot bar charts with optional
% pairwise significance overlays from compare_*_unpaired_all results.
%
% Name-Value options:
%   'chance'     : scalar chance level (defined by CTD) 
%   'band'       : diagonal half-width for DI (bins; default 1)
%   'tAxis'      : time vector for TGW 
%   'resDI'      : struct from compare_DI_unpaired_all (default [])
%   'resTGW'     : struct from compare_TGW_unpaired_all (default [])
%   'Alpha'      : alpha for significance (default 0.05)
%   'Adjust'     : 'none' | 'fdr' | 'bonferroni' (default 'fdr')
%   'UsePColumn' : which p-column from res.table ('p_two','p_right','p_left')

p = inputParser;
addParameter(p,'chance',0,@isnumeric);
addParameter(p,'band',1,@(x)isnumeric(x)&&isscalar(x));
addParameter(p,'tAxis',[],@isnumeric);
addParameter(p,'Colors',[],@isnumeric);
addParameter(p,'resDI',[],@(s)isstruct(s)||isempty(s));
addParameter(p,'resTGW',[],@(s)isstruct(s)||isempty(s));
addParameter(p,'Alpha',0.05,@(x)isnumeric(x)&&isscalar(x)&&x>0&&x<1);
addParameter(p,'Adjust','fdr',@(s)any(strcmpi(string(s),["none","fdr","bonferroni"])));
addParameter(p,'UsePColumn','p_two',@(s)any(strcmpi(string(s),["p_two","p_right","p_left"])));
addParameter(p,'ShowStars',true,@islogical);
addParameter(p,'MaxPairs',inf,@(x)isnumeric(x)&&isscalar(x));
parse(p,varargin{:});
chance = p.Results.chance;
band   = p.Results.band;
tAxis  = p.Results.tAxis;
C      = p.Results.Colors;
resDI  = p.Results.resDI;
resTGW = p.Results.resTGW;
alpha  = p.Results.Alpha;
adjust = lower(string(p.Results.Adjust));
pcol   = string(p.Results.UsePColumn);
showStars = p.Results.ShowStars;
maxPairs  = p.Results.MaxPairs;

R = numel(arrayList);
if isempty(C), C = lines(R); end

% ---------- Compute metrics ----------
DI = nan(R,1); DIsem = nan(R,1);
TGW = nan(R,1); TGWsem = nan(R,1);

for r = 1:R
    CTD = ctdat.validationScoreStruct.(arrayList{r}); % T x T x N
    di  = diagonality_index(CTD, 'band', band, 'chance', chance);
    tgw = temporal_generalization_width(CTD, 'chance', chance,'tAxis', tAxis);
    DI(r)    = di.mean;   DIsem(r)  = di.sem;
    TGW(r)   = tgw.mean;  TGWsem(r) = tgw.sem;
end

% ---------- Build pairwise p-matrices (optional) ----------
P_DI = []; pairs_DI = []; % [i j p_adj]
if ~isempty(resDI)
    [P_DI, pairs_DI] = build_pmatrix_from_result(resDI, arrayList, pcol, adjust);
end
P_TGW = []; pairs_TGW = []; % [i j p_adj]
if ~isempty(resTGW)
    [P_TGW, pairs_TGW] = build_pmatrix_from_result(resTGW, arrayList, pcol, adjust);
end

% ---------- Plot ----------
figure('Color','w','Name','CTD Dynamics Summary');

% ---- DI subplot ----
ax1 = subplot(1,2,1); hold(ax1,'on');
h1  = bar(ax1, 1:R, DI, 'FaceColor','flat', 'EdgeColor','none');
h1.CData = C;
for r = 1:R
    errorbar(ax1, r, DI(r), DIsem(r), 'k.', 'LineWidth', 2, 'CapSize', 3);
end
xticks(ax1, 1:R); xticklabels(ax1, arrayList); xtickangle(ax1, 30);
ylabel(ax1, 'Diagonality Index');
title(ax1, sprintf('Diagonality Index'));
box(ax1,'off'); 

% Overlay DI pairwise brackets
if showStars && ~isempty(P_DI)
    overlay_pairwise_brackets(ax1, DI, DIsem, P_DI, alpha, maxPairs);
end

% ---- TGW subplot ----
ax2 = subplot(1,2,2); hold(ax2,'on');
h2  = bar(ax2, 1:R, TGW, 'FaceColor','flat', 'EdgeColor','none');
h2.CData = C;
for r = 1:R
    errorbar(ax2, r, TGW(r), TGWsem(r), 'k.', 'LineWidth', 2, 'CapSize', 3);
end
xticks(ax2, 1:R); xticklabels(ax2, arrayList); xtickangle(ax2, 30);
if isempty(tAxis), ylabel(ax2, 'Width (bins)'); else, ylabel(ax2, 'Width (seconds)'); end
title(ax2, 'Temporal Generalization Width');
box(ax2,'off'); 

% Overlay TGW pairwise brackets
if showStars && ~isempty(P_TGW)
    overlay_pairwise_brackets(ax2, TGW, TGWsem, P_TGW, alpha, maxPairs);
end

% ---------- Pack outputs ----------
out.DI  = table(arrayList(:), DI, DIsem, 'VariableNames', {'Region','DI_mean','DI_sem'});
out.TGW = table(arrayList(:), TGW, TGWsem, 'VariableNames', {'Region','TGW_mean','TGW_sem'});
out.P_DI  = P_DI;   % R x R adjusted p-values 
out.P_TGW = P_TGW;  

end

% ===================== helpers =====================

function [P, pairs] = build_pmatrix_from_result(res, arrayList, pcol, adjust)
% Build symmetric R x R matrix of adjusted p-values using desired p column.
T = res.table;
regs = string(arrayList(:));
R = numel(regs);
idx = containers.Map(regs, num2cell(1:R));
P = nan(R);    % symmetric (NaN on diag)
raw = nan(R);  % store upper triangle raw ps to adjust together

% collect upper-tri p's in consistent region order
maskUT = false(R); 
plist = [];  ij = [];
for k = 1:height(T)
    A = string(T.RegionA(k)); B = string(T.RegionB(k));
    if ~(isKey(idx,A) && isKey(idx,B)), continue; end
    i = idx(A); j = idx(B);
    if i==j, continue; end
    if i>j, tmp=i; i=j; j=tmp; end
    plist(end+1,1) = T.(pcol)(k); 
    ij(end+1,:) = [i j]; 
end

% adjust across all collected pairs
padj = adjust_p(plist, adjust);
for e = 1:numel(padj)
    i = ij(e,1); j = ij(e,2);
    P(i,j) = padj(e);
    P(j,i) = padj(e);
end
pairs = ij; % indices of pairs used
end

function di = diagonality_index(CTD, varargin)
% DI = mass near diagonal / total mass (after baseline subtracting)
% CTD: T x T x N iterations; options: 'band' (bins), 'chance' (scalar)

p = inputParser;
addParameter(p,'band',1,@(x)isnumeric(x)&&isscalar(x));
addParameter(p,'chance',0,@isnumeric);
parse(p,varargin{:});
w = p.Results.band; chance = p.Results.chance;

A = max(CTD - chance, 0);                    % baseline subtract
T = size(A,1);
diagMask = abs((1:T)' - (1:T)) <= w;         
tot = squeeze(sum(sum(A,1),2));             
on  = squeeze(sum(sum(A .* diagMask,1),2)); 
di.iter = on ./ max(tot, eps);
di.mean = mean(di.iter,'omitnan');
di.sem  = std(di.iter,0,'omitnan')/sqrt(numel(di.iter));
end

function tgw = temporal_generalization_width(CTD, varargin)
% For each train time, FWHM across test times; return median per iter.
% Options: 'chance' (scalar), 'tAxis' (time vector, else bins)

p = inputParser;
addParameter(p,'chance',0,@isnumeric);
addParameter(p,'tAxis',[],@isnumeric);
parse(p,varargin{:});
chance = p.Results.chance; tAxis = p.Results.tAxis;

[T,~,N] = size(CTD);
widths = nan(N,1);

for n = 1:N 
    A = squeeze(CTD(:,:,n)) - chance; A(A<0)=0;
    w_i = nan(T,1);

    for i = 1:T % time points
        row = A(i,:);
        if all(~isfinite(row)) || max(row)==0, continue; end
        half = 0.5*max(row);
        above = row >= half;
        if ~any(above), continue; end
        idx = find(above);
        i0 = idx(1); i1 = idx(end);
        if isempty(tAxis)
            w_i(i) = i1 - i0 + 1;                         % bins
        else
            w_i(i) = tAxis(i1) - tAxis(i0);              % time units
        end
    end
    widths(n) = median(w_i,'omitnan');
end % end iterations 
tgw.iter = widths;
tgw.mean = mean(widths,'omitnan'); % over iterations 
tgw.sem  = std(widths,0,'omitnan')/sqrt(N); % over iterations 
end

function out = compare_DI_unpaired_all(diCell, names, varargin)
% Run unpaired label-shuffle permutation for DI across all 2-way region pairs.
% INPUTS
%   diCell : R-by-1 cell, each a vector of DI values (per iteration) for a region
%   names  : R-by-1 cellstr/string array of region names

p = inputParser;
addParameter(p,'NPerm',10000,@(x)isnumeric(x)&&isscalar(x)&&x>=100);
addParameter(p,'NBoot',5000,@(x)isnumeric(x)&&isscalar(x)&&x>=500);
addParameter(p,'Seed',[],@(x)isempty(x)||isscalar(x));
parse(p,varargin{:});
nPerm = p.Results.NPerm; nBoot = p.Results.NBoot; seed = p.Results.Seed;

if isstring(names) || ischar(names), names = cellstr(names); end
R = numel(diCell);
pairs = nchoosek(1:R,2);

rows = [];
pairStats = struct([]);
for k = 1:size(pairs,1)
    i = pairs(k,1); j = pairs(k,2);
    mA = diCell{i}; mB = diCell{j};

    stats = permtest_unpaired_mean_diff(mA, mB, ...
             'NPerm', nPerm, 'NBoot', nBoot, 'Seed', seed);

    rows = [rows; {names{i}, names{j}, ...
                   mean(mA,'omitnan'), mean(mB,'omitnan'), stats.diff, ...
                   stats.p_two, stats.p_right, stats.p_left, ...
                   stats.nA, stats.nB, ...
                   stats.ciA(1), stats.ciA(2), stats.ciB(1), stats.ciB(2), ...
                   stats.ciDiff(1), stats.ciDiff(2)}];

    pairStats(k).RegionA = names{i};
    pairStats(k).RegionB = names{j};
    pairStats(k).stats   = stats;
end

T = cell2table(rows, 'VariableNames', { ...
    'RegionA','RegionB','MeanA','MeanB','Diff_AminusB', ...
    'p_two','p_right','p_left','nA','nB', ...
    'CI_A_low','CI_A_high','CI_B_low','CI_B_high','CI_Diff_low','CI_Diff_high'});

out = struct();
out.table = T;            % tidy summary for all pairs
out.pairs = pairs;        % indices into 'names'
out.stats = pairStats;    % raw stats structs per pair
end

function out = compare_TGW_unpaired_all(tgwCell, names, varargin)
% Run unpaired permutation comparisons of TGW for all 2-way region pairs.
% INPUTS
%   tgwCell : R-by-1 cell, each a vector of TGW values (per iteration) for a region
%   names   : R-by-1 cellstr/string array of region names

p = inputParser;
addParameter(p,'NPerm',10000,@(x)isnumeric(x)&&isscalar(x)&&x>=100);
addParameter(p,'NBoot',5000,@(x)isnumeric(x)&&isscalar(x)&&x>=500);
addParameter(p,'Seed',[],@(x)isempty(x)||isscalar(x));
parse(p,varargin{:});
nPerm = p.Results.NPerm; nBoot = p.Results.NBoot; seed = p.Results.Seed;

if isstring(names) || ischar(names), names = cellstr(names); end
R = numel(tgwCell);
pairs = nchoosek(1:R,2);

rows = [];
pairStats = struct([]);
for k = 1:size(pairs,1)
    i = pairs(k,1); j = pairs(k,2);
    mA = tgwCell{i}; mB = tgwCell{j};

    % Unpaired permutation on mean difference (A - B)
    stats = permtest_unpaired_mean_diff(mA, mB, ...
             'NPerm', nPerm, 'NBoot', nBoot, 'Seed', seed);

    % Assemble a  row
    rows = [rows; {names{i}, names{j}, ...
                   mean(mA,'omitnan'), mean(mB,'omitnan'), stats.diff, ...
                   stats.p_two, stats.p_right, stats.p_left, ...
                   stats.nA, stats.nB, ...
                   stats.ciA(1), stats.ciA(2), stats.ciB(1), stats.ciB(2), ...
                   stats.ciDiff(1), stats.ciDiff(2)}]; %#ok<AGROW>

    pairStats(k).RegionA = names{i};
    pairStats(k).RegionB = names{j};
    pairStats(k).stats   = stats;
end

T = cell2table(rows, 'VariableNames', { ...
    'RegionA','RegionB','MeanA','MeanB','Diff_AminusB', ...
    'p_two','p_right','p_left','nA','nB', ...
    'CI_A_low','CI_A_high','CI_B_low','CI_B_high','CI_Diff_low','CI_Diff_high'});

out = struct();
out.table = T;            
out.pairs = pairs;        
out.stats = pairStats;  
end

function diCell = build_di_vectors(ctdat, arrayList, varargin)
% diCell{r} = vector of DI values (one per iteration) for region r
% ctdat.validationScoreStruct.(arrayList{r}) is T x T x N (accuracy)
% Options: 'Chance' (scalar), 'Band' (half-width in bins)

p = inputParser;
addParameter(p,'Chance',0,@isnumeric);
addParameter(p,'Band',1,@(x)isnumeric(x)&&isscalar(x)&&x>=0);
parse(p,varargin{:});
chance = p.Results.Chance;
band   = p.Results.Band;

R = numel(arrayList);
diCell = cell(R,1);
for r = 1:R
    CTD = ctdat.validationScoreStruct.(arrayList{r});   % T x T x N
    N   = size(CTD,3);
    m   = nan(N,1);
    for n = 1:N
        di = diagonality_index(CTD(:,:,n), 'band', band, 'chance', chance);
        m(n) = di.mean;   % scalar for that iteration
    end
    diCell{r} = m;
end
end

function tgwCell = build_tgw_vectors(ctdat, arrayList, varargin)
% tgwCell{r} = vector of TGW values (one per iteration) for region r
% ctdat.validationScoreStruct.(arrayList{r}) is T x T x N (accuracy)
% Options: 'Chance' (scalar), 'TAxis' (vector; [] -> widths in bins)

p = inputParser;
addParameter(p,'Chance',0,@isnumeric);
addParameter(p,'TAxis',[],@isnumeric);
parse(p,varargin{:});
chance = p.Results.Chance;
tAxis  = p.Results.TAxis;

R = numel(arrayList);
tgwCell = cell(R,1);
for r = 1:R
    CTD = ctdat.validationScoreStruct.(arrayList{r});   % T x T x N
    tgw = temporal_generalization_width(CTD, 'chance', chance, 'tAxis', tAxis);
    % tgw.iter is N×1, where each entry is the median FWHM across train times
    tgwCell{r} = tgw.iter(:);
end
end

function out = plot_pairwise_metric_matrix(res, varargin)
% Plot pairwise mean differences and significance from compare_*_unpaired_all().
%
% INPUT
%   res.table must contain columns:
%     RegionA, RegionB, MeanA, MeanB, Diff_AminusB, p_two  (plus CIs, etc.)
%
% OPTIONS
%   'MetricLabel' : e.g., 'DI (higher = more dynamic)'  
%   'Alpha'       : significance level (default 0.05)
%   'Adjust'      : 'none' | 'fdr' | 'bonferroni' (default 'none')
%
% OUTPUT (struct)
%   .regions  : cellstr of region names (plot order)
%   .DiffMat  : R x R matrix with A-B mean differences (diag = 0)
%   .PvalMat  : R x R matrix with (adjusted) p-values (two-sided)
%   .SigMat   : logical matrix where p < alpha
%   .p_table  : table of pairs with raw & adjusted p-values

p = inputParser;
addParameter(p,'MetricLabel','Metric',@(s)ischar(s)||isstring(s));
addParameter(p,'Alpha',0.05,@(x)isnumeric(x)&&isscalar(x)&&x>0&&x<1);
addParameter(p,'Adjust','none',@(s)any(strcmpi(string(s),["none","fdr","bonferroni"])));
addParameter(p,'Colormap',[],@(c)ischar(c)||isstring(c)||ismatrix(c));
addParameter(p,'ShowValues',true,@islogical);
parse(p,varargin{:});
metricLabel = string(p.Results.MetricLabel);
alpha = p.Results.Alpha;
adj   = lower(string(p.Results.Adjust));
cmap  = p.Results.Colormap;
showVals = p.Results.ShowValues;

T = res.table;
regs = unique([string(T.RegionA); string(T.RegionB)],'stable');
[regs, ord] = sort(regs);
R = numel(regs);

% Map region name -> index
idxOf = containers.Map(regs, num2cell(1:R));

% Build raw matrices
Diff = zeros(R); P = nan(R);
for k = 1:height(T)
    A = string(T.RegionA(k)); B = string(T.RegionB(k));
    i = idxOf(A); j = idxOf(B);
    d = T.Diff_AminusB(k);
    p2 = T.p_two(k);
    Diff(i,j) = d;
    Diff(j,i) = -d;
    P(i,j) = p2; P(j,i) = p2;
end

% Adjust p-values across all pairs (upper triangle only)
tri = find(triu(true(R),1));
raw = P(tri);
switch adj
    case "fdr"
        padj = fdr_bh(raw);
    case "bonferroni"
        padj = min(raw * numel(raw), 1);
    otherwise
        padj = raw;
end
P_adj = P;
P_adj(tri) = padj; P_adj = min(P_adj,1); P_adj = max(P_adj,0); % mirror
P_adj = P_adj + tril(P_adj.',-1); % symmetric; keep NaN on diag
Sig = false(R); Sig(tri) = padj < alpha; Sig = Sig | Sig.'; % no diag

% Colormap (diverging, centered at 0)
if isempty(cmap)
    % blue-white-red
    n=256; r = (0:n-1)'/(n-1);
    blu = [linspace(0,1,n/2)', linspace(0,1,n/2)', ones(n/2,1)];
    red = [ones(n/2,1), linspace(1,0,n/2)', linspace(1,0,n/2)'];
    cmap = [blu; red];
end

% Plot
figure('Color','w','Name',sprintf('Pairwise %s comparisons', metricLabel));
ax = axes; hold(ax,'on');

% Symmetric color limits around 0
mx = max(abs(Diff(:)));  if mx==0, mx = 1; end
imagesc(ax, Diff, [-mx mx]); colormap(ax, cmap);
cb = colorbar(ax); cb.Label.String = sprintf('Mean difference (A - B) of %s', metricLabel);
axis(ax,'equal','tight');
xticks(1:R); yticks(1:R);
xticklabels(regs); yticklabels(regs);
xtickangle(45);
title(ax, sprintf('Pairwise differences — %s | p<%.3g (%s)', metricLabel, alpha, upper(adj)));
set(ax,'TickLength',[0 0]); box(ax,'on');

% Zero diagonal line
for i=1:R, rectangle('Position',[i-0.5 i-0.5 1 1],'EdgeColor',[0.8 0.8 0.8]); end

% Annotations: numbers and significance stars (upper triangle)
for i=1:R
    for j=1:R
        if i==j, continue; end
        val = Diff(i,j);
        if showVals
            text(j, i, sprintf('%+.2f', val), 'HorizontalAlignment','center', ...
                 'VerticalAlignment','middle', 'FontSize',20, 'Color','k', ...
                 'FontWeight', 'bold', 'Interpreter','none');
        end
        if Sig(i,j)
            % Significance stars based on adjusted p
            p_ij = P_adj(i,j);
            stars = p_to_stars(p_ij);
            text(j, i, stars, 'HorizontalAlignment','right', ...
                 'VerticalAlignment','top', 'FontSize',24, 'Color','k', 'FontWeight','bold');
        end
    end
end

% Output
out = struct();
out.regions = cellstr(regs);
out.DiffMat = Diff;
out.PvalMat = P_adj;
out.SigMat  = Sig;

% adjusted p-values
pairs = nchoosek(1:R,2);
rows = cell(size(pairs,1), 6);
for k = 1:size(pairs,1)
    i = pairs(k,1); j = pairs(k,2);
    rows(k,:) = {regs(i), regs(j), Diff(i,j), P(i,j), P_adj(i,j), Sig(i,j)};
end
out.p_table = cell2table(rows, 'VariableNames', ...
    {'RegionA','RegionB','Diff_AminusB','p_two_raw','p_adj','is_sig'});

end

function out = plot_pairwise_TGW_matrix(res, varargin)
% Plot pairwise mean differences and significance for TGW results.
% INPUT: res = compare_TGW_unpaired_all(...) output
%
% OPTIONS
%   'Alpha'          : significance level (default 0.05)
%   'Adjust'         : 'none' | 'fdr' | 'bonferroni' (default 'none')
%   'Colormap'       : name or Mx3 (default diverging blue-red)
%   'ShowValues'     : true/false (default true)  % print A−B in each cell
%   'FlipForDynamic' : true/false (default false) % if true, plot (B−A) so + means A more dynamic
%
% OUTPUT
%   .regions, .DiffMat (A−B raw), .PvalMat (adj), .SigMat, .p_table

p = inputParser;
addParameter(p,'Alpha',0.05,@(x)isnumeric(x)&&isscalar(x)&&x>0&&x<1);
addParameter(p,'Adjust','none',@(s)any(strcmpi(string(s),["none","fdr","bonferroni"])));
addParameter(p,'Colormap',[],@(c)ischar(c)||isstring(c)||ismatrix(c));
addParameter(p,'ShowValues',true,@islogical);
addParameter(p,'FlipForDynamic',false,@islogical);
parse(p,varargin{:});
alpha = p.Results.Alpha;
adj   = lower(string(p.Results.Adjust));
cmap  = p.Results.Colormap;
showVals = p.Results.ShowValues;
flipDyn  = p.Results.FlipForDynamic;

T = res.table;
regs = unique([string(T.RegionA); string(T.RegionB)],'stable');
[regs,~] = sort(regs);
R = numel(regs);
idxOf = containers.Map(regs, num2cell(1:R));

% Build matrices: raw differences (A−B) and two-sided p
Diff = zeros(R); P = nan(R);
for k = 1:height(T)
    i = idxOf(string(T.RegionA(k)));
    j = idxOf(string(T.RegionB(k)));
    d = T.Diff_AminusB(k); p2 = T.p_two(k);
    Diff(i,j)=d; Diff(j,i)=-d; P(i,j)=p2; P(j,i)=p2;
end

% Adjust p-values across all pairs (upper triangle)
tri = find(triu(true(R),1));
raw = P(tri);
switch adj
    case "fdr",        padj = fdr_bh(raw);
    case "bonferroni", padj = min(raw*numel(raw), 1);
    otherwise,         padj = raw;
end
P_adj = P; P_adj(tri)=padj; P_adj = P_adj + tril(P_adj.',-1);

% What to plot: either raw A−B, or flipped so + means “A more dynamic”
PlotDiff = Diff;
if flipDyn
    PlotDiff = -Diff; % since smaller TGW = more dynamic
end

% Colormap
if isempty(cmap)
    n=256;
    blu = [linspace(0,1,n/2)', linspace(0,1,n/2)', ones(n/2,1)];
    red = [ones(n/2,1), linspace(1,0,n/2)', linspace(1,0,n/2)'];
    cmap = [blu; red];
end

% Plot
figure('Color','w','Name','Pairwise TGW comparisons');
ax = axes; hold(ax,'on');
mx = max(abs(PlotDiff(:))); if mx==0, mx=1; end
imagesc(ax, PlotDiff, [-mx mx]); colormap(ax, cmap);
cb = colorbar(ax);
if flipDyn
    cb.Label.String = 'Signed diff (B−A) so + = A more dynamic (TGW smaller)';
else
    cb.Label.String = 'Mean difference (A−B) TGW';
end
axis(ax,'equal','tight'); set(ax,'TickLength',[0 0]); box(ax,'on');
xticks(1:R); yticks(1:R); xticklabels(regs); yticklabels(regs); xtickangle(45);

title(ax, sprintf('Pairwise TGW (smaller = more dynamic) | p<%.3g (%s)', alpha, upper(adj)));

% grid boxes
for i=1:R, rectangle('Position',[i-0.5 i-0.5 1 1],'EdgeColor',[0.85 0.85 0.85]); end

% Significance stars + numeric values
Sig = false(R); Sig(tri) = (padj < alpha); Sig = Sig | Sig.';
for i=1:R
    for j=1:R
        if i==j, continue; end
        if showVals
            text(j,i,sprintf('%+.2f',Diff(i,j)), 'HorizontalAlignment','center', ...
                 'VerticalAlignment','middle','FontSize',20,'FontWeight','bold','Color','k');
        end
        if Sig(i,j)
            stars = p_to_stars(P_adj(i,j));
            text(j,i,stars,'HorizontalAlignment','right','VerticalAlignment','top', ...
                 'FontSize',24,'FontWeight','bold','Color','k');
        end
    end
end

% Outputs
out = struct();
out.regions = cellstr(regs);
out.DiffMat = Diff;          % raw A−B TGW
out.PvalMat = P_adj;         % adjusted p
out.SigMat  = Sig;

pairs = nchoosek(1:R,2);
rows = cell(size(pairs,1),6);
for k=1:size(pairs,1)
    i=pairs(k,1); j=pairs(k,2);
    rows(k,:) = {regs(i), regs(j), Diff(i,j), P(i,j), P_adj(i,j), Sig(i,j)};
end
out.p_table = cell2table(rows, 'VariableNames', ...
    {'RegionA','RegionB','Diff_AminusB','p_two_raw','p_adj','is_sig'});

end


function out = permtest_unpaired_mean_diff(mA, mB, varargin)
% Unpaired permutation test on mean(mA) - mean(mB).
% Returns two-sided and one-sided p-values, plus basic summaries.

p = inputParser;
addParameter(p,'NPerm',10000,@(x)isnumeric(x)&&isscalar(x)&&x>=100);
addParameter(p,'Seed',[],@(x)isempty(x)||isscalar(x));
addParameter(p,'NBoot',5000,@(x)isnumeric(x)&&isscalar(x)&&x>=500);
parse(p,varargin{:});
nPerm = p.Results.NPerm;
seed  = p.Results.Seed;
nBoot = p.Results.NBoot;

if ~isempty(seed), rng(seed); end

x = mA(:); x = x(isfinite(x));
y = mB(:); y = y(isfinite(y));
nA = numel(x); nB = numel(y);
if nA<2 || nB<2, error('Need >=2 iterations per group.'); end

obs = mean(x) - mean(y);
X = [x; y];
null = zeros(nPerm,1);
for i = 1:nPerm
    idx = randperm(nA+nB);
    XA  = X(idx(1:nA));
    XB  = X(idx(nA+1:end));
    null(i) = mean(XA) - mean(XB);
end

p_two   = (sum(abs(null) >= abs(obs)) + 1) / (nPerm + 1);
p_right = (sum(null >= obs) + 1) / (nPerm + 1);
p_left  = (sum(null <= obs) + 1) / (nPerm + 1);

% Bootstrap 95% CIs
ciA = boot_ci_mean(x, nBoot);
ciB = boot_ci_mean(y, nBoot);
ciD = boot_ci_diff_unpaired(x, y, nBoot);

out = struct();
out.nA = nA; out.nB = nB;
out.meanA = mean(x); out.meanB = mean(y);
out.diff  = obs;
out.p_two = p_two; out.p_right = p_right; out.p_left = p_left;
out.null  = null;
out.ciA   = ciA; out.ciB = ciB; out.ciDiff = ciD;
end

function ci = boot_ci_mean(v, B)
N = numel(v); boots = zeros(B,1);
for b=1:B
    idx = randi(N,[N,1]);
    boots(b) = mean(v(idx));
end
ci = prctile(boots,[2.5 97.5]);
end

function ci = boot_ci_diff_unpaired(x, y, B)
nA = numel(x); nB = numel(y);
boots = zeros(B,1);
for b=1:B
    xa = x(randi(nA,[nA,1]));
    yb = y(randi(nB,[nB,1]));
    boots(b) = mean(xa) - mean(yb);
end
ci = prctile(boots,[2.5 97.5]);
end

function overlay_pairwise_brackets(ax, means, sems, P, alpha, maxPairs)
% Draw brackets with stars above bars for all pairs with p<alpha (upper-tri).
R = numel(means);
[rows, cols] = find(triu(P < alpha, 1));
pairs = [rows cols];
if isempty(pairs), return; end
if isfinite(maxPairs) && size(pairs,1) > maxPairs
    % keep the most significant ones
    ps = arrayfun(@(k) P(pairs(k,1), pairs(k,2)), 1:size(pairs,1))';
    [~,ord] = sort(ps,'ascend');
    pairs = pairs(ord(1:maxPairs),:);
end

yl = ylim(ax); yspan = range(yl);
baseHeights = max(means(:) + sems(:), [], 'omitnan');
% per bar current "level" to avoid overlaps
lev = zeros(1,R);
for k = 1:size(pairs,1)
    i = pairs(k,1); j = pairs(k,2);
    pval = P(i,j);
    stars = p_to_stars(pval);

    % determine bracket level for this pair
    L = max(lev(i), lev(j)) + 1;
    lev(i) = L; lev(j) = L;

    % y-position: top of the taller bar + margin per level
    yTop = max(means([i j]) + sems([i j]));
    y = yTop + (0.05*yspan) * L;       % 5% of y-span per level
    h = 0.02*yspan;                     % tick height

    draw_bracket(ax, i, j, y, h);
    text(mean([i j]), y + 0.01*yspan, stars, ...
        'HorizontalAlignment','center', 'VerticalAlignment','bottom', ...
        'FontWeight','bold', 'Color','k', 'Parent', ax);
end
% expand y-lims slightly to ensure brackets are visible
ylim(ax, [yl(1), max(yl(2), max(means+sems) + 0.1*yspan*(max(lev)+1))]);
end

function draw_bracket(ax, x1, x2, y, h)
% Simple bracket: vertical ticks at x1/x2 and horizontal line at y
plot(ax, [x1 x1 x2 x2], [y-h y y y-h], 'k-', 'LineWidth', 1.5);
end

function stars = p_to_stars(p)
if p < 1e-3
    stars = '***';
elseif p < 1e-2
    stars = '**';
elseif p < 0.05
    stars = '*';
else
    stars = '';
end
end

function padj = fdr_bh(p)
% Benjamini-Hochberg FDR correction (vector p -> vector padj)
p = p(:); n = numel(p);
[ps, ix] = sort(p);
rank = (1:n)';
padj_s = ps .* n ./ rank;
% enforce monotonicity
for i = n-1:-1:1
    padj_s(i) = min(padj_s(i), padj_s(i+1));
end
padj = zeros(n,1); padj(ix) = min(padj_s, 1);
end

function padj = adjust_p(p, how)
switch lower(string(how))
    case "fdr"
        padj = fdr_bh(p);
    case "bonferroni"
        padj = min(p * numel(p), 1);
    otherwise
        padj = p;
end
end


