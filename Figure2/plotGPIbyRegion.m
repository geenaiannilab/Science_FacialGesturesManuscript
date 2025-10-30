
%%%%%%%%
%%%%%%%% PLOTTING Figure 2D 
%%%%%%%%  Single-Cell Activity and Selectivity in Cortical Face-Motor Regions
%%%%%%%%
%%%%%%%% GPI index by region with overlying statistics (no effect of cortical region on GPI)  
%%%%%%%% GPI index by region, only for selective cells (GPI > 0.5) 
%%%%%%%% Input data includes 'combinedResults', which contains
%%%%%%%%    'faceExpPrefIndex' (GPI) and total number of cells per region 
%%%%%%%%% GRI 09/11/2025 

clear all
set(0,'defaultAxesFontSize',20)

%% Control options
regions = {'S1','M1','PMv','M3'};
subjects = {'barney','thor'};
GPIthreshold = 0.5;

%% Load matfiles
pathtoData = 'matfiles/Fig2D_GPIvalues.mat';
load(pathtoData);

%% get the cortical regions per cell 
cortRegionOut = [];
GPIout = [];
meanFRsOut = [];

for ss = 1:length(subjects)
    if strcmpi(subjects{ss},'Barney')
        data = barneyRawData;
    elseif strcmpi(subjects{ss},'Thor')
        data = thorRawData;
    end

    % pool faceExpPrefs indices
    GPI = data.faceExpPref';
    cortRegion = data.cortRegion;
    meanFRs = data.meanFRs;

    cortRegionOut = [cortRegionOut; cortRegion];
    GPIout = [GPIout ; GPI];
    meanFRsOut = [meanFRsOut; meanFRs];

    clear cortRegion; clear GPI; clear meanFRs;
end

%% run 1 way anova for effect of cortical region on GPI 
% and run post-hoc tests 
[P_oneway,ANOVATAB_oneway,STATS_oneway] = anova1(GPIout, cortRegionOut,'off');

%% report it 
F = ANOVATAB_oneway{2,5};
df_between = ANOVATAB_oneway{2,3};
df_within = ANOVATAB_oneway{3,3};
p = ANOVATAB_oneway{2,6};
disp('Effect of cortical region on GPI (One way ANOVA):')
disp(['F(' num2str(df_between) ',' num2str(df_within) ') = ' num2str(F) ', p = ' num2str(p)]);

%% run 2 way anova for effect of cortical region and gesture type on meanFRs 
inputFR = [meanFRsOut(:,1) ; meanFRsOut(:,2); meanFRsOut(:,3)];
T    = cell(1, size(meanFRsOut,1)); T(:) = {'Threat'}; L = cell(1, size(meanFRsOut,1)); L(:) = {'Lipsmack'}; C = cell(1, size(meanFRsOut,1)); C(:) = {'Chew'};
inputExp = [T'; L'; C'];
inputCortRegion = [cortRegionOut ; cortRegionOut; cortRegionOut];
[P_twoway,TABLE_twoway,STATS_twoway] = anovan(inputFR,{inputExp inputCortRegion},'model','interaction','varnames',{'inputExp','inputCortRegion'},'display','off');
%[c_twoway,m_twoway, ~, names] = multcompare(STATS_twoway,'Dimension', [1 2]);

%% report it 
get_row = @(name) find(strcmp(TABLE_twoway(2:end,1), name)) + 1;  % +1 to offset header

row_gesture     = get_row('inputExp');
row_region      = get_row('inputCortRegion');
row_interaction = get_row('inputExp*inputCortRegion');
row_error       = size(TABLE_twoway,1)-1;   
df_err = TABLE_twoway{row_error,3};

F_g   = TABLE_twoway{row_gesture,6};     df_g   = TABLE_twoway{row_gesture,3};     p_g   = TABLE_twoway{row_gesture,7};
F_r   = TABLE_twoway{row_region,6};      df_r   = TABLE_twoway{row_region,3};      p_r   = TABLE_twoway{row_region,7};
F_int = TABLE_twoway{row_interaction,6}; df_int = TABLE_twoway{row_interaction,3}; p_int = TABLE_twoway{row_interaction,7};

disp('Effect of gesture (Two-way ANOVA):')
disp(['F(' num2str(df_g) ',' num2str(df_err) ') = ' num2str(F_g,'%.3f') ', p = ' num2str(p_g,'%.3g')])

disp('Effect of cortical region (Two-way ANOVA):')
disp(['F(' num2str(df_r) ',' num2str(df_err) ') = ' num2str(F_r,'%.3f') ', p = ' num2str(p_r,'%.3g')])

disp('Gesture × region interaction (Two-way ANOVA):')
disp(['F(' num2str(df_int) ',' num2str(df_err) ') = ' num2str(F_int,'%.3f') ', p = ' num2str(p_int,'%.3g')])

%% just the "selective" cells (GPI > 0.5) 
selectiveGPI = GPIout >= GPIthreshold;
meanFRs = meanFRsOut(selectiveGPI,:);

xRegions = categorical(cortRegionOut(selectiveGPI));
xRegions = reordercats(xRegions,{'S1' 'M1' 'PMv','M3'});
yGPIvalues = GPIout(selectiveGPI);
[cExpWinnerFR, cExpWinner] = max(meanFRs,[],2);

%% plot GPI by region, 
cmap1 = [         
    0.4940 0.1840 0.5560	
    0.6350 0.0780 0.1840
    0.8500 0.3250 0.0980	
    0.9290 0.6940 0.1250	
         ];
cmap2 = [         0.4940 0.1840 0.5560
    0.6350 0.0780 0.1840
    0.8500 0.3250 0.0980	
    0.9290 0.6940 0.1250	
        ];
    b(2).FaceAlpha = 0.5;

%% plots subjects' GPI by region, separately    

xall = [];
 yall = [];
 colorAll = [];
 faceAlphaAll = [];
 counter = 1;
for ii = 1:length(regions)
    for subj = 1:2
        if subj == 1
            thisRegion = barneyResults.(regions{ii});
            color = cmap1(ii,:);
            faceAlpha = 0.8;
            %color = [0 0.4470 0.7410];
        elseif subj == 2
            thisRegion = thorResults.(regions{ii});
            color = cmap1(ii,:);
            faceAlpha = 0.5;
            %color = [0.8500 0.3250 0.0980];
        end
        x = ones(thisRegion.totalNumberOfCells,1)*(counter);
        y = thisRegion.faceExpPrefIndex';

        xall = cat(1,xall, x);
        yall = cat(1,yall, y);
        colorAll = cat(1,colorAll, color);
        faceAlphaAll = cat(1,faceAlphaAll, faceAlpha);
        if subj == 1
            counter = counter +1;
        elseif subj == 2
            counter = counter +1.75;
        end
    end
end

plotFig = figure;
b = beeswarm(xall,yall,'Colormap', colorAll,'dot_size',1.2,'MarkerFaceAlpha',faceAlphaAll','overlay_style','ci','sort_style','fan');
xticks([1.5 4.25 7 9.75 ]);
xticklabels({'S1', 'M1','PMv', 'M3'})
xlabel('Region')
ylim([0 1.1]);
ylabel('Preference Index');
title('Gesture Preference Indices per Regions')
legend('Monkey B','','','Monkey T');

%% concatente both animals for plotting 
cmap1 = [         
    0.4940 0.1840 0.5560	
    0.6350 0.0780 0.1840
    0.8500 0.3250 0.0980	
    0.9290 0.6940 0.1250	];

 xall = [];
 yall = [];
 colorAll = [];
 faceAlphaAll = [];
 counter = 1;

 for ii = 1:length(regions)

     thisRegion = combinedResults.(regions{ii});
     color = cmap1(ii,:);
     faceAlpha = 0.8;

     x = ones(thisRegion.totalNumberOfCells,1)*(counter);
     y = thisRegion.faceExpPrefIndex';

     xall = cat(1,xall, x);
     yall = cat(1,yall, y);
     colorAll = cat(1,colorAll, color);
     faceAlphaAll = cat(1,faceAlphaAll, faceAlpha);
     counter = counter +2;

 end

 %% plot subjects together 
plotFig = figure;
b = beeswarm(xall,yall,'Colormap', colorAll,'dot_size',1.2,'MarkerFaceAlpha',faceAlphaAll','overlay_style','ci','sort_style','fan');
xticks([1 3 5 7]);
xticklabels({'S1', 'M1','PMv', 'M3'})
xlabel('Region')
ylim([0 1.1]);
ylabel('Preference Index');
title('Gesture Preference Indices per Regions')

%% selective cells 
figure;
swarmchart(xRegions,yGPIvalues, 200, cExpWinner,'filled','MarkerFaceAlpha',0.75,'MarkerEdgeAlpha',0.5);
expColorMap = [1 0 0
    0 1 0
    0 0 1]; % r,b,g
colormap(expColorMap)
yticks(0.5:0.1:1)
ylim([.45 1.05])
ylabel('Preference Index')
xlabel('Region')
title('Preference Indices of Selective Cells')
axes = gca;
axes.FontWeight = 'bold';

%% === APA-style reporting for two-way ANOVA + Tukey post-hoc ===

function print_anova_APA(T, varargin)
% PRINT_ANOVA_APA  Pretty-print one- or two-way ANOVA in APA style.
% Works with anovan's cell-array table or MATLAB's table.
%
% Usage:
%   print_anova_APA(TABLE);
%   print_anova_APA(TABLE, );
%   print_anova_APA(TABLE, 'IncludeEffectSize', true);  % adds partial eta^2
%
% Options:
%   'FactorLabels'      : struct or containers.Map mapping exact row names in T
%                         to friendly labels. If omitted, uses row names as-is.
%   'IncludeEffectSize' : logical (default: true). Reports ηp^2 = SS_effect/(SS_effect+SS_error).
%   'DigitsF'           : integer (default: 2) decimals for F.
%   'DigitsEta'         : integer (default: 2) decimals for ηp^2.

p = inputParser;
p.addParameter('FactorLabels', struct(), @(x)isstruct(x) || isa(x,'containers.Map'));
p.addParameter('IncludeEffectSize', true, @(x)islogical(x) && isscalar(x));
p.addParameter('DigitsF', 2, @(x)isnumeric(x)&&isscalar(x));
p.addParameter('DigitsEta', 2, @(x)isnumeric(x)&&isscalar(x));
p.parse(varargin{:});
opts = p.Results;

% --- Normalize table-like input to cell matrix & headers ---
[hdr, rows, C] = normalizeTable(T); % hdr: string row of headers; rows: string row names; C: cell data

% --- Locate key columns (robust to header variations) ---
colDF = findFirst(hdr, ["df","d.f","DF","D.F"], true);
colF  = findExact(hdr, "F");
colP  = findFirst(hdr, ["Prob>F","p","P","Pr(>F)"], true);
colSS = findFirst(hdr, ["SS","Sum Sq","SumSq","Sum of Squares","SumSq.","Sum"], true);

% --- Find error/residual row ---
rErr = find(contains(lower(rows), ["error","residual"]), 1, 'first');
if isempty(rErr)
    warning('Could not find Error/Residual row; effect sizes will be skipped.');
end

% --- Identify effect rows (exclude Total and Error) ---
isTotal = contains(lower(rows), "total");
effectIdx = find(~isTotal & (1:numel(rows))' ~= rErr);

% Partition into main effects vs interaction (if two-way)
isInteraction = contains(rows, "*") | contains(rows, "×") | contains(rows, "x");
rInt = effectIdx(isInteraction(effectIdx));
rMain = effectIdx(~isInteraction(effectIdx));

% Friendly label resolver
labelOf = @(rname) resolveLabel(rname, opts.FactorLabels);

% Header
fprintf('ANOVA (APA-style):\n');

% --- Print main effects ---
for r = rMain(:)'
    printOne(rows(r), r, C, colDF, colF, colP, colSS, rErr, labelOf, opts);
end

% --- Print interaction if present ---
for r = rInt(:)'
    % Replace * with × for APA typography
    lbl = strrep(labelOf(rows(r)), '*', ' × ');
    printOneWithCustomLabel(lbl, r, C, colDF, colF, colP, colSS, rErr, opts);
end

end % main

% ================= Helpers =================

function printOne(rname, r, C, colDF, colF, colP, colSS, rErr, labelOf, opts)
    lbl = labelOf(rname);
    printOneWithCustomLabel(lbl, r, C, colDF, colF, colP, colSS, rErr, opts);
end

function printOneWithCustomLabel(lbl, r, C, colDF, colF, colP, colSS, rErr, opts)
    [df1, Fv, pval] = grabStats(C, r, colDF, colF, colP);
    if isnan(df1) || isnan(Fv) || isnan(pval)
        return
    end
    df2 = grabDF2(C, rErr, colDF);

    % Build base string: F(df1, df2) = x.xx, p = .xxx
    Ffmt = sprintf(['%%.%df'], opts.DigitsF);
    base = sprintf('%s: F(%d, %d) = %s, %s', lbl, df1, df2, sprintf(Ffmt, Fv), formatP(pval));

    % Effect size (partial eta^2)
    es = '';
    if opts.IncludeEffectSize && ~isempty(colSS) && ~isnan(df2) && ~isempty(rErr)
        SSe = toNum(C{r, colSS});
        SSerr = toNum(C{rErr, colSS});
        if isfinite(SSe) && isfinite(SSerr) && (SSe+SSerr) > 0
            etap2 = SSe / (SSe + SSerr);
            Etafmt = sprintf(['%%.%df'], opts.DigitsEta);
            es = sprintf(', \x03B7p^2 = %s', sprintf(Etafmt, etap2)); % ηp^2
        end
    end

    fprintf('%s%s\n', base, es);
end

function [hdr, rows, C] = normalizeTable(T)
    if istable(T)
        % Convert to cell, keep VariableNames as headers, RowNames (if any) as first col
        hdr = string(['Row', T.Properties.VariableNames]);
        if ~isempty(T.Properties.RowNames)
            rowNames = string(T.Properties.RowNames(:));
        else
            % assume first column holds row labels if present in data
            % otherwise synthesize row labels
            if iscellstr(T{:,1}) || isstring(T{:,1})
                rowNames = string(T{:,1});
                T = T(:, 2:end);
                hdr = string(['Row', T.Properties.VariableNames]);
            else
                rowNames = string(compose("Row%d", (1:height(T))'));
            end
        end
        C = T{:,:};
        % Ensure everything is cell for uniform indexing
        C = num2cell(C);
        rows = rowNames;
        % Prepend fake first column (row names) into C for consistent getVal (not used further)
    else
        % Assume cell array like anovan output: first row = header, first col = row labels
        C = T;
        hdr = string(C(1,:));
        rows = string(C(2:end,1));
        C = C; 
    end

    % If T was a MATLAB table, rebuild a cell with header row to match downstream expectations
    if istable(T)
        % Compose a cell array with header row compatible with original user's code
        headerRow = cell(1, numel(hdr));
        for i = 1:numel(hdr), headerRow{i} = char(hdr(i)); end
        data = T{:,:};
        % Promote numeric to cell
        datac = num2cell(data);
        C = [headerRow; [cellstr(rows) datac]];
        hdr = string(C(1,:));
        rows = string(C(2:end,1));
    end
end

function idx = findFirst(hdr, keys, ignoreCase)
    if nargin < 3, ignoreCase = true; end
    idx = [];
    for k = 1:numel(keys)
        if ignoreCase, hit = find(contains(hdr, keys(k), 'IgnoreCase', true), 1);
        else,          hit = find(contains(hdr, keys(k)), 1);
        end
        if ~isempty(hit), idx = hit; return; end
    end
end

function idx = findExact(hdr, key)
    idx = find(strcmpi(hdr, key), 1);
end

function out = toNum(x)
    if isnumeric(x), out = x; return; end
    if isstring(x) || ischar(x), out = str2double(x); else, out = NaN; end
end

function [df1, Fv, p] = grabStats(C, r, colDF, colF, colP)
    df1 = NaN; Fv = NaN; p = NaN;
    try
        df1 = toNum(C{r+1, colDF});  % +1 because C includes header in row 1
        Fv  = toNum(C{r+1, colF});
        p   = toNum(C{r+1, colP});
    catch, end
end

function df2 = grabDF2(C, rErr, colDF)
    if isempty(rErr) || isempty(colDF), df2 = NaN; return; end
    try
        df2 = toNum(C{rErr+1, colDF});
    catch
        df2 = NaN;
    end
end

function s = resolveLabel(rname, labelMap)
    s = char(rname);
    if isstruct(labelMap)
        f = fieldnames(labelMap);
        for i = 1:numel(f)
            if strcmp(rname, f{i}), s = labelMap.(f{i}); return; end
        end
    elseif isa(labelMap, 'containers.Map')
        if isKey(labelMap, char(rname)), s = labelMap(char(rname)); end
    end
    % Beautify a bit
    s = strrep(s, '_', ' ');
end

function s = formatP(p)
% APA p-value formatting
    if ~isfinite(p) || isnan(p)
        s = 'p = NaN';
        return
    end
    if p < 0.001
        s = 'p < .001';
    else
        s = ['p = ' strip(erase(sprintf('%.3f', p), '0'), 'right')];
        % Ensure leading zero omitted per APA: .012 not 0.012
        if startsWith(s, 'p = 0.'), s = strrep(s, 'p = 0.', 'p = .'); end
        % Ensure exactly three decimals if needed
        if ~contains(s, '<') && ~contains(s, '.')
            s = sprintf('p = %.3f', p);
            s = strrep(s, 'p = 0.', 'p = .');
        end
    end
end
