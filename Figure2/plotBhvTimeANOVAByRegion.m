%%%%%%%%
%%%%%%%% PLOTTING Figure 2E,  
%%%%%%%%  Single-Cell Activity and Selectivity in Cortical Face-Motor Regions
%%%%%%%%
%%%%%%%% Plots results of 2-way anova (effect of gesture type, time and
%%%%%%%% their interaction on firing rate of individual neurons), over all
%%%%%%%%  regions. Then use chi-sq to test if there exists differences in
%%%%%%%%  fractions of significantly modulated cells by region
%%%%%%%%  Input data are  statistical table; output of chiSquaredGOF.m
%%%%%%%%  (FDR-corrected)

%%%%%%%%% GRI 09/11/2025 


clear all; close all 
set(0,'defaultAxesFontSize',24)

%% Load matfiles
S = load('matfiles/Fig2E.mat');   
Stats = S.Stats;

%% --- Config ---
regionOrder   = ["S1","M1","PMv","M3"];          % x-axis order
factorsToPlot = ["Behavior","Time","Interact"];  % must match Stats.Factor values (case-insensitive)
subtitleFactors = {'Gesture','Time','Gesture x Time'};

% Fixed colormap for regions: [S1; M1; PMv; M3]
barColors = [0.4940 0.1840 0.5560   % purple
             0.6350 0.0780 0.1840   % red
             0.8500 0.3250 0.0980   % orange
             0.9290 0.6940 0.1250]; % yellow

%% --- Figure & layout ---
figure('Color','w'); 
tl = tiledlayout(1,3,'TileSpacing','compact','Padding','compact');

for ff = 1:numel(factorsToPlot)
    ax = nexttile; hold(ax,'on');

    % ---------- Pull OVERALL row for this factor ----------
    isOverall = strcmpi(Stats.Level,'Overall') & strcmpi(Stats.Factor, factorsToPlot(ff));
    row = Stats(isOverall, :);

    vals = nan(1,numel(regionOrder));
    if ~isempty(row)
        for r = 1:numel(regionOrder)
            vn = "PropSig_" + regionOrder(r);
            if ismember(vn, string(row.Properties.VariableNames))
                vals(r) = double(row.(vn));
            end
        end
        % Convert to % if stored as fractions
        if nanmax(vals) <= 1, vals = 100*vals; end
    end


    % Draw bars with your fixed colors
    b = bar(ax, vals, 'FaceColor','flat');
    for r = 1:numel(regionOrder)
        b.CData(r,:) = barColors(r,:);
    end
    xticks(1:numel(regionOrder)); xticklabels(regionOrder);
    ylabel('% Significant Cells'); title(string(factorsToPlot(ff)));

    % Fix y-limits for all subplots
    ylim([0 80]); box off;

    % ---------- Significance overlays (close to bars) ----------
    placeTightSignificance(ax, Stats, factorsToPlot(ff), regionOrder, vals);


end

%% ===================== Helpers =====================

function placeTightSignificance(ax, Stats, factorName, regionOrder, vals)
    % Prefer FDR p if present, else raw
    rowO = Stats(strcmpi(Stats.Level,'Overall') & strcmpi(Stats.Factor,factorName),:);
    pG   = NaN;
    if ~isempty(rowO), pG = firstNonNan([rowO.p_adj, rowO.p_raw]); end

    TP = Stats(strcmpi(Stats.Level,'Pairwise') & strcmpi(Stats.Factor,factorName),:);

    yl   = [0 80];                 % fixed as requested
    span = yl(2)-yl(1);
    topBar = max(vals,[],'omitnan');

    % Start brackets just above tallest bar, but under the cap
    basePad   = 1.0;               % gap above tallest bar
    textPad   = 0.01;               % gap for label text
    barH      = 0.3;               % bracket vertical tab
    maxTop    = 79.0;              % never exceed 79 so text/bracket fit

    y0 = min(maxTop, topBar + basePad);

    % Compute pairwise count and spacing that fits under 80
    nPairs = height(TP);
    if nPairs==0
        % Only global
        drawSigBracket(ax, 1, numel(regionOrder), y0, p2stars(pG), ...
            'BarHeight', barH, 'TextOffset', textPad, 'LineWidth', 1.6);
        return;
    end

    % Try a comfortable spacing; shrink if needed to fit under the cap
    dy = 2.0;  % preferred spacing
    neededTop = y0 + nPairs*dy + barH + textPad;
    if neededTop > maxTop
        % squeeze spacing to fit
        dy = max(0.6, (maxTop - y0 - barH - textPad) / max(nPairs,1));
    end

    % Global bracket just above tallest bar
    drawSigBracket(ax, 1, numel(regionOrder), y0, p2stars(pG), ...
        'BarHeight', barH, 'TextOffset', textPad, 'LineWidth', 1.6);

    % Pairwise stacked tightly above
    level = 1;
    for r = 1:nPairs
        cmp = string(TP.Comparison(r));          % e.g., "S1 vs M1"
        [i,j] = parsePair(cmp, regionOrder);
        if isnan(i) || isnan(j), continue; end
        p = firstNonNan([TP.p_adj(r), TP.p_raw(r)]);

        yk = y0 + level*dy;
        if yk + barH + textPad > maxTop
            % If we still exceed, stop drawing further brackets
            break;
        end
        drawSigBracket(ax, min(i,j), max(i,j), yk, p2stars(p), ...
            'BarHeight', barH, 'TextOffset', textPad, 'LineWidth', 1.6);
        level = level + 1;
    end
end

function [a,b] = parsePair(cmp, regionOrder)
    parts = split(cmp,"vs");
    if numel(parts)~=2, a=NaN; b=NaN; return; end
    sA = strtrim(parts(1)); sB = strtrim(parts(2));
    a = find(regionOrder==sA,1,'first'); if isempty(a), a=NaN; end
    b = find(regionOrder==sB,1,'first'); if isempty(b), b=NaN; end
end

function lab = p2stars(p)
    if isnan(p), lab = 'n.s';
    elseif p < 0.001, lab = '***';
    elseif p < 0.01,  lab = '**';
    elseif p < 0.05,  lab = '*';
    else, lab = 'n.s';
    end
end

function p = firstNonNan(vec)
    p = NaN;
    for v = vec
        if ~isnan(v) && ~isempty(v)
            p = v; return;
        end
    end
end

function drawSigBracket(ax, x1, x2, y, label, varargin)
    args = struct('BarHeight',0.02,'TextOffset',0.01,'LineWidth',1.6);
    for k=1:2:numel(varargin), args.(varargin{k}) = varargin{k+1}; end
    hold(ax,'on');
    plot(ax, [x1 x1 x2 x2], [y y+args.BarHeight y+args.BarHeight y], ...
        'k','LineWidth',args.LineWidth,'Parent',ax);
    text(mean([x1 x2]), y+args.BarHeight+args.TextOffset, label, ...
        'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',14, 'FontWeight','bold','Parent',ax);
end
