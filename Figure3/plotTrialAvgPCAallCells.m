
%%%%%%%% PLOTTING Figure 3C, 3D 
%%%%%%%%  Facial gestures distinct, non-overlapping trajectories in neural space
%%%%%%%%
%%%%%%%%  Plots Time-resolved neural population activity trajectories (trial-averaged and smoothed, -1000 ms to +1000 ms, faint to full saturation) 
%%%%%%%     in the space of the first three principle components (all  cells across all regions).
%%%%%%%%  Each trajectory is the projection of trial-averaged neural activity during one gesture type onto the top three PCs
%%%%%%%%  
%%%%%%%% Also plots projection of trajectories from C onto first three principal components 
%%%%%%%  A and C: data from one typical day
%%%%%%%% Note that Figures 3C, 3D in the main manuscript are all cells, all regions together 
%%%%%%%% Input is data.AllDataOut.(region).score, which is nTimepoints x nPCs
%%%%%%%% 
%%%%%%%% This will also plot each region independently, as in Fig S9
% written GRI 250911

clear all; close all;
set(0,'defaultAxesFontSize',36)
set(0,'defaultAxesFontWeight','bold')

%% ============================
%% load data 
%% ============================

fileList = dir(fullfile('matfiles/Fig3C*.mat')); % change depending on which trajectory to plot 
dataFiles = {fileList.name};

for dd = 1:length(dataFiles)

    data = load(['MatFiles/' dataFiles{dd}]);
    if dd == 1
        regionLabel = 'All';
    else
        regionLabel = extractBetween(dataFiles{dd}, 'Fig3C_', '.mat') ; regionLabel = regionLabel{1};
    end
    bhvs2plot = data.bhvs2plot;
    taxis2take = data.taxis2take;
    
    latent = data.allDataOut.(regionLabel).latent;
    score = data.allDataOut.(regionLabel).score;
    maxDim = data.allDataOut.(regionLabel).nDimGreater90;

    %% visualize each PC, for each bhv
    f2 = figure; set(f2, 'Position',[476 92 1022 774]);
    sgtitle(['Region : ' regionLabel],'FontSize',36)
    colors = 'rbg';
    for PC = 1:3
        start = 1;
        for i = 1:length(bhvs2plot)
            stop = length(taxis2take)*i;
            PrincipleAxis = score(start:stop,PC);
            subplot(1,3,PC);
            plot(taxis2take,PrincipleAxis, colors(i), 'LineWidth', 6); xlabel('Time (sec)'); ylabel('a.u.'); title(['PC' num2str(PC)'])
            hold on
            start = start + length(taxis2take);
        end
    end

    %% plot perTimePt data projected onto PCs 1-3; labeled by bhvtype
    f1 = figure; set(f1, 'Position',[476 92 1022 774]);
    markerType = 'osd*.';
    markerSize = 300;
    [customColormap] = continuousColorMap(length(taxis2take));
    
    start = 1;
    for i = 1:length(bhvs2plot)
        stop = length(taxis2take)*i;
        scatter3(score(start:stop,1),score(start:stop,2),score(start:stop,3),markerSize, colormap(customColormap(:,:,i)),'filled',markerType(i))
        hold on
        start = start + length(taxis2take);
    end
    xlabel('Neural PC1')
    ylabel('Neural PC2')
    zlabel('Neural PC3')
    title(['Region: ' regionLabel ' Facial Gesture Neural Trajectories']);
    legend('Threat','Lipsmack','Chew')
    grid on; axes = gca; axes.LineWidth = 3; axes.GridLineStyle = ':';
%%     exportgraphics(gcf, ['~/Desktop/NeuralTrajs_' regions{rr}.label '.pdf'], 'ContentType', 'vector');
%     savefig(gcf, ['~/Desktop/NeuralTrajs_' regions{rr}.label '.fig']);

end

function plot_distance_shuffle_results(shuffleStats, shuffleNull, regionLabel, tAxis, varargin)
% Plot observed trajectory distances vs. shuffle null with significance shading.
%
% Required:
%   shuffleStats : output struct from test_distance_shuffle (stats)
%   shuffleNull  : output struct from test_distance_shuffle (nullStore)
%   regionLabel  : e.g., 'M1'
%   tAxis        : [T x 1] time vector corresponding to rows of distances
%
% Optional name/value:
%   'Pairs'        : cellstr of pairs to plot (default {'ThrVCh','ThrVLS','LSVCh'})
%   'ShowAverage'  : logical, add 4th subplot for 'Average' (default false)
%   'CIPrct'       : two-element vector for CI (default [2.5 97.5])
%   'ObsLineWidth' : scalar (default 2)
%   'NullLineWidth': scalar (default 1)
%
% Example:
%   plot_distance_shuffle_results(shuffleStats, shuffleNull, 'M1', taxis2take, ...
%       'ShowAverage', true, 'CIPrct', [5 95]);

p = inputParser;
addParameter(p,'Pairs', {'ThrVCh','ThrVLS','LSVCh'});
addParameter(p,'ShowAverage', false);
addParameter(p,'CIPrct', [2.5 97.5]);
addParameter(p,'ObsLineWidth', 2);
addParameter(p,'NullLineWidth', 1);
parse(p, varargin{:});
Pairs        = p.Results.Pairs;
ShowAverage  = p.Results.ShowAverage;
CIPrct       = p.Results.CIPrct;
lwObs        = p.Results.ObsLineWidth;
lwNull       = p.Results.NullLineWidth;

if strcmp(regionLabel, 'All') || strcmp(regionLabel, 'AllMotor')
    thisColor =   [0.25 0.25 0.25];
elseif strcmp(regionLabel, 'S1')
    thisColor = [0.4940 0.1840 0.5560];
elseif strcmp(regionLabel, 'M1')
    thisColor = [0.6350 0.0780 0.1840];
elseif strcmp(regionLabel,'PMv')
    thisColor = [0.8500 0.3250 0.0980];
elseif strcmp(regionLabel, 'M3')
    thisColor = [0.9290 0.6940 0.1250];
end

if ShowAverage && ~ismember('Average', Pairs)
    Pairs = [Pairs, {'Average'}];
end

nPlots = numel(Pairs);
T = numel(tAxis);

% Basic figure layout
tiledlayout(1, nPlots, 'Padding', 'compact', 'TileSpacing', 'compact');

for i = 1:nPlots
    pairName = Pairs{i};

    obs   = shuffleStats.(regionLabel).(pairName).obs(:);
    p_fdr = shuffleStats.(regionLabel).(pairName).p_fdr(:);
    sig   = shuffleStats.(regionLabel).(pairName).sig_fdr(:);
    nullM = shuffleNull.(regionLabel).(pairName);  % [T x nShuf]

    % Null summary
    muNull = mean(nullM, 2, 'omitnan');
    ciNull = prctile(nullM, CIPrct, 2);  % [T x 2]

    nexttile; hold on;

    % 1) Null CI band
    fill_between(tAxis, ciNull(:,1), ciNull(:,2), [0.85 0.85 0.85]); % light gray band

    % 2) Null mean (dashed)
    hNull = plot(tAxis, muNull, '--', 'LineWidth', lwNull, 'Color', [0.3 0.3 0.3]);

    % 3) Observed curve
    hObs = plot(tAxis, obs, 'LineWidth', lwObs,'Color',thisColor);

    % 4) Shade significant bins in the background (light tint of obs color)
    yl = get(gca,'YLim');
    shade_sig_runs(tAxis, sig, yl, hObs.Color);

    % Cosmetics
    grid on;
    xlabel('Time');
    ylabel('Euclidean Distance');
    title([regionLabel '  ' pairName], 'Interpreter','none');
    uistack(hObs,'top');  % keep obs on top of band/patches

    if i == nPlots
        legend([hObs hNull], {'Observed','Null'}, 'Location','best'); % keep legend simple
    end
end

end

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
