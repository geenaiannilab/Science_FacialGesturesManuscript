%%%%%%%
%%%%%%% Plot Figure S2
%%%%%%  Mutual information distributions by region 
%%%%%%  Imported data includes MI_corr_all which is bias corrected mutual
%%%%%%      information (see Methods), one value per region, and regions_all
%%%%%%      which is cortical region label 
%%%%%%  MI Computed between each neuron's spike counts and gesture (trial type) using empirical joint probability distributions.
%%%%%%  Data reflect the analysis window tmin = -500 ms to tmax = +500 ms.
%%%%%%
%%%%%   Also includes fracSig, which is fraction of neurons per region
%%%%%       whose p-value (for MIc being significant > null) was < 0.05
%%%% written GRI 


clear all ; close all ; 
set(0,'defaultAxesFontSize',20)
plotExtended = false; 

load('matfiles/FigS2.mat'); 
colorMap = [0.494 0.184 0.556; 0.635 0.078 0.184; 0.85 0.325 0.098; 0.929 0.694 0.125];

%% --- POOLED: Bias-corrected MI violin with overlays ---

% --- Bias-corrected MI by region ---
MI_corrPlot_all = '';
for rr = 1:length(uRegs) 
    thisRegion = uRegs{rr};
    theseCells = strcmp(thisRegion,regions);
    MI_corrPlot_all{rr} = MI_corr_all(theseCells)';
end

%% --- Plot Bias-corrected MI by region---
figure; 
violinplot(MI_corrPlot_all, regions_all,'facecolor', colorMap);          % returns one object per region
xticks(1:nRegs); xticklabels(uRegs);l = legend; 
ylabel('Bias-corrected MI (bits)');
title(['Per-region MI_{corr} (pooled across days)', 'tmin: ' num2str(time2keep(1)) ' , tmax:' num2str(time2keep(2))]);
axes = gca; axes.FontSize = 28; axes.FontWeight = 'bold';

[p_corr,~,stats_corr] = kruskalwallis(MI_corr_all, regions_all, 'off');
if p_corr < 0.05
    results = multcompare(stats_corr,'CType','dunn-sidak','Display','off');
    overlaySig(results);
end
l.String = l.String(1:2);

%% --- POOLED: fraction significant per region ---
figure('Name','Fraction significant (pooled)','NumberTitle','off');
b = bar(fracSig,'FaceColor','flat');
b.CData = colorMap;
xticklabels(uRegs);
ylabel('Fraction of significant neurons (p<0.05)');
title(['Significant neuron fractions per region (pooled) tmin: ' num2str(time2keep(1)) ' , tmax:' num2str(time2keep(2))]);

% --- Chi-square test of independence on pooled fractions ---
fprintf('Chi-square (pooled) across regions: χ²(%d) = %.2f, p = %.4f\n', ...
        stats_chi.df, chi2stat, p_chi);

%% --- EXTENDED PLOTTING not in paper 
if plotExtended

    %% --- POOLED: Raw MI violin with overlays ---
    
    % --- Raw MI by region ---
    MI_rawPlot_all = '';
    for rr = 1:length(uRegs) 
        thisRegion = uRegs{rr};
        theseCells = strcmp(thisRegion,regions);
        MI_rawPlot_all{rr} = MI_raw_all(theseCells)';
    end
    
    % --- Plot Raw MI by region---
    figure; 
    violinplot(MI_rawPlot_all, regions_all,'facecolor', colorMap);          % returns one object per region
    xticks(1:nRegs); xticklabels(uRegs);
    ylabel('Raw MI (bits)');
    title(['Per-region MI_raw (pooled across days), tmin: ' num2str(time2keep(1)) ' , tmax:' num2str(time2keep(2))]);
    
    [p_raw,~,stats_raw] = kruskalwallis(MI_raw_all, regions, 'off');
    if p_raw < 0.05
        results = multcompare(stats_raw,'CType','dunn-sidak','Display','off');
        overlaySig(results);
    end
    
    %% --- PER-DAY VIEW: fraction significant per region per day ---
    % make a matrix Day x Region of fractions, then heatmap + dot plot
    [uDays,~,idDay] = unique(PerDay.Day,'stable');
    [uRegsDay,~,idReg] = unique(PerDay.Region,'stable');
    F = accumarray([idDay,idReg], PerDay.FracSig, [], @mean, NaN);
    
    figure('Name','Fraction significant per day × region','NumberTitle','off');
    h = heatmap(uRegsDay, uDays, F);
    h.XLabel = 'Region'; h.YLabel = 'Day'; h.Title = 'Frac. significant (p<0.05)';
    colormap(parula); colorbar;
    
end

%% --- Helper function to overlay sig bars ---
function overlaySig(results)
    ylim_curr = ylim;
    y_top = ylim_curr(2);
    step = 0.05 * range(ylim_curr); % vertical spacing
    for r = 1:size(results,1)
        g1 = results(r,1); g2 = results(r,2); pval = results(r,6);
        if pval < 0.05
            addSigBar(g1,g2,y_top + r*step,pval);
        end
    end
    ylim([ylim_curr(1) y_top + (size(results,1)+2)*step]);
end

%% --- Draw a bar with stars ---
function addSigBar(x1,x2,y,p)
    line([x1 x2],[y y],'Color','k','LineWidth',1.5);
    line([x1 x1],[y-0.01 y],'Color','k','LineWidth',1.5);
    line([x2 x2],[y-0.01 y],'Color','k','LineWidth',1.5);
    if p < 0.001
        stars = '***';
    elseif p < 0.01
        stars = '**';
    else
        stars = '*';
    end
    text(mean([x1 x2]), y+0.01, stars, ...
        'HorizontalAlignment','center','FontSize',12);
end

function [chi2stat, p, stats] = chi2indep(counts)
    % counts: contingency table (rows = groups, cols = categories)
    % returns chi-square statistic, p-value, and stats struct

    rowSums = sum(counts,2);
    colSums = sum(counts,1);
    N       = sum(rowSums);

    expected = (rowSums * colSums) / N;
    chi2stat = sum((counts - expected).^2 ./ expected,'all');

    df = (size(counts,1)-1) * (size(counts,2)-1);
    p  = 1 - chi2cdf(chi2stat, df);

    stats.df = df;
    stats.expected = expected;
end
