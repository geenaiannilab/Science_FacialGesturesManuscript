%%%%%%%%
%%%%%%%%% PURPOSE: Analyze potential differences in latencies (timepoints
%%%%%%%%%   of max influence, of cell spiking --> movement) across all cells 
%%%%%%%%
%%%%%%%%  Will download previously generated/defined STAs *only within
%%%%%%%%  particular facial ROIs* &  calculate VARIANCE as well as MEAN
%%%%%%%%  FRAME, inside those ROIs, at each timepoint post-spike. Then plot
%%%%%%%%  it all, in order to examine potential differences across cells
%%%%%%%%
%%%%%%%%  (Does exact same thing as analyzeSTMALatencyWholeFrame.m, but
%%%%%%%%    within facial ROIs rather than across whole frame)
%%%%%%%%
%%%%%%%%  Will generate 4 plots per cortical region -- variance/mean frame
%%%%%%%%  for positive STAs, and variance/mean frame for negative STAs
%%%%%%%%   Each plot has 1 subplot per facial ROI
%%%%%%%%
%%%%%%%% OUTPUTS: none
%%%%%%%% DEPENDENCIES:  Assumes you already ran generateSTMAbyROI.m to generate 'allCellsSTMAbyROI_(null).mat'

clear all; 
close all;
set(0,'defaultAxesFontSize',20);

date = '171010';
subject ='Thor';
sess2run = 6;
collapseAreasFlag = 1;
saveFlag = 0;

%% STMA variables 
null = 'poiss';
cleanSTAflag = 1;

if strcmp(subject,'Thor')
    bhvSR = 120;
else
    bhvSR = 69.99; 
end

windowSize = 1; % seconds to examine STMA 
directionality = 'post'; %pre/post
zScoreThreshold = 3; % 
thresholdPositive = 0; 
thresholdNegative = 0; 
thresholdBoth = 1;
nFrames = ceil(bhvSR*windowSize);  % to eval post-spike
spatialFilterSize = 3; % size of 2D gauss filter for smoothing

ROIs = {'mouth','Lcheek','Rcheek','nose','eyes','Lear','Rear'};

channels2run =  [ 2 1; 9 1; 11 1; 22 1; 24 1; 32 1;...
    33 1; 34 1; 34 2; 36 1; 36 2; 38 1; 38 2; 39 1; 39 2; 40 1; 41 1; 41 2; 42 1; 42 2; 42 3; 43 1; 44 1; 48 1; 50 1 ; 52 2; 52 3; 53 1; 53 2; 54 1; 55 1; 56 1; 57 1; ...
    59 1; 59 2; 59 3; 60 1; 60 2; 61 1 ;  61 2; 62 1; 62 2; 64 1; 64 2; ...
    66 1; 66 2; 74 1; 74 2; 86 1; 86 2; 91 1; 91 2; 94 1; ...
    98 1; 100 1; 101 1; 102 1; 102 2; 103 1; 103 2; 104 1; 105 1; 105 2; 106 1; 108 1; 108 2; 109 1; 109 2; 110 1; 111 1; 112 1; 116 1; 123 1; 124 1; 125 1;...
    133 1; 137 1; 137 2; 138 1; 141 1; 141 2; 142 1; 142 2; 145 1; 145 2; 153 1; 157 1; 159 1; ...
    168 1; 179 1; 183 1; 183 2; 185 1; 186 1; 190 1; 191 1; 191 2; 191 3];
%%
% load the facial ROIs STAs
STMAdir = ['/Users/geena/Dropbox/FaceExpressionAnalysis/STMA2021/' subject '_' date];
bhvDir = (['/Users/geena/Dropbox/PhD/SUAinfo/' subject '_' date '/Data4Analysis']);

% load bhv data
filePattern = fullfile(bhvDir, ['bhvspikes_sub' num2str(sess2run) '*']);
bhvFile = dir(filePattern);
TDTsessID = cell2mat(extractBetween(bhvFile.name,['sub' num2str(sess2run) '_'],['_']));
bhvFile = bhvFile.name;
workDir = [STMAdir '/' subject '-' date '-' TDTsessID '/'];

load([workDir '/allCellsSTMAbyROI_' null '.mat']);
timeAx = 1/bhvSR:(1/bhvSR):(nFrames/bhvSR);

% for each cell & each ROI, calculate variance & mean frame, for positive &
% negative STAs
for unit = 1:length(channels2run)
    
    for ROI = 1:length(ROIs)
        disp(channels2run(unit,:))
        
        posROI= reshape((ROIsSTMA(unit).pos.(ROIs{ROI})), [size(ROIsSTMA(unit).pos.(ROIs{ROI}),1)*size(ROIsSTMA(unit).pos.(ROIs{ROI}),2), size(ROIsSTMA(unit).pos.(ROIs{ROI}),3)]);
        negROI = reshape((ROIsSTMA(unit).neg.(ROIs{ROI})), [size(ROIsSTMA(unit).neg.(ROIs{ROI}),1)*size(ROIsSTMA(unit).neg.(ROIs{ROI}),2), size(ROIsSTMA(unit).neg.(ROIs{ROI}),3)]);
        
        ROIstats.variance.pos(unit,ROI,:) = var(posROI);
        ROIstats.variance.neg(unit,ROI,:) = var(negROI);
        
        ROIstats.mean.neg(unit,ROI,:) = squeeze(mean(ROIsSTMA(unit).neg.(ROIs{ROI}),[1,2]));
        ROIstats.mean.pos(unit,ROI,:) = squeeze(mean(ROIsSTMA(unit).pos.(ROIs{ROI}),[1,2]));
        
        ROIstats.unit(unit,:) = channels2run(unit,:);
          
    end
   
end

% PLOT IT ALLLLLLLL 
[regions] = getChannel2CorticalRegionMapping(subject, collapseAreasFlag);
regions = regions(1:4);
%%% calculate & plot per cell, timept at which peak variance occurred 
%%%% one cell = nROIs x 2 (1 pos, 1 neg sta) 
%%
% cmap1 = [         0.8500 0.3250 0.0980	
%     0.9290 0.6940 0.1250	
%          0.4940 0.1840 0.5560	
%          0.6350 0.0780 0.1840];
cmap1 = [0.4940 0.1840 0.5560
    0.6350 0.0780 0.1840
    0.8500 0.3250 0.0980
    0.9290 0.6940 0.1250];

c = distinguishable_colors(length(regions));
xall = []; 
yall = []; 
colorAll = [];
faceAlphaAll = [];

stas = {'pos','neg'};

fig1 = figure; 

for ROI = 1:length(ROIs)
    
    for ii = 1:length(stas)
        
        for array = 1:length(regions)
             idx2use = ismember(ROIstats.unit(:,1), regions{array}.channels);
             %colorthis = c(array,:);
            colorthis = cmap1(array,:);

            thisArrayPosPeak = squeeze(ROIstats.variance.(stas{ii})(idx2use,ROI,1:end-1))'; %timepts x cells

            [~,peakIdx]= max(thisArrayPosPeak);

            % get rid of STAs that are always below threshold at all timepts
            nonzeroSTAs = any(thisArrayPosPeak,1);
            peakIdx = peakIdx .* nonzeroSTAs;
            peakIdx(peakIdx ==0) = [];

            thisArrayPosPeak = timeAx(peakIdx);

            x = ones(size(thisArrayPosPeak,2),1)*array;
            y = thisArrayPosPeak';

            xall = cat(1,xall, x);
            yall = cat(1,yall, y);
            colorAll = cat(1,colorAll, colorthis);
        end
    end

    
end

beeswarm(xall, yall, 'Colormap', colorAll,'sort_style','hex','corral_style','gutter', 'overlay_style','ci');
title('Per Cell Latencies in Spike-Triggered Maps');
ylim([-0.1 1.1]); ylabel('sec post-spike') 
xticks(1:4);
xticklabels({'S1','M1','PMv','M3'})
%%    
    
%%%% plot variances of positive STAs 
for array = 1:length(regions)
    
    fig1 = figure; 
    
    idx2use = ismember(ROIstats.unit(:,1), regions{array}.channels);
    c = distinguishable_colors(sum(idx2use));
    
    for ROI = 1:length(ROIs)
        thisArrayPosVar = squeeze(ROIstats.variance.pos(idx2use,ROI,:))';
        subplot(3,3,ROI);
        ax = plot(thisArrayPosVar,'LineWidth',3); 
        title(['Variance, PosSTA, array = ' regions{array}.label ' ROI = ' ROIs{ROI}]);
    end
    colororder(c);
    legendCell = cellstr((num2str(ROIstats.unit(idx2use,:))));
    legend(legendCell)
    
end

%%%% plot variances of negative STAs 
for array = 1:length(regions)
    
    fig1 = figure; 
    
    idx2use = ismember(ROIstats.unit(:,1), regions{array}.channels);
    c = distinguishable_colors(sum(idx2use));

    for ROI = 1:length(ROIs)
        
        thisArrayNegVar = squeeze(ROIstats.variance.neg(idx2use,ROI,:))';
        subplot(3,3,ROI);
        plot(thisArrayNegVar,'LineWidth',3);
        title(['Variance, NegSTA, array = ' regions{array}.label ' ROI = ' ROIs{ROI}]);
    end
    
    colororder(c);
    legendCell = cellstr((num2str(ROIstats.unit(idx2use,:))));
    legend(legendCell)
    
end


%%%% plot meanFrames of positive STAs 

for array = 1:length(regions)
    
    fig1 = figure; 
    
    idx2use = ismember(ROIstats.unit(:,1), regions{array}.channels);
    c = distinguishable_colors(sum(idx2use));

    for ROI = 1:length(ROIs)
        
        thisArrayPosMean = squeeze(ROIstats.mean.pos(idx2use,ROI,:))';
        subplot(3,3,ROI);
        plot(thisArrayPosMean,'LineWidth',3);
        title(['MeanFrame, PosSTA, array = ' regions{array}.label ' ROI = ' ROIs{ROI}]);
    end
    
    colororder(c);
    legendCell = cellstr((num2str(ROIstats.unit(idx2use,:))));
    legend(legendCell)
    
end

%%%% plot meanFrames of negative STAs 

for array = 1:length(regions)
    
    fig1 = figure; 
    
    idx2use = ismember(ROIstats.unit(:,1), regions{array}.channels);
    c = distinguishable_colors(sum(idx2use));

    for ROI = 1:length(ROIs)
        
        thisArrayNegMean = squeeze(ROIstats.mean.neg(idx2use,ROI,:))';
        subplot(3,3,ROI);
        plot(thisArrayNegMean,'LineWidth',3);
        title(['MeanFrame, NegSTA, array = ' regions{array}.label ' ROI = ' ROIs{ROI}]);
    end
    
    colororder(c);
    legendCell = cellstr((num2str(ROIstats.unit(idx2use,:))));
    legend(legendCell)
    
end


