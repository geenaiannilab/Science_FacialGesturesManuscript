function [latency,ROIs] = analyzeSTMALatencyROI(date, subject, workDir)

%%%%%%%%% PURPOSE: Analyze potential differences in latencies (timepoints
%%%%%%%%%   of max influence, of cell spiking --> movement) across all cells 
%%%%%%%%
%%%%%%%%  Will download previously generated/defined STAs *only within
%%%%%%%%  particular facial ROIs* &  calculate VARIANCE as well as MEAN
%%%%%%%%  FRAME, inside those ROIs, at each timepoint post-spike. 
%%%%%%%%
%%%%%%%%  (Does exact same thing as analyzeSTMALatencyWholeFrame.m, but
%%%%%%%%    within facial ROIs rather than across whole frame)

%%%%%%%%
%%%%%%%% OUTPUTS: none
%%%%%%%% DEPENDENCIES:  Assumes you already ran generateSTMAbyROI.m to generate 'allCellsSTMAbyROI_(null).mat'

% STMA variables 
null = 'poiss';

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

% load the facial ROIs STAs
load([workDir '/allCellsSTMAbyROI_' null '.mat']);
timeAx = 1/bhvSR:(1/bhvSR):(nFrames/bhvSR);

% for each cell & each ROI, calculate variance & mean frame, for positive &
% negative STAs
for unit = 1:length(ROIsSTMA)
    
    for ROI = 1:length(ROIs)

        posROI= reshape((ROIsSTMA(unit).pos.(ROIs{ROI})), [size(ROIsSTMA(unit).pos.(ROIs{ROI}),1)*size(ROIsSTMA(unit).pos.(ROIs{ROI}),2), size(ROIsSTMA(unit).pos.(ROIs{ROI}),3)]);
        negROI = reshape((ROIsSTMA(unit).neg.(ROIs{ROI})), [size(ROIsSTMA(unit).neg.(ROIs{ROI}),1)*size(ROIsSTMA(unit).neg.(ROIs{ROI}),2), size(ROIsSTMA(unit).neg.(ROIs{ROI}),3)]);
        
        ROIstats.variance.pos(unit,ROI,:) = var(posROI);
        ROIstats.variance.neg(unit,ROI,:) = var(negROI);
        
        ROIstats.mean.neg(unit,ROI,:) = squeeze(mean(ROIsSTMA(unit).neg.(ROIs{ROI}),[1,2]));
        ROIstats.mean.pos(unit,ROI,:) = squeeze(mean(ROIsSTMA(unit).pos.(ROIs{ROI}),[1,2]));
        
        ROIstats.unit(unit,:) = ROIsSTMA(unit).unit;
          
    end
   
end

%%% calculate per cell, timept at which peak variance occurred 
%%%% one cell = nROIs x 2 (1 pos, 1 neg sta) 
%%

for ROI = 1:length(ROIs)


    for unit = 1:length(ROIsSTMA)

        thisCellthisROIpos = squeeze(ROIstats.variance.pos(unit,ROI,1:end-1))'; %timepts x cells
        thisCellthisROIneg = squeeze(ROIstats.variance.neg(unit,ROI,1:end-1))'; %timepts x cells

        [~,peakIdxpos]= max(thisCellthisROIpos);
        [~,peakIdxneg]= max(thisCellthisROIneg);

        % get rid of STAs that are always below threshold at all timepts
%             nonzeroSTAs = any(thisCellthisROI,1);
%             peakIdx = peakIdx .* nonzeroSTAs;
%             peakIdx(peakIdx ==0) = [];

        latency(unit).pos.data(ROI) = timeAx(peakIdxpos);
        latency(unit).neg.data(ROI) = timeAx(peakIdxneg);
        latency(unit).unit =  ROIsSTMA(unit).unit;
       
    end

end
