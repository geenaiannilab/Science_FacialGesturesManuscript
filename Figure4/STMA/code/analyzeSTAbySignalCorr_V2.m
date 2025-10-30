% written GI 
%%%     
%%    only taking cells into account with adequately large STAs, this will calculate 
%% 1) their 2D spatial correlations 2) their signal correlations 3) the relationship between the two metrics 


close all;
clear all;
set(0,'defaultAxesFontSize',16);
screensize = get( groot, 'Screensize' );

date = '210704';
subject ='Barney';
workdir = (['/Users/geena/Dropbox/PhD/SUAinfo/' subject '_' date ]);
regions2analyze = {'All'};
bhvs2plot = [1 2 4];
bhvStrings = {'Threat','Lipsmack','Chew'};
pixThresh = 50 .^2; %5 .^2; % dont include cells with STAs areas smaller than this
option = 'descend'; %'descend' determines if you look at dissimilar or similar pairs first 

% load some pre-defined per-cell metrics 
corrMatrices = load(['/Users/geena/Dropbox/PhD/SUAinfo/Barney_210704/CorrMatricesV2/' regions2analyze{1} '.mat']);

% load two STA summaries / will need to update this 
workdir = '/Users/geena/Dropbox/FaceExpressionAnalysis/STMA2021/Barney_210704/Barney-210704-113250/';
thrROISTAs = load([workdir 'allCellsSTMAbyROI_poiss.mat']);
thrSTAs = load([workdir 'allCellsSTMAacrossTime_poiss_clean']);
thrSTAs.allCellsSTMAacrossTime.both = thrSTAs.allCellsSTMAacrossTime.pos + thrSTAs.allCellsSTMAacrossTime.neg;

workdir = '/Users/geena/Dropbox/FaceExpressionAnalysis/STMA2021/Barney_210704/Barney-210704-110137/';
lsROISTAs = load([workdir 'allCellsSTMAbyROI_poiss.mat']);
lsSTAs = load([workdir 'allCellsSTMAacrossTime_poiss_clean']);
lsSTAs.allCellsSTMAacrossTime.both = lsSTAs.allCellsSTMAacrossTime.pos + lsSTAs.allCellsSTMAacrossTime.neg;

% which cells they have in common
allCorrCells = convertCellLabels(corrMatrices.spikeIdxOriginal{1,1});
[~,indx,~] = intersect(allCorrCells, lsSTAs.allCellsSTMAacrossTime.unit,'rows');
allCorrCells = allCorrCells(indx, :);

% whole STA frames need to be recropped across behaviors 
[lsPixR, lsPixC,~] = size(lsSTAs.allCellsSTMAacrossTime.pos);
[thrPixR, thrPixC,~] = size(thrSTAs.allCellsSTMAacrossTime.pos);
newPixR = min(lsPixR, thrPixR);
newPixC = 530; %min(lsPixC, thrPixC); % to cut out chair artifact in lower right corner of Threat STAs

%  few parameters 
ROIs = lsROISTAs.params.ROIs;
nCells = length(lsROISTAs.ROIsSTMA);


%% first let's just take STAs with adequent #s of pixels above threshold 
counter = 1;
for cc = 1:nCells

    lsStaSize = sum(lsSTAs.allCellsSTMAacrossTime.both(1:newPixR,1:newPixC,cc) > 0,'all');
    thrStaSize = sum(thrSTAs.allCellsSTMAacrossTime.both(1:newPixR,1:newPixR,cc) > 0, 'all');

    if lsStaSize >= pixThresh && thrStaSize >= pixThresh % if enough pixels over threshold 
        trueSTAs.LS.both(:,:,counter) = lsSTAs.allCellsSTMAacrossTime.both(1:newPixR,1:newPixC,cc); % extract LS STA
        trueSTAs.Thr.both(:,:,counter) = thrSTAs.allCellsSTMAacrossTime.both(1:newPixR,1:newPixC,cc); % extract Thr STA 
        trueSTAs.unit(counter,:) = lsSTAs.allCellsSTMAacrossTime.unit(cc,:); % unit label 

        for rr = 1:length(ROIs) % also binarize by ROIs
            trueSTAs.LS.binarized(rr,counter) =  double(sum((sum(lsROISTAs.ROIsSTMA(cc).pos.(ROIs{rr}),3) ~= 0 ),'all') >= pixThresh);
            trueSTAs.Thr.binarized(rr,counter) =  double(sum((sum(thrROISTAs.ROIsSTMA(cc).pos.(ROIs{rr}),3) ~= 0 ),'all') >= pixThresh);
        end
        counter = counter +1;
    end

end

% now reset nCells to new #
nCells = length(trueSTAs.unit);

% for each cell pair, for each behavior, count the number of ROIs they have
% in common (ie., both 1s in 'mouth'; do NOT count two ROIs, both without
% STAs, as being "similar) 
% this is a super coarse measure of how similar the cell pairs' STAs are
for cc = 1:nCells
    for dd = 1:nCells
        corrLSROIbinarized(cc,dd) = sum((trueSTAs.LS.binarized(:,cc) + trueSTAs.LS.binarized(:,dd)) == 2) ./ length(ROIs);
        corrThrROIbinarized(cc,dd) =sum((trueSTAs.Thr.binarized(:,cc) + trueSTAs.Thr.binarized(:,dd)) == 2) ./ length(ROIs);
    end

end


% now only take cells left in common 
[~,indx,~] = intersect(allCorrCells, trueSTAs.unit,'rows');
allCorrCells = allCorrCells(indx, :);

% calc 2D spatial correlation for all cell pairs, for both bhvs 
% this is SPATIAL correlation only, ie., does not take magnitude into
% account (see " thisSTA > 0 " below 

corrLSSTAs = zeros(nCells, nCells);
corrThrSTAs = zeros(nCells, nCells);

for cc = 1:nCells
    for dd = 1:nCells
        thisCellLS = trueSTAs.LS.both(:,:,cc) > 0 ; 
        thatCellLS = trueSTAs.LS.both(:,:,dd) > 0 ;
        corrLSSTAs(cc,dd) = corr2(thisCellLS,thatCellLS );
        
        thisCellThr = trueSTAs.Thr.both(:,:,cc) > 0 ;
        thatCellThr = trueSTAs.Thr.both(:,:,dd) > 0 ;
        corrThrSTAs(cc,dd) = corr2(thisCellThr, thatCellThr);

    end

end

%% extract the signalCorr, noiseCorr, & normalized PETHs of all cellpairs in each behavior 
for bhv = 1:length(bhvs2plot)

    allCorr(:,:,bhv) = corrMatrices.signalCorrOverTime(bhv).corr(indx,indx) - eye(size(corrMatrices.signalCorrOverTime(bhv).corr(indx,indx)));        
    noiseCorr(:,:,bhv) = corrMatrices.noiseCorr(bhv).corr(indx,indx) - eye(size(corrMatrices.noiseCorr(bhv).corr(indx,indx)));        
    normTrialAvg(bhv).data = (corrMatrices.trialAvg(bhv).data(:,indx));

end


%% rank cell pairs by degree of dissimilarity in signal correlation 
diffAllCorr = abs(allCorr(:,:,2) - allCorr(:,:,1)); % 1 - 1 = 0 (0 if correlated in both bhvs)
                                               % -1 - -1 = 0 (0 if anticorrelated in both bhvs)
                                               % 1 - -1 = 2 (if switch correlations)
                                               % -1 - 1 = -2 (if switch correlations)

[diffAllCorrSorted, diffAllCorrSortedCellIDpairs, sortedOrigIndx] = rankPairs(diffAllCorr,allCorrCells,option);
sortedDiffSignalCorr.data = diffAllCorrSorted;
sortedDiffSignalCorr.cellIDPairs = diffAllCorrSortedCellIDpairs;

%% reorder cell pairs' value of whole STA correlation  AND 2) degree of STA ROI overlap 
%% for each behavior, according to rank-ordering determined by Diff signalCorr (LSsignalcorr - Thr Signal Corr) 

for i = 1:length(sortedOrigIndx)
    thisR = sortedOrigIndx(i,1);
    thisC = sortedOrigIndx(i,2);

    sortedCorrLSSTAs(i) = corrLSSTAs(thisR, thisC);
    sortedCorrLSROIbinarized(i) = corrLSROIbinarized(thisR, thisC);
    sortedCorrThrSTAs(i) = corrThrSTAs(thisR, thisC);
    sortedCorrThrROIbinarized(i) = corrThrROIbinarized(thisR,thisC);
    sortedThrNoiseCorr(i) = noiseCorr(thisR,thisC,1);
    sortedLSNoiseCorr(i) = noiseCorr(thisR,thisC,2);

end



%% plot most extreme cell-pairs' STAs (whole & ROI-based), along with their correlation coefficient measure 
plotFig3 = figure;  
set(plotFig3, 'Position', [screensize(3) screensize(4) screensize(3) screensize(4)]);

for pp = 1:2:length(sortedOrigIndx) % skips the cell A-B / cell B-A repeat

    cc = sortedOrigIndx(pp,1);
    dd = sortedOrigIndx(pp,2);

        plotFig3; 
        ha = tight_subplot(2,3,[.05 .05],[.02 .02],[.01 .01]);
        
        axes(ha(1)); 
        imagesc(lsSTAs.allCellsSTMAacrossTime.pos(1:newPixR,1:newPixC,cc )); colormap(jet)
        ha(1).XTick = []; ha(1).YTick = [];
        title(['Unit = ' num2str(lsSTAs.allCellsSTMAacrossTime.unit(cc,:))])

        axes(ha(2)); 
        imagesc(lsSTAs.allCellsSTMAacrossTime.pos(1:newPixR,1:newPixC,dd)); colormap(jet); 
        ha(2).XTick = []; ha(2).YTick = [];
        title(['Unit = ' num2str(lsSTAs.allCellsSTMAacrossTime.unit(dd,:)) '  LS STAs Corr = ' num2str(corrLSSTAs(cc,dd))]);
 
  axes(ha(3));
        plot(corrMatrices.taxis, normTrialAvg(2).data(:,cc),'linew',6); hold on;
        plot(corrMatrices.taxis, normTrialAvg(2).data(:,dd),'linew',6); hold off;
        title([' DiffCorr  = ' num2str(sortedDiffSignalCorr.data(pp)) '; LS PETH: Corr = ' num2str(allCorr(cc,dd,2))]);
        %m = sprintf(['LS NoiseCorr=\n' num2str(sortedLSNoiseCorr(pp))]); text(0, 0,m,'FontSize',12);

        axes(ha(4));
        imagesc(thrSTAs.allCellsSTMAacrossTime.pos(1:newPixR,1:newPixC,cc )); colormap(jet)
        ha(4).XTick = []; ha(4).YTick = [];
        title(['Unit = ' num2str(thrSTAs.allCellsSTMAacrossTime.unit(cc,:))])

        axes(ha(5));
        imagesc(thrSTAs.allCellsSTMAacrossTime.pos(1:newPixR,1:newPixC,dd)); colormap(jet); 
        ha(5).XTick = []; ha(5).YTick = [];
        title(['Unit = ' num2str(thrSTAs.allCellsSTMAacrossTime.unit(dd,:)) '  Thr STAs Corr = ' num2str(corrThrSTAs(cc,dd))]);
        
       
        axes(ha(6))
        plot(corrMatrices.taxis, normTrialAvg(1).data(:,cc),'linew',6); hold on;
        plot(corrMatrices.taxis, normTrialAvg(1).data(:,dd),'linew',6); hold off;
        title([' DiffCorr  = ' num2str(sortedDiffSignalCorr.data(pp)) '; Thr PETH: Corr = ' num2str(allCorr(cc,dd,1))]);
        %m = sprintf(['Thr NoiseCorr=\n' num2str(sortedThrNoiseCorr(pp))]); text(0, 0,m,'FontSize',12);

        pause; clf;
end


%%



function [cellLabelsOut] = convertCellLabels(cellLabelsIn)
    
   
    cellLabelsOut = cellfun(@(x) str2double(strsplit(x, '-')), cellLabelsIn, 'UniformOutput', false);
    cellLabelsOut = vertcat(cellLabelsOut{:});

end


     
