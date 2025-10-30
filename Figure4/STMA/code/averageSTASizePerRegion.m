% written GI 
%%%     
%%    count fraction of cells per region that have STAs vs not (based on one whole-frame metric) 
%  & calculate per cell STA size (LS and thr separate) 

close all;
clear all;
set(0,'defaultAxesFontSize',18);
screensize = get( groot, 'Screensize' );

%% Load matfiles
pathToMatFiles = '/Users/geena/Dropbox/FaceExpressionAnalysis/STMA2021';

thorData.Thr = load(fullfile(pathToMatFiles, '/Thor_171010/Thor-171010-150131/allCellsSTMAacrossTime_poiss_highThresh_clean.mat'));
thorData.LS = load(fullfile(pathToMatFiles, '/Thor_171010/Thor-171010-154538/allCellsSTMAacrossTime_poiss_highThresh_clean.mat'));
[thorRegions] = getChannel2CorticalRegionMapping('Thor', 1);
thorPixThresh = 25 .^2;% roughly 25 x 25 p is a reasonable activation in thor's camera 

barneyData.Thr = load(fullfile(pathToMatFiles, '/Barney_210704/Barney-210704-113250/allCellsSTMAacrossTime_poiss_clean.mat'));
barneyData.LS = load(fullfile(pathToMatFiles, '/Barney_210704/Barney-210704-110137/allCellsSTMAacrossTime_poiss_clean.mat'));
[barneyRegions] = getChannel2CorticalRegionMapping('Barney', 1);
barneyPixThresh = 50 .^2; % roughly 50 x 50 p is a reasonable activation in barneys's camera 

%%  few parameters 
nLSCells = length(barneyData.LS.allCellsSTMAacrossTime.unit);
nThrCells = length(barneyData.Thr.allCellsSTMAacrossTime.unit);
%ROIs = barneyData.barneyThrROIs.params.ROIs;
if nLSCells ~= nThrCells
    disp('LS & Thr STAs not computed for all same neurons!')
    return;
else
    nCells = nLSCells;
end

%% count number of cells with significantly (>thresh) positive/negative STAs
% at a given pixel; 

bhvs = {'Thr','LS'};


for bhv = 1:length(bhvs)

    maskPos = barneyData.(bhvs{bhv}).allCellsSTMAacrossTime.pos ~= 0 ; % all above thresh ++ pixels
    maskNeg = barneyData.(bhvs{bhv}).allCellsSTMAacrossTime.neg ~= 0 ; % all above thresh -- pixels 

    cellLabels = barneyData.(bhvs{bhv}).allCellsSTMAacrossTime.unit;

    for rr = 1:length(barneyRegions)

        regionIndices = find(cellLabels(:,1) >= min(barneyRegions{rr}.channels) & cellLabels(:,1) <= max(barneyRegions{rr}.channels));
        
        sizePosTmp = squeeze(sum(maskPos(:,:,regionIndices),[1 2]));
        sizeNegTmp = squeeze(sum(maskNeg(:,:,regionIndices),[1 2]));
        sizeAnyTmp = sizePosTmp + sizeNegTmp;

        barneyResults.(barneyRegions{rr}.label).(bhvs{bhv}).maskPos = maskPos(:,:,regionIndices);
        barneyResults.(barneyRegions{rr}.label).(bhvs{bhv}).sizePos =  sizePosTmp(sizePosTmp >= barneyPixThresh);
        barneyResults.(barneyRegions{rr}.label).(bhvs{bhv}).maskNeg = maskNeg(:,:,regionIndices);
        barneyResults.(barneyRegions{rr}.label).(bhvs{bhv}).sizeNeg =  sizeNegTmp(sizeNegTmp >= barneyPixThresh);
        barneyResults.(barneyRegions{rr}.label).(bhvs{bhv}).maskAny = maskPos(:,:,regionIndices) | maskNeg(:,:,regionIndices);
        barneyResults.(barneyRegions{rr}.label).(bhvs{bhv}).sizeAny = sizeAnyTmp(sizeAnyTmp >= barneyPixThresh);

        barneyResults.(barneyRegions{rr}.label).(bhvs{bhv}).totalNumberOfCells = length(regionIndices);

        
    end
end


%% repeat for Thor 

for bhv = 1:length(bhvs)

    maskPos = thorData.(bhvs{bhv}).allCellsSTMAacrossTime.pos ~= 0 ; % all above thresh ++ pixels
    maskNeg = thorData.(bhvs{bhv}).allCellsSTMAacrossTime.neg ~= 0 ; % all above thresh -- pixels 

    cellLabels = thorData.(bhvs{bhv}).allCellsSTMAacrossTime.unit;

    for rr = 1:length(thorRegions)

        regionIndices = find(cellLabels(:,1) >= min(thorRegions{rr}.channels) & cellLabels(:,1) <= max(thorRegions{rr}.channels));
        
        sizePosTmp = squeeze(sum(maskPos(:,:,regionIndices),[1 2]));
        sizeNegTmp = squeeze(sum(maskNeg(:,:,regionIndices),[1 2]));
        sizeAnyTmp = sizePosTmp + sizeNegTmp;

        thorResults.(thorRegions{rr}.label).(bhvs{bhv}).maskPos = maskPos(:,:,regionIndices);
        thorResults.(thorRegions{rr}.label).(bhvs{bhv}).sizePos =  sizePosTmp(sizePosTmp >= barneyPixThresh);
        thorResults.(thorRegions{rr}.label).(bhvs{bhv}).maskNeg = maskNeg(:,:,regionIndices);
        thorResults.(thorRegions{rr}.label).(bhvs{bhv}).sizeNeg =  sizeNegTmp(sizeNegTmp >= barneyPixThresh);
        thorResults.(thorRegions{rr}.label).(bhvs{bhv}).maskAny = maskPos(:,:,regionIndices) | maskNeg(:,:,regionIndices);
        thorResults.(thorRegions{rr}.label).(bhvs{bhv}).sizeAny = sizeAnyTmp(sizeAnyTmp >= barneyPixThresh);

        thorResults.(thorRegions{rr}.label).(bhvs{bhv}).totalNumberOfCells = length(regionIndices);

        
    end
end

%%
% plot STA count & size across regions, per subject  
regions = {'S1','M1','PMv','M3','All'};
bhvs = {'Thr','LS'};

% concatente for plotting 
cmap1 = [         
    0.4940 0.1840 0.5560	
    0.6350 0.0780 0.1840
    0.8500 0.3250 0.0980	
    0.9290 0.6940 0.1250	
         0 0 0 ];

% count fraction of cells in each region w/ STAs 
for bhv = 1:length(bhvs)

    for ii = 1:length(regions)
    
            thisRegionB = barneyResults.(regions{ii});
            thisRegionT = thorResults.(regions{ii});

            barneyCellsWithSTAs(bhv,ii) = length(thisRegionB.(bhvs{bhv}).sizeAny) ./ thisRegionB.(bhvs{bhv}).totalNumberOfCells;
            thorCellsWithSTAs(bhv,ii) = length(thisRegionT.(bhvs{bhv}).sizeAny) ./ thisRegionT.(bhvs{bhv}).totalNumberOfCells;
            bothAnimalsCellsWithSTAs(bhv,ii) = (length(thisRegionB.(bhvs{bhv}).sizeAny) + length(thisRegionT.(bhvs{bhv}).sizeAny)) ./ ...
                (thisRegionB.(bhvs{bhv}).totalNumberOfCells + thisRegionT.(bhvs{bhv}).totalNumberOfCells);
    end

end

for ii = 1:length(regions)

    thisRegionB = barneyResults.(regions{ii});
    thisRegionT = thorResults.(regions{ii});

    barneyCombined(ii,:) = (length(thisRegionB.Thr.sizeAny) + length(thisRegionB.LS.sizeAny)) ./ (thisRegionB.Thr.totalNumberOfCells + thisRegionB.LS.totalNumberOfCells);
    thorCombined(ii,:) = (length(thisRegionT.Thr.sizeAny) + length(thisRegionT.LS.sizeAny)) ./ (thisRegionT.Thr.totalNumberOfCells + thisRegionT.LS.totalNumberOfCells);
   
    allCombined(ii,:) = (length(thisRegionB.Thr.sizeAny) + length(thisRegionT.Thr.sizeAny) + length(thisRegionB.LS.sizeAny) + length(thisRegionT.LS.sizeAny)) ./ ...
        (thisRegionB.Thr.totalNumberOfCells + thisRegionT.Thr.totalNumberOfCells + thisRegionB.LS.totalNumberOfCells + thisRegionT.LS.totalNumberOfCells);
end

figure; 
b1 = bar([barneyCombined thorCombined],'FaceColor','flat') ;

for k = 1:length(regions)
    b1(1).CData(k,:) = cmap1(k,:);
    b1(2).CData(k,:) = cmap1(k,:);
    b1(2).FaceAlpha = 0.5;
end
ylabel('Fraction of Cells '); ylim([0 1]);
title('Fraction Cells Per Region w/ STAs');
set(gca, 'XTick', 1:5); set(gca, 'XTickLabel', regions);

figure; 
b1 = bar(allCombined','FaceColor','flat') ;

for k = 1:length(regions)
    b1(1).CData(k,:) = cmap1(k,:);
    
end
ylabel('Fraction of Cells '); ylim([0 1]);
title('Fraction Cells Per Region w/ STAs');
set(gca, 'XTick', 1:5); set(gca, 'XTickLabel', regions);


figure; 
b1 = bar(bothAnimalsCellsWithSTAs','FaceColor','flat'); 

for k = 1:length(regions)
    b1(1).CData(k,:) = cmap1(k,:);
    
    b1(2).CData(k,:) = cmap1(k,:);
end
b1(1).EdgeColor = 'r'; b1(1).LineWidth = 2;
b1(2).EdgeColor = 'b';b1(2).LineWidth = 2;
ylabel('Fraction Cells Per Region w/ STAs '); ylim([0 1]);
title(['Fraction Cells Per Region w/ STAs in Threat (red) or LS (Blue)']);
set(gca, 'XTick', 1:5);
set(gca, 'XTickLabel', regions);


figure; 
b1 = bar(thorCellsWithSTAs','FaceColor','flat'); 

for k = 1:length(regions)
    b1(1).CData(k,:) = cmap1(k,:);
    
    b1(2).CData(k,:) = cmap1(k,:);
end
b1(1).EdgeColor = 'r'; b1(1).LineWidth = 2;
b1(2).EdgeColor = 'b';b1(2).LineWidth = 2;
ylabel('Fraction Cells Per Region w/ STAs '); ylim([0 1]);
title(['Thor: Fraction Cells Per Region w/ STAs in Threat (red) or LS (Blue)']);
set(gca, 'XTick', 1:5);
set(gca, 'XTickLabel', regions);


figure; 
b1 = bar(barneyCellsWithSTAs','FaceColor','flat'); 

for k = 1:length(regions)
    b1(1).CData(k,:) = cmap1(k,:);
    
    b1(2).CData(k,:) = cmap1(k,:);
end
b1(1).EdgeColor = 'r'; b1(1).LineWidth = 2;
b1(2).EdgeColor = 'b';b1(2).LineWidth = 2;
ylabel('Fraction Cells Per Region w/ STAs '); ylim([0 1]);
title(['Barney: Fraction Cells Per Region w/ STAs in Threat (red) or LS (Blue)']);
set(gca, 'XTick', 1:5);
set(gca, 'XTickLabel', regions);

% barney, size
regions = {'S1','M1','PMv','M3'};
for bhv = 1:length(bhvs)
    
    xall = []; yall = []; colorAll = []; faceAlphaAll = [];
    counter = 1; 

    for ii = 1:length(regions)
    
            thisRegion = barneyResults.(regions{ii});
            color = cmap1(ii,:);
            faceAlpha = 0.8;
        
            x = ones(length(thisRegion.(bhvs{bhv}).sizeAny),1)*(counter);
            y = thisRegion.(bhvs{bhv}).sizeAny;
    
            xall = cat(1,xall, x);
            yall = cat(1,yall, y);
            colorAll = cat(1,colorAll, color);
            faceAlphaAll = cat(1,faceAlphaAll, faceAlpha);
           
            counter = counter +1;
    
    end
    
    plotFig = figure;
    b = beeswarm(xall,yall,'Colormap', colorAll,'corral_style','rand','sort_style','fan','dot_size',1.2,'MarkerFaceAlpha',faceAlphaAll','overlay_style','ci');
    xticks([1 2 3 4 ]); xlim([0 5])
    xticklabels({'S1', 'M1','PMv', 'M3'})
    xlabel('Region')
    ylabel('STA sizes');
    title([ bhvs{bhv} ' STA Size, all Cells per Regions; BARNEY'])

     plotFig = figure;
    boxplot(yall,xall)
    xticks([1 2 3 4 ]); xlim([0 5])
    xticklabels({'S1', 'M1','PMv', 'M3'})
    xlabel('Region')
    ylabel('STA sizes');
    title([ bhvs{bhv} ' STA Size, all Cells per Regions; BARNEY'])
end

% repeat for thor, size

for bhv = 1:length(bhvs)
    
    xall = []; yall = []; colorAll = []; faceAlphaAll = [];
    counter = 1; 

    for ii = 1:length(regions)
    
            thisRegion = thorResults.(regions{ii});
            color = cmap1(ii,:);
            faceAlpha = 0.8;
        
            x = ones(length(thisRegion.(bhvs{bhv}).sizeAny),1)*(counter);
            y = thisRegion.(bhvs{bhv}).sizeAny;
    
            xall = cat(1,xall, x);
            yall = cat(1,yall, y);
            colorAll = cat(1,colorAll, color);
            faceAlphaAll = cat(1,faceAlphaAll, faceAlpha);
           
            counter = counter +1;
    
    end
    
    plotFig = figure;
    b = beeswarm(xall,yall,'Colormap', colorAll,'corral_style','rand','sort_style','fan','dot_size',1.2,'MarkerFaceAlpha',faceAlphaAll','overlay_style','ci');
    xticks([1 2 3 4 ]); xlim([0 5])
    xticklabels({'S1', 'M1','PMv', 'M3'})
    xlabel('Region')
    ylabel('STA sizes');
    title([ bhvs{bhv} ' STA Size, all Cells per Regions; THOR'])

     plotFig = figure;
    boxplot(yall,xall)
    xticks([1 2 3 4 ]); xlim([0 5])
    xticklabels({'S1', 'M1','PMv', 'M3'})
    xlabel('Region')
    ylabel('STA sizes');
    title([ bhvs{bhv} ' STA Size, all Cells per Regions; THOR'])
end

regions = {'S1','M1','PMv','M3'};
for bhv = 1:length(bhvs)
    
    xall = []; yall = []; colorAll = []; faceAlphaAll = [];
    counter = 1; 

    for ii = 1:length(regions)
    
            thisRegionB = barneyResults.(regions{ii});
            thisRegionT = thorResults.(regions{ii});

            allData = [thisRegionB.(bhvs{bhv}).sizeAny ; thisRegionT.(bhvs{bhv}).sizeAny];
            color = cmap1(ii,:);
            faceAlpha = 0.8;
        
            x = ones(length(allData),1)*(counter);
            y = allData;
    
            xall = cat(1,xall, x);
            yall = cat(1,yall, y);
            colorAll = cat(1,colorAll, color);
            faceAlphaAll = cat(1,faceAlphaAll, faceAlpha);
           
            counter = counter +1;
    
    end
    
    plotFig = figure;
    b = beeswarm(xall,yall,'Colormap', colorAll,'corral_style','rand','sort_style','fan','dot_size',1.2,'MarkerFaceAlpha',faceAlphaAll','overlay_style','ci');
    xticks([1 2 3 4 ]); xlim([0 5])
    xticklabels({'S1', 'M1','PMv', 'M3'})
    xlabel('Region')
    ylabel('STA sizes');
    title([ bhvs{bhv} ' STA Size, all Cells per Regions; Monkey B and T'])
end
