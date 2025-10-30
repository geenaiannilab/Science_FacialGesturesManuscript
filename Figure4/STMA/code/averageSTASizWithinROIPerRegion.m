% written GI 
%%%     
%%    count fraction of cells per region that have STAs vs not (based on one ROI metric) 
%%  & calculate per cell STA size (LS and thr separate), then plot 

close all;
clear all;
set(0,'defaultAxesFontSize',18);
screensize = get( groot, 'Screensize' );

%% Load matfiles
pathToMatFiles = '/Users/geena/Dropbox/FaceExpressionAnalysis/STMA2021';

thorData.Thr = load(fullfile(pathToMatFiles, '/Thor_171010/Thor-171010-150131/allCellsSTMAbyROI_poiss.mat'));
thorData.LS = load(fullfile(pathToMatFiles, '/Thor_171010/Thor-171010-154538/allCellsSTMAbyROI_poiss.mat'));
[thorRegions] = getChannel2CorticalRegionMapping('Thor', 1);
thorPixThresh = 25 .^2;% roughly 25 x 25 p is a reasonable activation in thor's camera but this is entirely arbitrary

barneyData.Thr = load(fullfile(pathToMatFiles, '/Barney_210704/Barney-210704-113250/allCellsSTMAbyROI_poiss.mat'));
barneyData.LS = load(fullfile(pathToMatFiles, '/Barney_210704/Barney-210704-110137/allCellsSTMAbyROI_poiss.mat'));
[barneyRegions] = getChannel2CorticalRegionMapping('Barney', 1);
barneyPixThresh = 50 .^2; % roughly 50 x 50 p is a reasonable activation in barneys's camera, but entirely arbitrary 

%%  a few parameters 
ROIs = barneyData.Thr.params.ROIs;
nCellsBarney = length(barneyData.LS.ROIsSTMA);
nCellsThor = length(thorData.LS.ROIsSTMA);

% for each cell, and each behavior, count the number of above threshold
% pixels inside each ROI ("bhvROIsum" is pixel count inside each ROI (continuous), "bhvROIbinarized" is 0s/1s if ROI had adequate # pixels) 
for cc = 1:nCellsBarney

    for rr = 1:length(ROIs)

        barneyResults.lsROIsum(rr,cc) =  sum(barneyData.LS.ROIsSTMA(cc).pos.(ROIs{rr}),'all'); % for each ROI & each cell, sum of all positive pixels 
        barneyResults.lsROIbinarized(rr,cc) = double(barneyResults.lsROIsum(rr,cc) >= barneyPixThresh); % for each ROI & each cell, 1/0 if adequate pixels inside or not
        barneyResults.thrROIsum(rr,cc) = sum(barneyData.Thr.ROIsSTMA(cc).pos.(ROIs{rr}),'all'); % for each ROI & each cell, sum of all positive pixels 
        barneyResults.thrROIbinarized(rr,cc) = double(barneyResults.thrROIsum(rr,cc) >= barneyPixThresh);

    end

end
ids = [barneyData.LS.ROIsSTMA(:).unit]; % cell channel IDs
barneyResults.cellID = reshape(ids,2,(length(ids)/2))';

% repeat for thor
for cc = 1:nCellsThor

    for rr = 1:length(ROIs)

        thorResults.lsROIsum(rr,cc) =  sum(thorData.LS.ROIsSTMA(cc).pos.(ROIs{rr}),'all'); % for each ROI & each cell, sum of all positive pixels 
        thorResults.lsROIbinarized(rr,cc) = double(thorResults.lsROIsum(rr,cc) >= thorPixThresh); % for each ROI & each cell, 1/0 if adequate pixels inside or not
        thorResults.thrROIsum(rr,cc) = sum(thorData.Thr.ROIsSTMA(cc).pos.(ROIs{rr}),'all'); % for each
        thorResults.thrROIbinarized(rr,cc) = double(thorResults.thrROIsum(rr,cc) >= thorPixThresh);

    end

end
ids = [thorData.LS.ROIsSTMA(:).unit]; % cell IDs
thorResults.cellID = reshape(ids,2,(length(ids)/2))';

% for each behavior and each region, calc fraction of cells with 
bhvs = {'Thr','LS'};

for bhv = 1:length(bhvs)

    cellLabels = barneyResults.cellID;
    for rr = 1:length(barneyRegions)

        regionIndices = find(cellLabels(:,1) >= min(barneyRegions{rr}.channels) & cellLabels(:,1) <= max(barneyRegions{rr}.channels));
        barneyResults.(barneyRegions{rr}.label).(bhvs{bhv}).fractionCellsPosSTAs = sum(barneyResults.lsROIbinarized(:,regionIndices),2)' ./ length(regionIndices); % for each ROI, fraction of cells in region w/ positive STA at that location 
        barneyCount2plot(:,rr,bhv) = barneyResults.(barneyRegions{rr}.label).(bhvs{bhv}).fractionCellsPosSTAs; % for each ROI, fraction of cells in region w/ positive STA at that location 

    end
end


% repeat for Thor 
for bhv = 1:length(bhvs)

    cellLabels = thorResults.cellID;
    for rr = 1:length(thorRegions)

        regionIndices = find(cellLabels(:,1) >= min(thorRegions{rr}.channels) & cellLabels(:,1) <= max(thorRegions{rr}.channels));
        thorResults.(thorRegions{rr}.label).(bhvs{bhv}).fractionCellsPosSTAs = sum(thorResults.lsROIbinarized(:,regionIndices),2)' ./ length(regionIndices); % % for each ROI, fraction of cells in region w/ positive STA at that location 
        thorCount2plot(:,rr,bhv) = thorResults.(thorRegions{rr}.label).(bhvs{bhv}).fractionCellsPosSTAs; % for each ROI, fraction of cells in region w/ positive STA at that location 
    end
end


% collect size of STA inside each ROI, per behavior 
cellLabels = barneyResults.cellID;
for rr = 1:length(barneyRegions)

    regionIndices = find(cellLabels(:,1) >= min(barneyRegions{rr}.channels) & cellLabels(:,1) <= max(barneyRegions{rr}.channels));
    
    for ROI = 1:length(ROIs)
        idx = barneyResults.lsROIsum(ROI,regionIndices) >= thorPixThresh;
        data = barneyResults.lsROIsum(ROI,regionIndices);
        barneySize2plotLS{ROI,rr} = data(idx);
        
        idx = barneyResults.thrROIsum(ROI,regionIndices) >= thorPixThresh;
        data = barneyResults.thrROIsum(ROI,regionIndices);  
        barneySize2plotThr{ROI,rr} = data(idx);
    end

end

% repeat for thor 
cellLabels = thorResults.cellID;
for rr = 1:length(thorRegions)

    regionIndices = find(cellLabels(:,1) >= min(thorRegions{rr}.channels) & cellLabels(:,1) <= max(thorRegions{rr}.channels));
    
    for ROI = 1:length(ROIs)
        idx = thorResults.lsROIsum(ROI,regionIndices) >= thorPixThresh;
        data = thorResults.lsROIsum(ROI,regionIndices);
        thorSize2plotLS{ROI,rr} = data(idx);
        
        idx = thorResults.thrROIsum(ROI,regionIndices) >= thorPixThresh;
        data = thorResults.thrROIsum(ROI,regionIndices);  
        thorSize2plotThr{ROI,rr} = data(idx);
    end

end

% for each ROI, plot fraction of cells in region w/ positive STA at that location 
figure; 
subplot 121; b1 = bar(barneyCount2plot(:,1:4,1)); title('B, Thr');
axes= gca; axes.XTick = 1:length(ROIs); axes.XTickLabel = ROIs; legend('S1','M1','PMv','M3');
subplot 122; b2 = bar(barneyCount2plot(:,1:4,2)); title('B, LS');
axes= gca; axes.XTick = 1:length(ROIs); axes.XTickLabel = ROIs; 

figure; 
subplot 121; b3 = bar(thorCount2plot(:,1:4,1)); title('T, Thr');
axes= gca; axes.XTick = 1:length(ROIs); axes.XTickLabel = ROIs; legend('S1','M1','PMv','M3');
subplot 122; b4 = bar(thorCount2plot(:,1:4,2)); title('T, LS');
axes= gca; axes.XTick = 1:length(ROIs); axes.XTickLabel = ROIs; 


% plot size of STA within each ROI, per behavior, per animal 
% barney
cmap1 = linspecer((length(barneyRegions)-1),'sequential'); %   'sequential'

for ROI = 1:length(ROIs)
     xall = [];
    yall = [];
    colorAll = [];
    counter = 1;

    for ii = 1:length(barneyRegions)-1
    
    color = cmap1(counter,:);
    
    y = barneySize2plotThr{ROI,ii}';
    x = ones(length(y),1)*(counter);

        xall = cat(1,xall, x);
        yall = cat(1,yall, y);
        colorAll = cat(1,colorAll, color);
       

       counter = counter +1;
       
    end
    plotFig = figure;
    b = beeswarm(xall,yall,'Colormap', colorAll,'dot_size',1.2,'overlay_style','ci');
    axes = gca; axes.XTick=1:4; axes.XTickLabel={'S1','M1','PMv','M3'}; xlabel('region')
    title(['Barney - STA size in ' ROIs(ROI) ])

end

% plot size of STA within each ROI, per behavior, per animal 
% thor
for ROI = 1:length(ROIs)
     xall = [];
    yall = [];
    colorAll = [];
    counter = 1;

    for ii = 1:length(thorRegions)-1
    
    color = cmap1(counter,:);
    
    y = thorSize2plotThr{ROI,ii}';
    x = ones(length(y),1)*(counter);

        xall = cat(1,xall, x);
        yall = cat(1,yall, y);
        colorAll = cat(1,colorAll, color);
       

       counter = counter +1;
       
    end
    plotFig = figure;
    b = beeswarm(xall,yall,'Colormap', colorAll,'dot_size',1.2,'overlay_style','ci');
    axes = gca; axes.XTick=1:4; axes.XTickLabel={'S1','M1','PMv','M3'}; xlabel('region')
    title(['Thor - STA size in ' ROIs(ROI) ])

end


