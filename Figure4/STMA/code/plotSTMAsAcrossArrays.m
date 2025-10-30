%%%%%%%%%
%%%%%%%%% PURPOSE: once you have run collapseSTMAsAcrossCells., 
%%%%%%%%%   this script will plot the *spatial STAs* of all cells,
%%%%%%%%    in order to test if cells of a given cortical region have spatially-distinct STAs. 
%%%%%%%% Per array, this will generate plots where the color indicates % of cells at that spatial location,
%%%%%%%%        with a modulated STA 
%%%%%%%%  Per array, this will plot both the excitatory/positive &
%%%%%%%%            inhibitory/negative aspects the STAs 
%%%%%%%%
%%%%%%%%  DEPENDENCIES: Assumes you already ran collapseSTMAsAcrossCells.m to generate 'allCellsSTMAAcrossTime_(null).mat'

clear all; 
close all;
set(0,'defaultAxesFontSize',20)

date = '210704';
subject ='Barney';
sess2run = 5;
null = 'poiss'; % poiss, rand, or randLong 
cleanSTAflag = 1;

rawCount = 0; %  plot raw # of cells per location, or % of total 
collapseArrays = 1; %  plot all M1s, all PMvs, and all M3s together, or each array separately 
saveFlag = 1;

% load bhv data
bhvDir = (['/Users/geena/Dropbox/PhD/SUAinfo/' subject '_' date '/Data4Analysis']);
filePattern = fullfile(bhvDir, ['bhvspikes_sub' num2str(sess2run) '*']);
bhvFile = dir(filePattern);
TDTsessID = cell2mat(extractBetween(bhvFile.name,['sub' num2str(sess2run) '_'],['_']));
bhvFile = bhvFile.name;

% set up files to load 
STMAdir = ['/Users/geena/Dropbox/FaceExpressionAnalysis/STMA2021/' subject '_' date];
workDir = [STMAdir '/' subject '-' date '-' TDTsessID '/'];

if cleanSTAflag
    f2load = [workDir 'allCellsSTMAAcrossTime_' null '_clean.mat'];
elseif ~cleanSTAflag
    f2load = [workDir 'allCellsSTMAacrossTime_' null '.mat'];
end

if exist(f2load,'file')
    load(f2load);
else
    return
end

% count number of cells with significantly (>thresh) positive/negative STAs
% at a given pixel; then plot 
maskPos = allCellsSTMAacrossTime.pos ~= 0 ;
maskNeg = allCellsSTMAacrossTime.neg ~= 0 ;
allCellsSTMAacrossTime.pos(maskPos) = 1;
allCellsSTMAacrossTime.neg(maskNeg) = 1;

[regions] = getChannel2CorticalRegionMapping(subject, collapseArrays);

for array = 1:length(regions)
    fig = figure('units','normalized','outerposition',[0 0 1 1]); 
     if array == 1 
         colorAxE = [0 0.6];
     else
         colorAxE = [0 0.6];
     end
     colorAxI = [0 0.2];
    
    idx2use = ismember(allCellsSTMAacrossTime.unit(:,1), regions{array}.channels);
    
    if rawCount
        allNeuronSTMApos2plot = sum(allCellsSTMAacrossTime.pos(:,:,idx2use),3);
        allNeuronSTMAneg2plot = sum(allCellsSTMAacrossTime.neg(:,:,idx2use),3);
        %colorAx = [0 max(allNeuronSTMApos2plot,[],'all')];
    else
         allNeuronSTMApos2plot = (sum(allCellsSTMAacrossTime.pos(:,:,idx2use),3)) ./ sum(idx2use);
        allNeuronSTMAneg2plot =  (sum(allCellsSTMAacrossTime.neg(:,:,idx2use),3)) ./ sum(idx2use) ;
        cMaxE = max(allNeuronSTMApos2plot,[],'all');
        cMaxI = max(allNeuronSTMAneg2plot, [],'all');
        colorAxE = [0 cMaxE];
        colorAxI = [0 cMaxI];
    end
    
    subplot 121;
    imagesc(allNeuronSTMApos2plot); colormap(flipud(pink)); colorbar; caxis(colorAxE);
    title('Excitatory STAs');
    axes = gca; axes.XTick = [];    axes.YTick = [];
    subplot 122;
    imagesc(allNeuronSTMAneg2plot); colormap(flipud(pink)); colorbar; caxis(colorAxI)
    title('Inhibitory STAs')
    sgtitle([regions{array}.label ' All Cells STAs  : Session= ' num2str(sess2run) ],'FontSize',24);%' Null: ' params.null])
    axes = gca; axes.XTick = [];    axes.YTick = [];
    
    if saveFlag & params.cleanSTAflag & collapseArrays
        savefig(fig, [workDir 'allCellsSTMAs_' null '_clean_' regions{array}.label 'collapsed.fig']);
        saveas(fig, [workDir 'allCellsSTMAs_' null '_clean_' regions{array}.label 'collapsed.png']);
    elseif saveFlag & ~params.cleanSTAflag & ~collapseArrays
        savefig(fig, [workDir 'allCellsSTMAs_' null '_' regions{array}.label '.fig']);
        saveas(fig, [workDir 'allCellsSTMAs_' null '_' regions{array}.label '.png']);
    elseif saveFlag & params.cleanSTAflag & ~collapseArrays
        savefig(fig, [workDir 'allCellsSTMAs_' null '_clean_' regions{array}.label '.fig']);
        saveas(fig, [workDir 'allCellsSTMAs_' null '_clean_' regions{array}.label '.png']);
    elseif saveFlag & ~params.cleanSTAflag & collapseArrays
         savefig(fig, [workDir 'allCellsSTMAs_' null '_' regions{array}.label 'collapsed.fig']);
        saveas(fig, [workDir 'allCellsSTMAs_' null '_' regions{array}.label 'collapsed.png']);
    end
end

edges = 0:500:65000;

fig = figure;

for array = 1:length(regions)
    idx2use = ismember(allCellsSTMAacrossTime.unit(:,1), regions{array}.channels);
    
    allNeuronsSTMAposSize = squeeze(sum(allCellsSTMAacrossTime.pos(:,:,idx2use),[1,2]));
    allNeuronsSTMAnegSize =  squeeze(sum(allCellsSTMAacrossTime.neg(:,:,idx2use),[1,2])); 
    
    subplot 211; histogram(allNeuronsSTMAposSize,edges); hold on
    subplot 212; histogram(allNeuronsSTMAnegSize,edges); hold on
    sgtitle(['each cell STA size in Pixels, array=' regions{array}.label]);
    %pause;
    %disp(['max PosSTA size = ' num2str(max(allNeuronsSTMAposSize))])
    %disp(['max NegSTA size = ' num2str(max(allNeuronsSTMAnegSize))])
    
end
hold off

