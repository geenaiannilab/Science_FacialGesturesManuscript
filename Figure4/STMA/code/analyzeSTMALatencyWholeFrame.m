

%%%%%%%%%
%%%%%%%%% PURPOSE: analyze potential differences in latencies (timepoints
%%%%%%%%%   of max influence, of cell spiking --> movement) across cells 
%%%%%%%%  Per cell, will download real & shuffled STA, convert to zscore,
%%%%%%%%  threshold (or not), separate into negative/positive STAs, and then calculate VARIANCE as well as MEAN FRAME,
%%%%%%%%        across *entire frame* at each timepoint post-spike. 
%%%%%%%   Will end up with allCellsSTMALatency.MeanperFrame(pos/neg) &
%%%%%%%   allCellsSTMALatency.VarLatency(pos/neg), which are both of dimensions nCells x
%%%%%%%         nTimepoints post-spike
%%%%%%%%  OUTPUTS: will generate allCellsVarlatency_(null).mat, which is
%%%%%%%%  per cell variance, & per cell mean frame, for all post-spike
%%%%%%%%  timebins
%%%%%%%%
%%%%%%%%
%%%%%%%%  DEPENDENCIES: real & shuffled STA per cell 
clear all; 
close all;

date = '210704';
subject ='Barney';
sess2run = 5;
null = 'poiss'; % 'poiss','rand', or 'randLong' 
collapseArraysFlag = 0;    %  plot all M1s, all PMvs, and all M3s together, or each array separately 
plotFlag  = 1;
saveFlag = 1;

channels2run =   [ 14 1; 17 1; 18 1; 21 1; 22 1; 24 1; 24 2; 25 1; 29 1; ...
38 1; 38 2; 41 1; 41 2; 41 3; 46 1; 52 1; 53 1; 56 1; 56 2; 57 1; 60 1; 60 2; 61 1; ... 
   66 1; 74 1; 75 1; 88 1; 96 1; ...
   97 1; 100 1; 101 1; 101 2; 101 3; 101 4; 102 1; 102 2; 104 1; 108 1; 108 2; 109 1; 112 1; 112 2; 113 1; 119 1; 119 2; 122 1; 127 1; 128 1; 128 2;...
129 2; 130 1; 134 1; 135 1; 136 1; 136 2; 140 1; 140 2; 145 1; 151 1; 151 2; 151 3; 152 1; 154 1; 157 1; 158 1; 159 1; 160 1; 160 2;...
    168 1; 171 1; 174 1; 177 1; 178 1; 178 2; 179 1; 190 2;...
    201 1; 201 2; 204 1; 204 2; 209 1; 210 1; 215 1; 217 1; 219 1; 219 2; 220 1; 221 1; 222 1;...
    225 1; 226 1; 226 2; 227 1; 228 1; 228 2; 233 1; 233 2; 238 1; 240 1];

% STMA variables 
bhvSR = 69.99;
windowSize = 1; % seconds to examine STMA 
directionality = 'post'; %pre/post

zScoreThreshold = 0; % 
thresholdPositive = 0; 
thresholdNegative = 0; 
thresholdBoth = 1;
nFrames = ceil(bhvSR*windowSize);  % to eval post-spike
spatialFilterSize = 3; % size of 2D gauss filter for smoothing

% set up files to load 
bhvDir = (['D:\Dropbox\Thor_FTExperiment\SUAinfo\' subject '_' date '/Data4Analysis']);
STMAdir = ['D:\Dropbox\FaceExpressionAnalysis\STMA2021\' subject '_' date];
vidDirOut = STMAdir; 

% load bhv data
filePattern = fullfile(bhvDir, ['bhvspikes_sub' num2str(sess2run) '*']);
bhvFile = dir(filePattern);
TDTsessID = cell2mat(extractBetween(bhvFile.name,['sub' num2str(sess2run) '_'],['_']));
bhvFile = bhvFile.name;
workDir = [STMAdir '\' subject '-' date '-' TDTsessID '\'];

obhv = load([bhvDir '/' bhvFile]);
obhv = obhv.obhv;

% camera variables different across subjects
if strcmpi(subject, 'Barney')
    faceCameraStream = 'FT1dat';
elseif strcmpi(subject, 'Thor')
    faceCameraStream = 'FT2dat';
end

%mouthCrop = load([workDir '\mouthROI.mat']);

for jj = 1:length(channels2run)

    ch = channels2run(jj,1);
    unit=channels2run(jj,2);

    rSTMAfile = ([workDir 'STMA_ch' num2str(ch) 'unit' num2str(unit) '.mat']);
    sSTMAfile = ([workDir 'STMA' null '_ch' num2str(ch) 'unit' num2str(unit) '.mat']);

    if exist(rSTMAfile,'file')
        rSTMA = load(rSTMAfile,'STMA');
    else
        disp([rSTMAfile '--> rSTMA doesnt exist. Quitting']);
        continue
    end

    if exist(sSTMAfile,'file')
        sSTMA = load(sSTMAfile);
    else
        disp([sSTMAfile '--> rSTMA doesnt exist. Quitting']);
        continue
    end

    %convert rSTMA to zScore
    nSpikes = rSTMA.STMA.nSpikes;
    try
        meanSim = sSTMA.sSTMA.meanSim;
        stdSim = sSTMA.sSTMA.stdSim;
    catch ME
         meanSim = sSTMA.rSTMA.meanSim;
        stdSim = sSTMA.rSTMA.stdSim;
    end
    STMAzscore = (rSTMA.STMA.pixels_mean - meanSim) ./ stdSim; 

    
    %spatially smooth with 1x1 gaussian, threshold 
    STMAzscoreFilt = imgaussfilt(STMAzscore,spatialFilterSize, 'FilterDomain', 'spatial');
    
    %threshold data 
    if thresholdPositive
        STMAzscoreFiltThresh = STMAzscoreFilt >= zScoreThreshold;
    elseif thresholdNegative 
        STMAzscoreFiltThresh = STMAzscoreFilt <= zScoreThreshold;
    elseif thresholdBoth
        STMAzscoreFiltThresh = abs(STMAzscoreFilt) >= zScoreThreshold;
    end

    % mask out non-threshold crossing pixels
    STMAzscoreFiltThresh = (STMAzscoreFilt .* STMAzscoreFiltThresh);

    % Remove any pixels *not* above zscore threshold in >=2

    ispos = STMAzscoreFiltThresh .* (STMAzscoreFiltThresh > 0);
    isneg = STMAzscoreFiltThresh .* (STMAzscoreFiltThresh  < 0);
    
    STMApos = reshape(ispos, [size(ispos,1)*size(ispos,2), size(ispos,3)]);
    STMAneg = reshape(isneg, [size(isneg,1)*size(isneg,2), size(isneg,3)]);
     
    allCellsSTMALatency.MeanperFrame.pos(jj,:) = squeeze(mean(ispos,[1,2]))';
    allCellsSTMALatency.MeanperFrame.neg(jj,:) = squeeze(mean(isneg,[1,2]))';
    allCellsSTMALatency.Varlatency.pos(jj,:) = var(STMApos);
    allCellsSTMALatency.Varlatency.neg(jj,:) =  var(STMAneg);
    allCellsSTMALatency.unit(jj,:) = [ch unit];
    
end

if saveFlag
    save([workDir 'allCellsVarlatency_' null '.mat'],'allCellsSTMALatency','null','zScoreThreshold','-v7.3');
end

if plotFlag
    
    [regions] = getChannel2CorticalRegionMapping(subject, collapseArraysFlag);
    
    timeAx = (0:nFrames) ./ bhvSR;

    for array = 1:length(regions)
        fig(1) = figure; 

        idx2use = ismember(allCellsSTMALatency.unit(:,1), regions{array}.channels);

        thisArrayPos = squeeze(allCellsSTMALatency.MeanperFrame.pos(idx2use,:)');
        thisArrayNeg = squeeze(allCellsSTMALatency.MeanperFrame.neg(idx2use,:)');

        thisArrayPosVar = allCellsSTMALatency.Varlatency.pos(idx2use,:)';
        thisArrayNegVar = allCellsSTMALatency.Varlatency.neg(idx2use,:)';

        subplot 221; 
        plot(timeAx, thisArrayPos); xlabel('cell'); ylabel('meanFrame'); 
        title('Mean Frame value, Positive STAs'); 
        subplot 222;
        plot(timeAx, abs(thisArrayNeg)); xlabel('cell'); ylabel('meanFrame'); 
        title('Mean Frame value, Negative STAs')

        subplot 223;
        plot(timeAx, thisArrayPosVar);xlabel('frame post-spike'); ylabel('var');
        title('variance Frame value, Positive STAs')

        subplot 224;
        plot(timeAx, thisArrayNegVar);xlabel('frame post-spike'); ylabel('var');
        title('variance Frame value, Negative STAs')

        sgtitle(['Latencies of Cells in  ' regions{array}.label ' : session= ' num2str(sess2run)])
       % savefig(fig, [workDir '\' regions{array}.label '_latency_' null '.fig'])
    end

end