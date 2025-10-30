%%%%%%%%%
%%%%%%%%% PURPOSE: once you have run STMA2021wrapper & nullSTMAwrapper,
%%%%%%%%%   this will download each cell's real/shuffled STAs, convert them
%%%%%%%%%   to zscore, and then flatten across time, in order to evaluate if
%%%%%%%%%   cells of a given cortical region have spatially-distinct STAs. 
%%%%%%%%  OUTPUTS: allCellsSTMAAcrossTime_(null).mat -- is nPixelWidth x
%%%%%%%%                nPixelHeight x nCells array, where values are
%%%%%%%%                summed z-scores at that location, across
%%%%%%%%                time/latencies 
%%%%%%%%  
%%%%%%%%    
%%%%%%%%  DEPENDENCIES: real & shuffled STA per cell 

clear all; 
close all;
set(0,'defaultAxesFontSize',24)

date = '210704';
subject ='Barney';
sess2run = 2;

%% STMA variables 
if strcmp(subject,'Thor')
    bhvSR = 120;
else
    bhvSR = 69.99; 
end

cleanSTAflag = 1;
null = 'poiss'; % 'poiss','rand', or 'randLong' 

plotFlag = 1;
saveFlag = 0;

windowSize = 1; % seconds to examine STMA 
directionality = 'post'; %pre/post
zScoreThreshold = 3; % all pixel w/ absvalue < this, will be masked out 
thresholdPositive = 0; 
thresholdNegative = 0; 
thresholdBoth = 1;
nFrames = ceil(bhvSR*windowSize);  % to eval post-spike
spatialFilterSize = 3; % size of 2D gauss filter for smoothing

%% channels to examine 
channels2run =  [     38 1; 38 2; 41 1; 41 2; 
41 3; 46 1; 52 1; 53 1; 56 1; 56 2; 57 1; 60 1; 60 2; 61 1; ... 
   66 1; 74 1;75 1; 88 1; 96 1; ...
   97 1; 101 1; 101 2; 101 3; 101 4; 102 1; 102 2; 104 1; 108 1; 108 2; 109 1; 112 1; 112 2;113 1;
 119 1; 119 2; 122 1; 127 1; 128 1; 128 2;...
129 2; 130 1; 134 1; 135 1; 136 1;  136 2; 140 1; 140 2; 145 1; 151 1; 151 2; 151 3; 152 1; 154 1; 157 1; 158 1; 159 1; 160 1; 160 2;...
    168 1; 171 1; 174 1; 177 1; 178 1; 178 2; 179 1; 190 2;...
    201 1; 201 2; 204 1; 204 2; 209 1; 210 1; 215 1; 217 1; 219 1; 219 2; 220 1; 221 1; ...
    225 1; 226 1; 226 2; 227 1; 228 1; 228 2; 233 1; 233 2; 238 1; 240 1];

   %
%%  

% set up files to load 
bhvDir = (['/Users/geena/Dropbox/Thor_FTExperiment/SUAinfo/' subject '_' date '/Data4Analysis']);
STMAdir = ['/Users/geena/Dropbox/FaceExpressionAnalysis/STMA2021/' subject '_' date];
vidDirOut = STMAdir; 

% load bhv data
filePattern = fullfile(bhvDir, ['bhvspikes_sub' num2str(sess2run) '*']);
bhvFile = dir(filePattern);
TDTsessID = cell2mat(extractBetween(bhvFile.name,['sub' num2str(sess2run) '_'],['_']));
bhvFile = bhvFile.name;
workDir = [STMAdir '/' subject '-' date '-' TDTsessID '/'];

obhv = load([bhvDir '/' bhvFile]);
obhv = obhv.obhv;

% camera variables different across subjects
if strcmpi(subject, 'Barney')
    faceCameraStream = 'FT1dat';
elseif strcmpi(subject, 'Thor')
    faceCameraStream = 'FT2dat';
end

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
    
    % rotate by 90 if Thor's data
    if strcmp(subject,'Thor')
        rSTMA.STMA.pixels_mean = rot90(rSTMA.STMA.pixels_mean);
        rSTMA.STMA.pixels_std = rot90(rSTMA.STMA.pixels_std);
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

   if  strcmp(subject,'Thor')
        STMAzscoreFiltThresh = STMAzscoreFiltThresh(:,:,1:end-1);
    end
    

    % Remove any pixels *not* above zscore threshold in >=2
    % temporally consecutive frames in STMA 
    if cleanSTAflag
        [STMAzscoreFiltThresh] = cleanSTMA(STMAzscoreFiltThresh,0);
    end
    
    ispos = STMAzscoreFiltThresh .* (STMAzscoreFiltThresh > 0);
    isneg = STMAzscoreFiltThresh .* (STMAzscoreFiltThresh  < 0);
    
    allCellsSTMAacrossTime.pos(:,:,jj) = sum(ispos,3);
    allCellsSTMAacrossTime.neg(:,:,jj) =  sum(isneg,3);
    allCellsSTMAacrossTime.unit(jj,:) = [ch unit];
    
    if plotFlag
        figure('units','normalized','outerposition',[0 0 1 1]);
        
        cmax1 = max(sum(ispos,3),[],'all');
        cmin1 = min(sum(ispos,3),[],'all');
        cmax2 = max(sum(isneg,3),[],'all');
        cmin2 = min(sum(isneg,3),[],'all');

        subplot 121; imagesc(allCellsSTMAacrossTime.pos(:,:,jj)); title('Excitatory STA'); colormap(pink); %caxis([cmin1 cmax2]); 
        axes = gca; axes.XTick = []; axes.YTick = [];

        subplot 122; imagesc(allCellsSTMAacrossTime.neg(:,:,jj)*-1); colorbar;  title('Inhibitory STA'); colormap(pink); %caxis([cmin1 cmax2]); 
        axes = gca; axes.XTick = []; axes.YTick = [];

        sgtitle(['Cell: B-' num2str(ch) '-' num2str(unit) ],'FontSize',28)
        pause; clf
    end
    
end

%%%%%%%%% FIX THIS ! 

if saveFlag 
    params.date = date;
    params.subject = subject; 
    params.sess2run = sess2run;
    params.cleanSTAflag = cleanSTAflag;
    params.null = null; % 'poiss','rand', or 'randLong' 
    params.windowSize = windowSize; % seconds to examine STMA 
    params.directionality = directionality; %pre/post
    params.zScoreThreshold = zScoreThreshold; % all pixel w/ absvalue < this, will be masked out 
    params.thresholdPositive = thresholdPositive; 
    params.thresholdNegative = thresholdNegative; 
    params.thresholdBoth = thresholdBoth;
    params.nFrames = nFrames;  % to eval post-spike
    params.spatialFilterSize = spatialFilterSize; % size of 2D gauss filter for smoothing
end

if saveFlag & cleanSTAflag
    save([workDir 'allCellsSTMAacrossTime_' null '_clean_highThres.mat'],'allCellsSTMAacrossTime','params', '-v7.3')
elseif saveFlag & ~cleanSTAflag
    save([workDir 'allCellsSTMAacrossTime_' null '.mat'],'allCellsSTMAacrossTime','params', '-v7.3')
end



