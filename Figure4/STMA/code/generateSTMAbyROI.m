%%%%%%%%
%%%%%%%%% PURPOSE: for each cell, compute zscored negative/positive STAs, and
%%%%%%%%  then active a GUI, prompting you to define facial ROIs -- ( 7 for now, 'mouth','Lcheek','Rcheek','nose','eyes','Lear','Rear');
%%%%%%%%
%%%%%%%%  The program will then save both an ROI facial template (faceROI_template.mat), as well as a .mat file per ROI 
%%%%%%%%     which gives the coordinates of a rectangle which defines the ROI
%%%%%%%%
%%%%%%%%  The program will finally pull out the zscored STA for each cell,
%%%%%%%%  *only within* the given ROIs, so that statistics can be done per
%%%%%%%%  cell, per ROI 
%%%%%%%%
%%%%%%%% OUTPUTS: Saves 'allCellsSTMAbyROI_(null).mat', which is nCells x 1
%%%%%%%%    structure containing 14 ROIs per cell (7 computed on positive STA,
%%%%%%%%    7 on negative; each ROI is nPixels x nPixels x nTimepoints
%%%%%%%%    post-spike
%%%%%%%%
%%%%%%%% DEPENDENCIES: real & shuffled STAs
%%%%%%%%               defineFaceROIs.m
%%%%%%%%               applyFaceROIs.m

clear all; 
close all;

date = '171010';
subject ='Thor';
sess2run = 6;

defineROIsFlag = 0; % if not already defined 
saveFlag = 1;

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
    33 1; 34 1; 34 2; 36 1; 36 2; 38 1; 38 2; 39 1; 39 2; 40 1; 41 1; 41 2; 42 1; 42 2; 42 3; 43 1; 44 1; 48 1; 50 1 ; 52 2; 52 3; 53 1; 53 2; 54 1; 55 1; 55 2; 56 1; 57 1; ...
    59 1; 59 2; 59 3; 60 1; 60 2; 61 1 ;  61 2; 62 1; 62 2; 64 1; 64 2; ...
    66 1; 66 2; 74 1; 74 2; 86 1; 86 2; 91 1; 91 2; 94 1; ...
    98 1; 100 1; 101 1; 102 1; 102 2; 103 1; 103 2; 104 1; 105 1; 105 2; 106 1; 108 1; 108 2; 109 1; 109 2; 110 1; 111 1; 112 1; 116 1; 123 1; 124 1; 125 1;...
    133 1; 137 1; 137 2; 138 1; 141 1; 141 2; 142 1; 142 2; 145 1; 145 2; 153 1; 157 1; 159 1; 160 1;...
    168 1; 179 1; 183 1; 183 2; 185 1; 186 1; 190 1; 191 1; 191 2; 191 3];
%%

% set up files to load 
bhvDir = (['D:/Dropbox/Thor_FTExperiment/SUAinfo/' subject '_' date '/Data4Analysis']);
STMAdir = ['D:/Dropbox/FaceExpressionAnalysis/STMA2021/' subject '_' date];
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

% define ROIs by hand with GUI if necessary 
if defineROIsFlag
    defineFaceROIs('date',date, 'subject', subject, 'TDTsessID', TDTsessID);
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
    
    % Separate excitatory & inhibitory effects in STMA
    STMApos = STMAzscoreFiltThresh .* (STMAzscoreFiltThresh > 0);
    STMAneg = STMAzscoreFiltThresh .* (STMAzscoreFiltThresh  < 0);
    
    ROIsSTMA(jj).pos = applyFaceROIs('date',date, 'subject', subject, 'STMA', STMApos, 'TDTsessID', TDTsessID);
    ROIsSTMA(jj).neg = applyFaceROIs('date',date, 'subject', subject, 'STMA', STMAneg, 'TDTsessID', TDTsessID);
    ROIsSTMA(jj).unit = [ch unit];
    
end


if saveFlag

    if strcmp(subject,'Thor')
        params.bhvSR = 120;
    else
        params.bhvSR = 69.99; 
    end
    params.null = null;
    params.cleanSTAflag = cleanSTAflag;
    params.windowSize = windowSize; % seconds to examine STMA 
    params.directionality = directionality; %pre/post
    params.zScoreThreshold = zScoreThreshold; % 
    params.thresholdPositive = thresholdPositive; 
    params.thresholdNegative = thresholdNegative; 
    params.thresholdBoth = thresholdBoth;
    params.nFrames =nFrames;  % to eval post-spike
    params.spatialFilterSize = spatialFilterSize; % size of 2D gauss filter for smoothing

    params.ROIs = ROIs;
    
    save([workDir 'allCellsSTMAbyROI_' null '.mat'],'ROIsSTMA','params','-v7.3');
end