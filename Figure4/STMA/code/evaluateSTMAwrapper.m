%%%
%%% written GI 220201
%%% PURPOSE: once you have computed a cell's STMA with STMAwrapper2021.m, &
%%% computed it's null population shuffles via nullSTMAwrapper.m, this
%%% script will for each cell, download real & shuffled simulated STAs,
%%% transform to z-score, plot & save a movie of each cell's STA

%%% setup 
date = '210704';
subject ='Barney';
sess2run = 2;
channels2run =  [     38 1];
%[ 2 1; 9 1; 11 1; 22 1; 24 1; 32 1; ...
%    33 1; 34 1; 34 2; 36 1; 36 2; 38 1; 38 2; 39 1; 39 2; 40 1; 41 1; 41 2; 42 1; 42 2; 42 3; 43 1; 44 1; 48 1; 50 1 ; 52 2; 52 3; 53 1; 53 2; 54 1; 55 1;]
% 56 1; 57 1; ...
%     59 1; 59 2; 59 3; 60 1; 60 2; 61 1 ;  61 2; 62 1; 62 2; 64 1; 64 2; ...
%     66 1; 66 2; 74 1; 74 2; 86 1; 86 2; 91 1; 91 2; 94 1; ...
%     98 1; 100 1; 101 1; 102 1; 102 2; 103 1; 103 2; 104 1; 105 1; 105 2; 106 1; 108 1; 108 2; 109 1; 109 2; 110 1; 111 1; 112 1; 116 1; 123 1; 124 1; 125 1;...
%     133 1; 137 1; 137 2; 138 1; 141 1; 141 2; 142 1; 142 2; 145 1; 145 2; 153 1; 157 1; 159 1; ...
    
% STMA variables 
bhvSR = 69.99;
windowSize = 1; % seconds to examine STMA 
directionality = 'post'; %pre/post

zScoreThreshold = 3; % NOTE: IF YOU WANT TO EVAL NEG STMAS, THIS NEEDS TO BE NEGATIVE
thresholdPositive = 0; 
thresholdNegative = 0; 
thresholdBoth = 1;
nFrames = ceil(bhvSR*windowSize);  % to eval post-spike
spatialFilterSize = 3; % size of 2D gauss filter for smoothing


% set up files to load 
bhvDir = (['/Users/geena/Dropbox/PhD/SUAinfo/' subject '_' date '/Data4Analysis']);
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
    sSTMAfile = ([workDir 'STMApoiss_ch' num2str(ch) 'unit' num2str(unit) '.mat']);

    if exist(rSTMAfile,'file')
        rSTMA = load(rSTMAfile);
    else
        disp([rSTMAfile '--> rSTMA doesnt exist. Quitting']);
        return
    end

    if exist(sSTMAfile,'file')
        sSTMA = load(sSTMAfile);
    else
        disp([sSTMAfile '--> rSTMA doesnt exist. Quitting']);
        return
    end

    % rotate by 90 if Thor's data
    if strcmp(subject,'Thor')
        rSTMA.STMA.pixels_mean = rot90(rSTMA.STMA.pixels_mean);
        rSTMA.STMA.pixels_std = rot90(rSTMA.STMA.pixels_std);
    end
    
    %convert rSTMA to zScore
    nSpikes = rSTMA.STMA.nSpikes;
    meanSim = sSTMA.sSTMA.meanSim;
    stdSim = sSTMA.sSTMA.stdSim;
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
    % temporally consecutive frames in STMA 
    %[cleanSTMAzscoreFiltThresh] = cleanSTMA(STMAzscoreFiltThresh,1);

    if  strcmp(subject,'Thor')
        STMAzscoreFiltThresh = STMAzscoreFiltThresh(:,:,1:end-1);
    end
    
    % prep movie 
    loops = size(STMAzscoreFiltThresh,3);
    F(loops) = struct('cdata',[],'colormap',[]);
    figure('units','normalized','outerposition',[0 0 1 1]);
    cmax = max(STMAzscoreFiltThresh(find(STMAzscoreFiltThresh)));
    cmin = min(STMAzscoreFiltThresh(find(STMAzscoreFiltThresh)));
    if isempty(cmax)
        cmax = 3;
    end
    if isempty(cmin)
        cmin = -3;
    end
    if cmin == cmax 
        cmin = cmin - 1;
    end
    
    for i = 1:loops
        ax(1) = subplot(1,2,1);
        ax(1).XTick = [];
        ax(1).YTick = [];
        imagesc(STMAzscoreFiltThresh(:,:,i)) ; caxis([cmin cmax]); colorbar
        ax(2) = subplot(1,2,2);
        ax(2).XTick = [];
        ax(2).YTick = [];
        imagesc(STMAzscoreFilt(:,:,i));        caxis([cmin cmax]);
        tt = i/bhvSR; tt = round(tt*100)/100;
        sgtitle(['ch' num2str(ch) '-' num2str(unit) ' ~' num2str(tt) '  s post-spike ; avg ' num2str(nSpikes) ' spikes; Thres =' num2str(zScoreThreshold)]);
        colorbar
        colormap(flipud(pink));

        F(i) = getframe(gcf);    

    end

    %write video out to STMA results directory
    myVideo = VideoWriter([vidDirOut '\' subject '-' date '-' TDTsessID '\STMAPoissmap_ch' num2str(ch) 'unit' num2str(unit) '.avi']);
    myVideo.FrameRate = 30;  % reduced from default of 30 
    open(myVideo);
    writeVideo(myVideo, F);
    close(myVideo)

    close all;

end

