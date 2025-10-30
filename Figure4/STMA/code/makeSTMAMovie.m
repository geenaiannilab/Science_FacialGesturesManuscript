%%%
%%% written GI 230815
%%% PURPOSE: for specified cell & subsession, this script will download 
%%% real & shuffled simulated STAs,
%%% transform to z-score, plot & save a movie (and stills) of STMA

%%% setup 
date = '210704';
subject ='Barney';
sess2run = 2;
channels2run =  [     38 1];
   
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
    [STMAzscoreFiltThresh] = cleanSTMA(STMAzscoreFiltThresh,0);

    if  strcmp(subject,'Thor')
        STMAzscoreFiltThresh = STMAzscoreFiltThresh(:,:,1:end-1);
    end
    
    % prep movie 
    loops = size(STMAzscoreFiltThresh,3);
    F(loops) = struct('cdata',[],'colormap',[]);
    figure('outerposition',[476    74   851   792]);
    cmax = 10; %max(STMAzscoreFiltThresh(find(STMAzscoreFiltThresh)));
    cmin = 2; %min(STMAzscoreFiltThresh(find(STMAzscoreFiltThresh)));
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
        
        imagesc(STMAzscoreFiltThresh(:,:,i));
        ax(1) = gca;
        ax(1).XTick = [];
        ax(1).YTick = [];
        caxis([cmin cmax]); colorbar
        
        tt = i/bhvSR; tt = round(tt*100)/100;
        title(['Neuron #' num2str(ch) '-' num2str(unit) ' ~' num2str(tt) '  s post-spike']);
        colormap(flipud(pink));

        F(i) = getframe(gcf);    

    end

    %write video out to STMA results directory
    %myVideo = VideoWriter([vidDirOut '\' subject '-' date '-' TDTsessID '\STMAPoissmap_ch' num2str(ch) 'unit' num2str(unit) '.avi']);
    myVideo = VideoWriter([ subject '-' date '-' TDTsessID '\STMAPoissmap_ch' num2str(ch) 'unit' num2str(unit)], 'MPEG-4');
    myVideo.FrameRate = 10;  % reduced from default of 30 
    open(myVideo);
    writeVideo(myVideo, F);
    close(myVideo)

    close all;

end

