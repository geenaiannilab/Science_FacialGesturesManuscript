%%%%%%%%%
%%%%%%%% PURPOSE: compute spike-triggered movement average for each cell
%%%%%%%% DEPENDENCIES: ExtractMovieFrames2021.m
%%%%%%%%               frameBounds4STMA2021.m
%%%%%%%%               getSpTrigAvg2021.m 

date = '210704';
subject ='Barney';
sess2run = 6;
nFramesExported = 57997;
channels2run = [62 2; 64 1; 64 2; ...
    66 1; 66 2; 74 1; 74 2; 86 1; 86 2; 91 1; 91 2; 94 1; ...
    98 1; 100 1; 101 1; 102 1; 102 2; 103 1; 103 2; 104 1; 105 1; 105 2; 106 1; 108 1; 108 2; 109 1; 109 2; 110 1; 111 1; 112 1; 116 1; 123 1; 124 1; 125 1;...
    133 1; 137 1; 137 2; 138 1; 141 1; 141 2; 142 1; 142 2; 145 1; 145 2; 153 1; 157 1; 159 1; 160 1;...
    168 1; 179 1; 183 1; 183 2; 185 1; 186 1; 190 1; 191 1; 191 2; 191 3];%...%[33 1; 34 1; 34 2; 36 1; 36 2; 38 1; 38 2; 39 1; 39 2; 40 1; 41 1; 41 2; 42 1; 42 2; 42 3; 43 1; 44 1; 48 1; 50 1  52 2; 52 3; 53 1; 53 2; 54 1; 55 1; 55 2; 56 1; 57 1; ...
%    59 1; 59 2; 59 3; 60 1; 60 2; 61 1 ;  61 2; 62 1; ; 

% behavioral video variables
extractFramesFlag = 0;
downsampleScale = 1;
cropOption = 'applyCrop';

% STMA variables
bhvSR = 120;%69.99;
windowSize = 1; % seconds to examine STMA 
directionality = 'post'; %pre/post
postSpikeFrames = ceil(bhvSR*windowSize);

% set up files to load 
workDir = (['K:/' subject '_' date '/']);
if strcmpi(subject,'Barney') || strcmpi(subject,'Dexter')
    vidDir = [workDir '/FT1/'];
elseif strcmpi(subject,'Thor') 
    vidDir = [workDir '/FT2/'];
end

infoDir = ([workDir 'TDTdata/Info']);

% info & video files 
infoname = [infoDir '/' date '_param.mat'];
if exist(infoname,'file')
    load(infoname);
else
    disp([infoname '--> info file doesnt exist. Quitting']);
    return
end

videoFiles = [workDir '/' subject '_' date '_fnames_relation.mat'];
if exist(videoFiles,'file')
    vfiles = load(videoFiles);
    vfiles = vfiles.files;
else
    disp([videoFiles '--> video fileRelation doesnt exist. Quitting']);
    return
end

i2take = ones(size(infodat.SubsetTable,1),1);
i2take(strcmpi(infodat.SubsetTable(:,2),'test'))=0;
i2take(strcmpi(infodat.SubsetTable(:,2),'noise'))=0;
%vfiles = vfiles(logical(i2take),:);
sessions = infodat.SubsetTable(logical(i2take),:);

for i = 1:length(sessions)
    if ~(strcmpi(vfiles{i,1},sessions{i,1}))
        disp(['vfile= ' vfiles{i,1} ' and tdtfile= ' sessions{i,1} ' do not match; check!'])
        return
    end
end

%disp(sessions);
%input = input('indicate which subsession # to run STMA shuffle ---');
%session2analyze = vfiles(input,:);
session2analyze = vfiles(sess2run,:);

TDTfile = session2analyze{1};
TDTsessID = extractAfter(TDTfile, [subject '-' date '-']);

if strcmpi(subject,'Barney') || strcmpi(subject,'Dexter')
    FT1file = extractBefore(session2analyze{3},'.avi');
elseif strcmpi(subject,'Thor')
    FT1file = extractBefore(session2analyze{2},'.avi');
end
% load bhv data
TDTdir = [workDir '/Data4Analysis/'];
filePattern = fullfile(TDTdir, ['*' TDTsessID '*']);
bhvFile = dir(filePattern);
bhvFile = bhvFile.name;
obhv = load([TDTdir '/' bhvFile]);
obhv = obhv.obhv;

% load cropped/extracted frames 
filePattern = fullfile(vidDir,['cropped_' FT1file '*']);
videoFile = dir(filePattern);
videoFile = videoFile.name;

% camera variables different across subjects
if strcmpi(subject, 'Barney') || strcmpi(subject, 'Dexter')
    faceCameraStream = 'FT1dat';
    cameraSync = obhv.(faceCameraStream).time;
elseif strcmpi(subject, 'Thor')
    faceCameraStream = 'FT2dat';
    cameraSync = obhv.(faceCameraStream).time(1:nFramesExported);
end

%%% apply cropping & extract movie frames into .mat file, eg. ('catFrames.mat')
if extractFramesFlag
    ExtractMovieFrames2021([vidDir FT1file],[vidDir videoFile],cropOption,downsampleScale);
end
%

%%% take absolute difference between frames as proxy for movement & reshape
%%% to 2D array for processing 

catFrames = load([vidDir videoFile]);
catFrames = catFrames.catFrames;
pixelRows = size(catFrames,1);
pixelColumns = size(catFrames,2);

diffFrames = ProcessPixels(catFrames); 
clear catFrames;


for jj = 1:length(channels2run)
   
    clear STMA spikeTimesAdj framesLowerBound framesUpperBound i pixelRow pixelColumn currPixel avg_r std_r nSpikes

            ch = channels2run(jj,1); %channel & unit # of neuron to examine
            unit = channels2run(jj,2); %channel & unit # of neuron to examine
            disp(channels2run(jj,:));

            % preallocate
            STMA = struct('pixels_mean', [], 'pixels_std', [],'nSpikes', []);
            STMA.pixels_mean = zeros(pixelRows,pixelColumns,postSpikeFrames+1);
            STMA.pixels_std = zeros(pixelRows,pixelColumns,postSpikeFrames+1);

            % return frame #s that correspond to windowSize prior to each spiketime
                % new inputs = spikeTimes, cameraSync
            [spikeTimesAdj, framesLowerBound, framesUpperBound] = frameBounds4STMA2021(obhv.spikes.el(ch).sp(unit).times, cameraSync, ...
                bhvSR, directionality, windowSize);

            % for each pixel in the video, generate "average stimulus" vector & store
             tic 
             for i = 1:pixelColumns*pixelRows
                 
                    [pixelRow, pixelColumn] = ind2sub([pixelRows pixelColumns], i);
                    currPixel = squeeze(diffFrames(pixelRow,pixelColumn,:))'; 
                    [avg_r, std_r, nSpikes] = getSpTrigAvg2021(currPixel, spikeTimesAdj, framesLowerBound, framesUpperBound, directionality, windowSize, postSpikeFrames);

                    STMA.pixels_mean(pixelRow, pixelColumn, :) = avg_r;
                    STMA.pixels_std(pixelRow,pixelColumn,:) = std_r;

             end
             toc

            % store relevant info in output structure 
            STMA.nSpikes = nSpikes;
            STMA.spikeTimes = spikeTimesAdj';
            STMA.ch = ch;
            STMA.unit = unit;
            STMA.directionality = directionality;

            % save it
            outdir = ['D:/Dropbox/FaceExpressionAnalysis/STMA2021/' subject '_' date '/' subject '-' date '-' TDTsessID '/' ];
            if ~exist(outdir, 'dir')
               mkdir(outdir)
            end

            save([outdir 'STMA_ch' num2str(ch) 'unit' num2str(unit)],'STMA','-v7.3');
end





