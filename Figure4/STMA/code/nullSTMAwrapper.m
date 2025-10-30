%%%%%%%%%%
%%%%%%%% PUPROSE: generates a null distribution of shuffled STAs, thru
%%%%%%%  using parforloop to engage multiple CPUs
%%%%%% MOST up to date version of this script & all of its functions live on the
%%%%%% server, under /Freiwald/gianni/STMA/Programs
%%%%%% Assumes you have already run 1) STMA2021wrapper.m, to crop & extract
%%%%%% video frames, and compute real STA per cell 
%%%%%% 
%%%%% DEPENDENCIES: shuffledSTMA.m
%%%%%%%%            generateShuffledSpikeTrain.m 
%%%%%%%%            generatPoissonSpikeTrain.m
%%%%%%%%            frameBounds4STMA2021.m
%%%%%%%%            getSpTrigAvg2021.m


channels2run = [ 8 1; 14 1; 17 1; 18 1; 21 1; 22 1; 24 1; 24 2; 25 1; 29 1; ...
   38 1; 38 2; 41 1; 41 2; 41 3; 46 1; 52 1; 53 1; 56 1; 56 2; 57 1; 60 1; 60 2; 61 1; ...
   66 1; 74 1; 75 1; 88 1; 96 1; ...
   97 1; 100 1; 101 1; 101 2; 101 3; 101 4; 102 1; 102 2; 104 1; 108 1; 108 2; 109 1; 112 1; 112 2; 113 1; 119 1; 119 2; 122 1; 127 1; 128 1; 128 2;...
129 2; 130 1; 134 1; 135 1; 136 1; 136 2; 140 1; 140 2; 145 1; 151 1; 151 2; 151 3; 152 1; 154 1; 157 1; 158 1; 159 1; 160 1; 160 2;...
    168 1; 171 1; 174 1; 177 1; 178 1; 178 2; 179 1; 184 1; 190 1; 190 2;...
    201 1; 201 2; 204 1; 204 2; 209 1; 210 1; 215 1; 217 1; 219 1; 219 2; 220 1; 221 1; 222 1;...
    225 1; 226 1; 226 2; 227 1; 228 1; 228 2; 233 1; 233 2; 238 1; 240 1];
sess2run = 5;

for jj = 1:length(channels2run)
    
clearvars -except jj channels2run sess2run

%%% setup 
date = '210704';
subject ='Barney';
ch = channels2run(jj,1);
unit=channels2run(jj,2);

% STMA variables 
bhvSR = 69.99;
windowSize = 1; % seconds to examine STMA 
directionality = 'post'; %pre/post
postSpikeFrames = ceil(bhvSR*windowSize);

% shuffled STMA variables
nCPUs = 40;
nIterations = 100;
shuffleWindow = 1; % in sec (NOT necessarily same as windowSize!)
plotSpShuffleFlag = 0; % plot the SDF of shuffled vs. real spikes 

% set up files to load 
workDir = (['/Freiwald/gianni/EphysBackup/' subject '_' date '/']);
vidDir = ['/Freiwald/gianni/STMA/' subject '_' date '/FT1/'];
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
vfiles = vfiles(logical(i2take),:);
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
TDTsessID = extractAfter(TDTfile, ['Barney-' date '-']);
FT1file = extractBefore(session2analyze{3},'.avi');
vidDirOut = [workDir subject '-' date '-' TDTsessID];

if ~exist(vidDirOut,'dir')
    mkdir(vidDirOut)
end

% load bhv data
TDTdir = [workDir '/Data4Analysis/'];
filePattern = fullfile(TDTdir, ['*' TDTsessID '*']);
bhvFile = dir(filePattern);
bhvFile = bhvFile.name;

obhv = load([TDTdir '/' bhvFile]);
obhv = obhv.obhv;
spikes = obhv.spikes.el(ch).sp(unit).times;

% load cropped/extracted frames 
filePattern = fullfile(vidDir,['cropped*' FT1file '*']);
videoFile = dir(filePattern);
videoFile = videoFile.name;

% camera variables different across subjects
if strcmpi(subject, 'Barney')
    faceCameraStream = 'FT1dat';
elseif strcmpi(subject, 'Thor')
    faceCameraStream = 'FT2dat';
end

%%% take absolute difference between frames as proxy for movement & reshape
%%% to 2D array for processing 

catFrames = load(videoFile);
catFrames = catFrames.catFrames;
pixelRows = size(catFrames,1);
pixelColumns = size(catFrames,2);

diffFrames = ProcessPixels(catFrames); 
clear catFrames;

diffFrames = reshape(diffFrames, [pixelRows*pixelColumns, size(diffFrames,3)]);

% % % NULL POPULATION CREATION 
% % % shuffle the spike times within a specified window, & re-compute the STMA x100
% % % 
% % % this portion of the script performs nIterations of a spike-triggered average analysis on shuffled spike data, 
% % % in order to
% % % create a null distribution against which to measure significance of the
% % % "real" spike-triggered average
% % % 
% % % we will use a parfor loop (parallel computing toolbox) to speed this up,
% % % by running up to 4 iterations simultaneously on 4 CPUs (my computer,
% % % yours may have more & use the max # the function detects in your machine)
% % % 
% % % OK SOME NOTES NOW: *takes mic*
% % % 0) diffIntensity is a "broadcast variable" -- it gets passed to every
% % % iteration of the function (see: it's an input to the parforloop below);
% % % this means it is copied to each CPU; if your broadcast variables are too
% % % large, you will run out of memory (ie, 2GB broadcast variable x 4 CPUs =
% % % 8GB RAM needed)
% % % 
% % % 1) by putting the parfor loop inside a function, and then adding that
% % % ',0', we force the desktop matlab to act as the first CPU, and we can set
% % % breakpoints/debug inside the parforloop (something that is not possible
% % % otherwise)
% % % 
% % % 2) each iteration acts on a separate CPU, they do not communicate, so
% % % there cannot be any variables inside the parforloop that depend on a
% % % previous iteration 
% % % 
% % % 3) indexing back into the structure we're filling with results
% % % ('tmpSTMA_fill') requries some creativity due to the limitations in how
% % % parfor loops allow indexing -- see "Sliced Variables" matlab help
% % % in general, you can only have ONE loop variable ('ii') 

tic 
parpool('local', nCPUs)

parfor ii=1:nIterations

      [tmpSTMA_fill(:,:,:,ii), simSpikeTimesAdj] = shuffledSTMA(ii, diffFrames, pixelRows, pixelColumns, spikes, obhv.(faceCameraStream).time, bhvSR, windowSize, directionality, shuffleWindow)


      fprintf('\bIteration\n');
      disp([num2str(ch) '-' num2str(unit)])

end

save([vidDirOut 'STMArand_ch' num2str(ch) 'unit' num2str(unit)],'ch','unit',...
    'directionality','FT1file','postSpikeFrames','shuffleWindow','TDTsessID','windowSize','-v7.3');
delete(gcp);

end