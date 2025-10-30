%%%%% Analyzes data for  Figures 1 B-F; 
%%%%%%% Facial gestures are distinguishable during a naturalistic social paradigm
%%%%%  Generates trial-averaged, and individual trial behavioral trajectories from DLC markers
%%%%%%
%%%%%%%  markerTrajectories.mat contains results of behavioral trajectories
%%%%%%%         in face-marker space; contained in bhvTraj structure
%%%%%%%         (nTimepoints x nTrials x nPCs)
%%%%%%%  markerTrajectoriesTSNE.mat contains dimensionality-reduced
%%%%%%%         results (tsneResults.Y is nTimepoints x nDims) 
%%%%% 
%%%%%%%%%
%%%%%%%%% GRI 09/11/2025 

close all; clear all;
set(0,'defaultAxesFontSize',24);
set(0,'defaultAxesFontWeight','bold');

%% data params 
subject = 'Barney';
date = '210704';
runTSNE = false;
subsessions2plot = [2:5 7:10];
bhvs2plot = [1 2 4];
gestureNames = {'threat','lipsmack','chew'};
markers2use =  'all'; excludeTongue = 1; plotDLCFlag = false; 
newBinSize = 0.01;
tmin = 1;
tmax = 2;
centerNormalizeFlag = 0;
centerFlag = 1;
smoothSize = 5; 

minRestFlag = 1;
minRest = abs(tmin); % in sec; minimal rest prior to move onset (trials to include)

%% set up files to load 

workDir = ['/Users/geena/Dropbox/PhD/SUAinfo/' subject '_' date '/Data4Analysis'];
infoDir = ([workDir 'TDTdata/Info']);

%% concatenate subsession files 
% get relevant files 
flist = natsortfiles(dir([workDir '/bhvspikes_sub*.mat']));
flist = fullfile({flist.folder}, {flist.name}); 

% determine which type of subsession to plot
if strcmpi(subsessions2plot,'face expression')
    subs2use = flist(contains(flist, 'face expression'));
elseif strcmpi(subsessions2plot,'visual expression')
    subs2use = flist(contains(flist, 'visual expression'));
elseif strcmpi(subsessions2plot,'chew expression')
    subs2use = flist(contains(flist, 'chew expression'));
elseif strcmpi(subsessions2plot,'chew')
    subs2use = flist(contains(flist, 'chew'));
elseif strcmpi(subsessions2plot,'rest')
    subs2use = flist(contains(flist, 'rest'));
elseif strcmpi(subsessions2plot, 'all')
    subs2use = flist;
else 
    for i = 1:length(subsessions2plot)
       subs2use(i) = flist(contains(flist, ['sub' num2str(subsessions2plot(i)) '_']));
    end
end

for session = 1:length(subs2use)    
    obhvIn = load(subs2use{session});
    obhvIn = obhvIn.obhv;
    obhv(session) = obhvIn;
end
clear obhvIn;

subsessionInds = [];
allBhv.markers = [];
allBhv.trialTypeData= [];
allBhv.s = []; 
nTotalTrials = [];

for ss = 1:length(subs2use)

    sessionID = obhv(ss).Info.blockname;

    % camera variables different across subjects
    if strcmpi(subject, 'Barney') || strcmpi(subject, 'Dexter')
        faceCameraStream = 'FT1dat';
        bhvSR = 70;

    elseif strcmpi(subject, 'Thor')
        faceCameraStream = 'FT2dat';
        bhvSR = 120;
    end
    
    pulseLength = length(obhv(ss).(faceCameraStream).time);
    nFramesExported = size(obhv(ss).(faceCameraStream).Markers,1);
    if nFramesExported > pulseLength
        disp(['Warning -- sessionID ' num2str(obhv(ss).Info.blockname) ' dropped frames'])
        disp(['nFramesExported from DLC = ' num2str(nFramesExported)]);
        disp(['pulse Length = ' num2str(pulseLength)]);
    end

    %%%%%%% extract DLC marker positions & labels from obhv
    %%%%%
    % 
    DLCin = obhv(ss).(faceCameraStream);
    DLCtimeBase = DLCin.time(1:nFramesExported);
    [DLCmarkers, markerLabels] = processDLCMarkers('DLCin', DLCin, 'markers2use', markers2use,'excludeTongue', true,'smoothSize',smoothSize,'plotFlag',plotDLCFlag);
    subsessionInds(ss+1) = size(DLCmarkers,1); % length of each subsession

    %%%% extract trials 
    obhvThis = takeSpacedBhvs(obhv(ss),minRest,minRestFlag); % remove trials w/ <500 ms rest prior move onset
    allTrials = ismember(cell2mat(obhvThis.evScore(:,5)),bhvs2plot);

    % get rid of chew trials inside faceExp sets
    if strcmpi(subject, 'Thor')  && cell2mat(obhvThis.evScore(1,6)) == 3
        bhv2lose = ~(cell2mat(obhvThis.evScore(:,5)) == 4); 
        allTrials = logical(allTrials.*bhv2lose);
    end

    % remove events whose trial duration is outside of framesExported
    fmin = tmin*bhvSR; % convert to frames
    fmax = tmax*bhvSR;

    allOnsets = cell2mat(obhvThis.evScore(:,1));

    if strcmpi(sessionID, 'Thor-171010-154538') 
        disp(['sessionID = ' sessionID '  -- adjusting onsets for video export']);
        allOnsets = allOnsets - 2000; 
        %theseOnsets = theseOnsets((theseOnsets-fmin) >= 0);
    end

    theseTrials = (allOnsets >= (fmin )) & (allOnsets + fmax <= nFramesExported);
    theseTrials = logical(allTrials .* theseTrials);
        
    trialTypes = cell2mat(obhvThis.evScore(theseTrials,5));
    nTotalTrials(ss,:) = hist(trialTypes,bhvs2plot);

    theseOnsets = allOnsets(theseTrials); % handscored index, in subsession

    % extract trialData
    for tt = 1:length(theseOnsets)
        onsetMin = theseOnsets(tt) - fmin;
        onsetMax = theseOnsets(tt) + fmax;
        thisData = DLCmarkers(onsetMin:onsetMax,:); % DLCmarkers is nFrames x (nMarkersx2)

        %%%% downsample the DLC tracking data
        DLCtimeBase = -tmin:(1/bhvSR):tmax;
        [thisData, downsampledTimeAxis] = rebinBhv(thisData, DLCtimeBase, newBinSize);

        thisTrialTypeData = repmat(trialTypes(tt), [size(thisData,1)],1);

        allBhv.markers = [allBhv.markers; thisData];
        allBhv.trialTypeData = [allBhv.trialTypeData; thisTrialTypeData];
        allBhv.s = [allBhv.s ; obhv(ss).Info.blockname];
        
    end
    
   

end

%%% perform decomposition 
%% Input is nTimepoints x nMarkers,
if centerFlag
     [coeff,score,latent, ~, explained] = pca(allBhv.markers); % centering (subtract column/cell means) done automatically by pca()
elseif centerNormalizeFlag
     [coeff,score,latent, ~, explained] = pca(allBhv.markers,'VariableWeights', 'variance'); 
else
    [coeff,score,latent, ~, explained] = pca(allBhv.markers,'Centered','off');
    disp('performing PCA on uncentered, non-normalized data!')
end


%% per expression type, per trial, extract X&Y positions from -tmin/tmax around onset 
%  bhvTraj(bhv).data = nTimepoints x nMarkers x nTrials 
for bhv = 1:length(bhvs2plot)

    theseTrials = (allBhv.trialTypeData) == bhvs2plot(bhv);
    trialDur = length(downsampledTimeAxis);
    nTrials = sum(nTotalTrials);

    bhvTraj(bhv).data = reshape(score(theseTrials,:),[trialDur, nTrials(bhv),size(DLCmarkers,2)]);
    
    bhvTraj(bhv).avg = squeeze(mean(bhvTraj(bhv).data,2));
    bhvTraj(bhv).sem = squeeze(std(bhvTraj(bhv).data,[],2) ./ sqrt(nTrials(bhv)));

end

%%% do the TSNE calculations
%%
if runTSNE

    nPCs2use = find(cumsum(latent) ./sum(latent) >= 0.95,1,'first');
    data = score(:,1:nPCs2use);
    tsneResults = tsne_reproducible(data, 123, 'NumDimensions', 2, 'Perplexity', 30);
   
    % Save result for later
    figure;
    h = gscatter(tsneResults.Y(:,1), tsneResults.Y(:,2),allBhv.trialTypeData,'rbg',[],12);
    legend('Threat','Lipsmack','Chew'); legend boxoff
    xlabel('t-SNE 1');
    ylabel('t-SNE 2');
    title('Characterization of Three Facial Expressions');
end
%%
%% per expression type, per trial, extract new projection of markers 
%  bhvTraj(bhv).data = nTimepoints x nMarkers x nTrials 
for bhv = 1:length(bhvs2plot)

    theseTrials = (allBhv.trialTypeData) == bhvs2plot(bhv);
    trialDur = length(downsampledTimeAxis);
    nTrials = sum(nTotalTrials);

    bhvTraj(bhv).data = reshape(score(theseTrials,:),[trialDur, nTrials(bhv),size(DLCmarkers,2)]);
    
    bhvTraj(bhv).avg = squeeze(mean(bhvTraj(bhv).data,2));
    bhvTraj(bhv).sem = squeeze(std(bhvTraj(bhv).data,[],2) ./ sqrt(nTrials(bhv)));

end

% original marker data, separated by bhv 
for bhv = 1:length(bhvs2plot)

    theseTrials = (allBhv.trialTypeData) == bhvs2plot(bhv);
    trialDur = length(downsampledTimeAxis);
    nTrials = sum(nTotalTrials);

    markerData2plot(bhv).data = reshape(allBhv.markers(theseTrials,:),[trialDur, nTrials(bhv),size(DLCmarkers,2)]);
    
    markerData2plot(bhv).avg = squeeze(mean(markerData2plot(bhv).data,2));
    markerData2plot(bhv).sem = squeeze((std(markerData2plot(bhv).data,[],2)) ./ sqrt(nTrials(bhv)));

end

meanPos = mean(allBhv.markers,1);
taxis = downsampledTimeAxis;
PCweights = coeff(:,1:10);
percentVar = explained; 
[customColormap] = continuousColorMap(length(downsampledTimeAxis));

clear DLCin; clear obhv; clear DLCmarkers;

%% save tSNE embedding 
%    save('tsne_run123.mat','R');

%% save files 
outDir = ['/Users/geena/Dropbox/PhD/SUAinfo/' subject '_' date '/DLC'];
save([outDir '/markerTrajectories.mat'], ...
    'bhvTraj','PCweights','meanPos', 'coeff','latent','explained','taxis','bhvs2plot','gestureNames','markerData2plot','markerLabels','customColormap');

%%
function result = tsne_reproducible(X, seed, varargin)
%TSNE_REPRODUCIBLE Run t-SNE with a fixed RNG seed and save state
%
% result = tsne_reproducible(X, seed, ...)
%   X       : data matrix [N x D]
%   seed    : integer seed for reproducibility

    % Save RNG state before running
    prevRng = rng;             % store to restore later if you like
    rng(seed);                 % set deterministic seed
    
    rngState = rng;            % record exact RNG state used
    Y = tsne(X, varargin{:});  % run tsne with user options
    
    % Pack results
    result.Y    = Y;
    result.rng  = rngState;
    result.seed = seed;
    
    % Restore previous RNG state so this call doesn’t affect your workflow
    rng(prevRng);
end


function [bhvOut,downsampledTimeBaseVector] = rebinBhv(bhvIn, timeBaseOld, newBinsize)
    
downsampledTimeBaseVector = timeBaseOld(1):newBinsize:timeBaseOld(end);
bhvOut = zeros(length(downsampledTimeBaseVector), size(bhvIn,2));

for j = 1:size(bhvIn,2)
    bhvOut(:,j) = interp1(timeBaseOld, bhvIn(:,j), downsampledTimeBaseVector);
end

end

