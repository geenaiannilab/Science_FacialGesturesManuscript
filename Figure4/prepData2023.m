
%%% 
%%% PURPOSE: PSID 2023; prepping neural & behavioral trial data 
%%% Will bin spikes & behavioral (DLC face marker) timesseries at
%%% same resolution, into trials 

%%% return bhvTrialData binned from -tmin:analysisWin:tmax +/- smoothed
%%% returns neuralTrialData binned from -tmin:analysisWin:tmax, +/- smoothed
%%%    ***note, you are NOT smoothing bhv in this incantation (because
%%%    taking PCs does effectively the same, and DLC already smoothed it)
%%%
%%% Pools all trials across all days, all subjects 
%%% Saves output (as input to makePseudopopPSID2023.m) 
%%%
%%% OUTPUTS:[region]_allTrials.mat, which contains:
% allBhvTrials is nTrials x 1 cell array; each contains nTimepoints x
% nBhvPCs from a single trial
%
% allNeuralTrials is nTrials x 1 cell array; each contains nTimepoints x
% nCells from a single trial
%
% allLabels is nTrials x 3 cell array ;
% 1st column is nCells x1 cell array of channel #s
% 2nd column is nCells x1 cell array of regionLabels
% 3rd column is subject name
%
% cellListOut is pooledList of cellIDs (similar to allLabels, 1st column,
% but explicitly separated by date & subject) 

clear all; 
subjects = {'Barney','Thor'};
dates = {'Barney_210704', 'Barney_210706','Thor_171010','Thor_171027'};
regionList = {'S1','M1','PMv','M3','All'};

saveFlag = true; 

bhvs2train = [ 1 ];
bhvs2test = [1 ];
outputTrials = 'ThrOnly'; % for naming the output files 
markers2use =  'all'; excludeTongue = 1; 
markerPCA = true; 
bhvPCs = 1:10; velocityFlag = false; 
tmin = 0.5; tmax = 1.5; 

newBinSize = 0.01; % 0.01; in s; for binning neural data initially
neuralSmoothSize = 50; % for neural data  
analysisWin = 0.05; % common timestep for behavior and neural data 

% pool all data across days/sessions
counter = 1;
for ss = 1:length(subjects)
    theseDates = dates(contains(dates,subjects{ss}));
 
    for dd = 1:length(theseDates)
        workDir = ['/Users/geena/Dropbox/PhD/SUAinfo/' theseDates{dd} '/Data4Analysis'];

        % get relevant files
        flist = natsortfiles(dir([workDir '/bhvspikes_sub*.mat']));
        flist = fullfile({flist.folder}, {flist.name});

        for ff = 1:length(flist)
            obhvIn = load(flist{ff}); obhvIn = obhvIn.obhv;
            DLCdata(counter).subject = extractBefore(theseDates{dd},'_'); 
            DLCdata(counter).date = extractAfter(theseDates{dd},'_');

            DLCdata(counter).Info = obhvIn.Info;
            DLCdata(counter).FT1dat = obhvIn.FT1dat;
            DLCdata(counter).FT2dat = obhvIn.FT2dat;
            DLCdata(counter).evScore = obhvIn.evScore;
            DLCdata(counter).spikes = obhvIn.spikes;
            counter = counter + 1;
        end
        clear obhvIn;
    end % end theseDates
end % end subjects 

%%%%% Prep bhv PCs & spike data 

for rr = 1:length(regionList)
    allNeuralTrials = [];
    allNeuralDat = [];
    allBhvTrials = [];
    allBhvLabels = [];
    allLabels = cell(1,1);

    startC = 1;
    stopC = 0;

    for dd = 1:length(DLCdata)

        subject = DLCdata(dd).subject;
        thisDate = DLCdata(dd).date;

        [regionMap] = getChannel2CorticalRegionMapping(subject, 1);
        chls = regionMap{1,rr}.channels;

        % remove pesky channels
        if strcmpi(subject, 'Barney') && strcmpi(thisDate, '210704') 
            if strcmpi(regionList{rr},'PMv') 
                chls = chls(chls~= 129);
            elseif strcmpi(regionList{rr},'M1') 
                chls = chls(chls~= 46);
            elseif strcmpi(regionList{rr},'M3') 
                chls = chls(chls~= 240);
            elseif strcmpi(regionList{rr},'All')
                chls([46 129 240]) = [];
            end
        elseif strcmpi(subject, 'Barney') && strcmpi(thisDate, '210706') 
            if strcmpi(regionList{rr},'M1') 
                chls = chls(chls~= 86);
            elseif strcmpi(regionList{rr},'All')
                chls = chls(chls~= 86);
            end
        end

        % camera variables different across subjects
        if strcmpi(subject, 'Barney') || strcmpi(subject, 'Dexter')
            faceCameraStream = 'FT1dat'; cameraSync = DLCdata(dd).(faceCameraStream).time; bhvSRorig = 70;
        elseif strcmpi(subject, 'Thor')
            faceCameraStream = 'FT2dat'; cameraSync = DLCdata(dd).(faceCameraStream).time; bhvSRorig = 120;
        end

        pulseLength = length(DLCdata(dd).(faceCameraStream).time);
        nFramesExported = size(DLCdata(dd).(faceCameraStream).Markers,1);
        if nFramesExported > pulseLength
            disp(['Warning -- sessionID ' num2str(DLCdata(dd).Info.blockname) ' dropped frames'])
            disp(['nFramesExported from DLC = ' num2str(nFramesExported)]);
            disp(['pulse Length = ' num2str(pulseLength)]);
        end

        %%%%%%% PREPARE the neural & behavioral data for PSID
        %%%%%
        if isempty(DLCdata(dd).evScore(:,5))
            continue
        end

        % process DLC data, chunk into trials, &
        % take only trials preceded by stillness & within bounds of camera pulse
        origTimeBase = cameraSync; DLCwin = mean(diff(origTimeBase));
        [DLCmarkers, markerLabels] = processDLCMarkers('DLCin', DLCdata(dd).(faceCameraStream), 'markers2use', markers2use, 'excludeTongue',true);

        DLCdata(dd) = spacedBhvsPSID(DLCdata(dd),tmin,1);
        [onsets, evindx] = identifyTrials2Take2023(DLCdata(dd), nFramesExported, bhvs2train, origTimeBase, tmin, tmax);
        if isempty(onsets)
            continue
        end
        % return bhvTrialData binned from -tmin:analysisWin:tmax +/- smoothed
        % at END of timebin (shifted one forward relative to neural) 
        [bhvTrials] = prepBhvTrialData2023('DLCdata', DLCdata(dd), 'bhvMeasure', DLCmarkers, 'markerPCA', markerPCA, 'bhvPCs', bhvPCs, ...
            'onsets', onsets, 'evindx', evindx, ...
            'tmin', tmin, 'tmax', tmax, 'win', DLCwin, 'commonBinSize', analysisWin, 'smoothFlag', false);

        % returns neuralTrialData binned from -tmin:analysisWin:tmax, +/- smoothed
        [neuralTrials,neuralDat] = prepNeuralTrialData2023('DLCdata', DLCdata(dd),'chls', chls, ...
            'onsets', onsets, 'evindx', evindx, ...
            'newBinSize', newBinSize, 'tmin', tmin, 'tmax', tmax, 'commonBinSize', analysisWin,'smoothFlag', true, 'smoothSize', neuralSmoothSize);

        % pool all trials across days, subjects
        allNeuralTrials = cat(1,allNeuralTrials,neuralTrials);
        allNeuralDat = cat(1,allNeuralDat, neuralDat);
        allBhvTrials = cat(1,allBhvTrials, bhvTrials.bhvData);
        allBhvLabels = cat(1,allBhvLabels, bhvTrials.trialType);

        cellLabels = {neuralDat.cellID}; cellLabels = repmat(cellLabels,1,length(neuralTrials))';
        arrayLabels = {neuralDat.array}; arraylabels = repmat(arrayLabels,1,length(neuralTrials))';
        subjectLabels =  {DLCdata(dd).subject}; subjectLabels = repmat(subjectLabels,1,length(neuralTrials))';
        dateLabels = {DLCdata(dd).date}; dateLabels = repmat(dateLabels,1,length(neuralTrials))';

        stopC = stopC + length(neuralTrials);

        allLabels(startC:stopC,1) = cellLabels;
        allLabels(startC:stopC,2) = arrayLabels;
        allLabels(startC:stopC,3) = subjectLabels;
        allLabels(startC:stopC,4) = dateLabels;

        startC = startC + length(neuralTrials);
    end % end subsessions 

    % pool single list of cellID labels, per subject, across days
    cellListOut = cell(1,length(subjects));

    for ss = 1:length(subjects)

        subject = subjects{ss};
        dates = unique(allLabels(strcmpi(allLabels(:,3),subject),4),'stable');

        for dd = 1:length(dates)
            datesIndx = strcmpi(dates{dd},allLabels(:,4));
            dayCellLabels = allLabels(datesIndx,1); dayCellLabels = vertcat(dayCellLabels{:});
            cellListOut{1,ss}{1,dd} = unique(vertcat(dayCellLabels{:}),'stable');
            cellListOut{1,ss}{2,dd} = repmat(dates(dd),length(cellListOut{1,ss}{1,dd}),1);
        end

        cellListOut{2,ss}{1,1} = subject;

    end


    if saveFlag
        save(['/Users/geena/Dropbox/PhD/SUAinfo/PSIDResults2023/' regionList{rr} '_' outputTrials 'Trials.mat'],...
            'allBhvTrials', 'allBhvLabels','allLabels','allNeuralTrials', 'cellListOut',...
            'analysisWin','tmin','tmax','bhvs2train','bhvs2test','markers2use','markerPCA','bhvPCs');
    end

    clear cellListOut
end

