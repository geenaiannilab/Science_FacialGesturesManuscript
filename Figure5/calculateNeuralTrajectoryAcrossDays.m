% written GI updated 220104
%
% PUPRPOSE: Calculate reduced neural "state space"
%   where each point in reduced neural state-space represents a BIN in TIME

% https://www.sciencedirect.com/science/article/pii/S0896627313002237#bib21
% https://elifesciences.org/articles/54313.pdf
% https://www2.tntech.edu/leap/murdock/books/v1chap3.pdf

close all;
clear all;
set(0,'defaultAxesFontSize',28)

subject2analyze ={'combined'};
nIterations = 50;
popType = 'balanced';
saveResultsFlag = false; 
plotResultsFlag = true; 
workdir = (['/Users/geena/Dropbox/PhD/SUAinfo/Pseudopopulations/' popType '/N_50cells/']);

% DO NOT CHANGE!
win = 0.02; % in sec
tmin = 0.8;
tmax = 0.8;
taxis = -tmin:win:tmax;  

centerOnlyFlag = 0; % ONLY mean center input to PCA (dont divide by std)
centerNormalizeFlag = 1; % normalize cells/variables by their variance in addition to mean-centering (PCA on correlation matrix)    
threshlowFRFlag = 1; % 0/1;  remove low FR neurons 
threshlowFR = 0.1; %  sp/s threshold 
%%% NOTE RE: PREPROCESSING INPUT
%%% pca() in matlab does mean-centering *by default*, but does not divide by stdev by default 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate & then downsample gaussKern for convolution w/ spikes binned at "win" resolution
smParams.smoothSize = 50; %std in milisec
smParams.filterType = 1; % 1= gauss, 2= caus half-gauss, 3=box
smWin = genSmWin(smParams);
gaussKern = dwnsampleKernel(smWin,win);

extraBins = length(gaussKern); % extra bins to discard due to conv edges 
tmin = tmin +(extraBins*win); 
tmax = tmax + (extraBins*win);
edges = -tmin:win:tmax;  % use 'edges' to rebin (before conv)

% get rid of the edges where convolution results in artifacts & only use
% those bins later 
tmin2take = tmin - (extraBins*win); 
tmax2take = tmax - (extraBins*win);
taxis2take = -tmin2take:win:tmax2take;
taxis2takeIdx = edges >= taxis2take(1) & edges <= taxis2take(end);
        
% 
originalTaxis = -1:0.001:1.2;

%%%%%%%%%%%%%%%%%%
for ss = 1:length(subject2analyze)

    thisSubj = subject2analyze{ss};

    % repeat for each pseudopopulation
    for iter = 1:nIterations
        
        % load the pseudopopulation
        if strcmpi(subject2analyze{ss},'combined')  
            pseudoPop = load([workdir 'twoSubjPseudopopulation_N' num2str(iter) '.mat']);
        else
            pseudoPop = load([workdir 'pseudopopulation_N' num2str(iter) '.mat']);

        end

        % define some variables
        arrayList = pseudoPop.regionList;
        bhvs2plot = unique(pseudoPop.pseudoPopulation.(subject2analyze{:}).bhvTrials);
        trialLabels = pseudoPop.pseudoPopulation.(thisSubj).bhvTrials;

        for aa = 1:length(arrayList)

            nCells = size(pseudoPop.pseudoPopulation.(thisSubj).(arrayList{aa}),3);

            % dimensions are ntrials x nbins x ncells
            allCellsPSTH.(arrayList{aa}) = zeros(length(trialLabels),length(edges)-1,nCells);

            for trial = 1:size(allCellsPSTH.(arrayList{aa}),1) %per trial
                for cell = 1:size(allCellsPSTH.(arrayList{aa}),3) % per cell
                        
                    raster = pseudoPop.pseudoPopulation.(thisSubj).(arrayList{aa})(trial,:,cell);
                    rasterTimes = originalTaxis(find(raster));
                    newBinnedSpikes = histcounts(rasterTimes, edges); % bin w/ extras for discarding conv artifact

                    allCellsPSTH.(arrayList{aa})(trial, :, cell) = conv(newBinnedSpikes, gaussKern, 'same') ./win; %divide by win for FR
                end
            end
            
             % remove low FR neurons
            if threshlowFRFlag
                meanFRs = squeeze(mean(allCellsPSTH.(arrayList{aa}),[1 2]));
                allCellsPSTH.(arrayList{aa}) = allCellsPSTH.(arrayList{aa})(:,:,meanFRs' >= threshlowFR);
            end

            %  average FRs; this is the trial-averaged response PER cell
            for i = 1:length(bhvs2plot)
                theseTrials = (pseudoPop.pseudoPopulation.(thisSubj).bhvTrials == bhvs2plot(i));
                allCellsTrialAvgResponse(i).bhv = squeeze(mean(allCellsPSTH.(arrayList{aa})(theseTrials, taxis2takeIdx, :)));
                allCellsTrialAvgResponse(i).nReps = sum(pseudoPop.pseudoPopulation.(thisSubj).bhvTrials == bhvs2plot(i));
            end

            % concatente all averaged responses for PCA input 
            neuralInput = vertcat(allCellsTrialAvgResponse(:).bhv);
            nFeatures = size(neuralInput,2);

            % Trial-Averaged PCA 
            % Input is nBins x nCells
            if centerOnlyFlag
                [coeff,score,latent, ~] = pca(neuralInput); % centering (subtract column/cell means) done automatically by pca()
                input2plot = neuralInput - mean(neuralInput); %for plotting only, later
            elseif centerNormalizeFlag
                [coeff,score,latent, ~] = pca(neuralInput,'VariableWeights', 'variance');
                input2plot = (neuralInput - mean(neuralInput)) ./ std(neuralInput,[],1);
            else
                [coeff,score,latent, ~] = pca(neuralInput,'Centered','off');
                input2plot = neuralInput;
                disp('performing PCA on uncentered, non-normalized data!')
            end


            % explained variance
            explained = latent/sum(latent);

            % pull out per bhv neural Traj & calc velocity 
            start = 1;
            for bhv = 1:length(bhvs2plot)
                stop = length(taxis2take)*bhv;
                neuralTraj.(arrayList{aa}).bhv(bhv).score{iter} = score(start:stop,:);
                start = start + length(taxis2take);
            end

            clearvars -except neuralTraj pseudoPop popType subject2analyze thisSubj nIterations workdir...
            win tmin tmax tmin2take tmax2take taxis taxis2takeIdx taxis2take edges originalTaxis minRestFlag minRest ...
            centerOnlyFlag centerNormalizeFlag threshlowFRFlag threshlowFR...
            gaussKern ss iter aa arrayList bhvs2plot trialLabels plotResultsFlag saveResultsFlag;


        end % end arrayList 

    end % end iterations 
end % end subjects 


%
% calculate & plot how far away in neural space starting positions of each trajectory are 

% get rid of extra dimensions 
for aa = 1:length(arrayList)
    for bhv = 1:length(bhvs2plot)
        minDim.(arrayList{aa}).dim(bhv) = min( cellfun(@(c) size(c,2), neuralTraj.(arrayList{aa}).bhv(bhv).score));

    end
end

for aa = 1:length(arrayList)
    for bhv = 1:length(bhvs2plot)
        thisMinDim = minDim.(arrayList{aa}).dim(bhv);
        
        for iter = 1:nIterations 
            neuralTraj.(arrayList{aa}).bhv(bhv).score{1,iter} = neuralTraj.(arrayList{aa}).bhv(bhv).score{1,iter}(:,1:thisMinDim);
        
        end

    end
end

%%
for aa = 1:length(arrayList)
    
    nDims = min(minDim.(arrayList{aa}).dim,[],'all');
    maxDim = nDims;

        for iter = 1:nIterations

                % position coordinates 
                ThrCoords = neuralTraj.(arrayList{aa}).bhv(1).score{1,iter};
                LSCoords = neuralTraj.(arrayList{aa}).bhv(2).score{1,iter};
                ChewCoords = neuralTraj.(arrayList{aa}).bhv(3).score{1,iter};
    
                %  velocity 
                velocity.(arrayList{aa}).Thr = (diff(ThrCoords)) ./ (taxis2take(end) - taxis2take(1));
                velocity.(arrayList{aa}).LS = (diff(LSCoords)) ./ (taxis2take(end) - taxis2take(1));
                velocity.(arrayList{aa}).Chew = (diff(ChewCoords)) ./ (taxis2take(end) - taxis2take(1));
    
               % velocity as a function of time (velocity being composed of
               % 1 to nMaxDimensions in neural space)

                for dd = 1:maxDim

                    % nNeuralDimensional distance vector, per bhv
                    sosThrChDistance = 0;
                    sosThrLSDistance = 0;
                    sosLSChDistance = 0 ;

                    % nNeuralDimensional velocity vector, per bhv
                    sosThrVel = 0;
                    sosLSVel = 0;
                    sosChewVel = 0;
                    
                    for dim = 1:dd
                        
                        % nNeuralDimensional distance vector, per bhv

                        sosThrChDistance = sosThrChDistance + ((ThrCoords(:,dim) - ChewCoords(:,dim)) .^2);
                        sosThrLSDistance = sosThrLSDistance + ((ThrCoords(:,dim) - LSCoords(:,dim)) .^2);
                        sosLSChDistance = sosLSChDistance + ((LSCoords(:,dim) - ChewCoords(:,dim)) .^2);

                        % nNeuralDimensional velocity vector, per bhv
                        sosThrVel = sosThrVel + (velocity.(arrayList{aa}).Thr(:,dim) .^2);
                        sosLSVel = sosLSVel + (velocity.(arrayList{aa}).LS(:,dim) .^2);
                        sosChewVel = sosChewVel + (velocity.(arrayList{aa}).Chew(:,dim) .^2);

                    end

                    % nNeuralDimensional distance vector, per bhv
                    distanceThrVCh.(arrayList{aa}).data(:,dd,iter) = sqrt(sosThrChDistance);
                    distanceThrVLS.(arrayList{aa}).data(:,dd,iter) = sqrt(sosThrLSDistance);
                    distanceLSVCh.(arrayList{aa}).data(:,dd,iter) = sqrt(sosLSChDistance);

                    % nNeuralDimensional velocity vector, per bhv
                    velThrTotal.(arrayList{aa}).data(:,dd,iter) = sqrt(sosThrVel);
                    velLSTotal.(arrayList{aa}).data(:,dd,iter) = sqrt(sosLSVel);
                    velChewTotal.(arrayList{aa}).data(:,dd,iter) = sqrt(sosChewVel);

                    
                end % end dimensionality 

        end % end iterations 

        % view the results
        if plotResultsFlag
            for dd = [1 3 12 20]%1:maxDim
                
                fig1 = figure;
                for ii = 1:nIterations
                    plot(taxis2take,distanceThrVCh.(arrayList{aa}).data(:,dd,ii),'Color',[0.3010 0.7450 0.9330]); hold on;
                    plot(taxis2take,distanceLSVCh.(arrayList{aa}).data(:,dd,ii),'Color',[0.4660 0.6740 0.1880])
                    plot(taxis2take,distanceThrVLS.(arrayList{aa}).data(:,dd,ii),'Color',[0.6350 0.0780 0.1840]	)
                end
                title([arrayList{aa} ' Distance as fxn of Dim=' num2str(dd)])
                legend('Thr-Ch','LS-Ch','Thr-LS')
                %saveas(fig1, ['~/Desktop/' arrayList{aa} '_trajDist_dim' num2str(dd) '.png'])

                fig2 = figure; 
                for ii = 1:nIterations
                    plot(taxis2take(1:end-1),velThrTotal.(arrayList{aa}).data(:,dd,ii),'r'); hold on;
                    plot(taxis2take(1:end-1),velLSTotal.(arrayList{aa}).data(:,dd,ii),'b')
                    plot(taxis2take(1:end-1),velChewTotal.(arrayList{aa}).data(:,dd,ii),'g')
                end
                title([arrayList{aa} ' velocity as fxn of Dim=' num2str(dd)])
                %saveas(fig2, ['~/Desktop/' arrayList{aa} '_trajSpeed_dim' num2str(dd) '.png'])
                

            end
        end

end % end arrayList 

if saveResultsFlag
    save([workdir 'neuralTrajectory_' subject2analyze{:} '.mat'],'velChewTotal','velLSTotal','velThrTotal',...
        'distanceThrVCh','distanceLSVCh','distanceThrVLS','subject2analyze','nIterations','popType','win','tmin2take', 'tmax2take');
    
end

    
function binnedSpikes = binSpikes(evtime, sptimes, win, tmin, tmax)
%%%%%%%%%
ax = -tmin:win:tmax;
binnedSpikes = zeros(length(evtime),length(ax)-1); %initialize

for i=1:length(evtime)
    curT = evtime(i);
    spt = sptimes(find(sptimes >curT-tmin & sptimes < curT + tmax))-curT;%find spikes times after tmin & before tmax
                                                                         % subtracting curT converts spiketimes to -tmin:tmax scale 
                                                                         % (ie, first possible spktime is -tmin)
                                                                          
    binnedSpikes(i,:) = histcounts(spt, ax); %returns integer of # of spikes occurring in that bin


end
end
