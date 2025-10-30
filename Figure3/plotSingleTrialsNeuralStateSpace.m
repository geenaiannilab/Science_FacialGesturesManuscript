
%%%%%%%% PLOTTING Figure 3A
%%%%%%% Single-event population activity vectors (all cells from all regions, summed spike counts, -500 to +1000ms) 
%%%%%%% in the space of first three principal components, cluster by gesture-type 
%%%%%%% Input is nTrials x nCells , 
%%%%%%%   where each entry is # of spikes within a trial
%%%%%%%  Each point in state-space is a trial 
%%%%%%%
%%%%%%% data.allCellsPerTrialSpikeCount contains summed spike counts;  nTrials x nCells 
%%%%%%  data.allTrials contains which categorical gesture occured, nTrials
%%%%%% 
% written GRI 250911

clear all; close all;
set(0,'defaultAxesFontSize',36)
set(0,'defaultAxesFontWeight','bold')

data = load('matfiles/Fig3A.mat'); % change input file to examine plots from Supplementary Figures
thisRegion = 'All';

centerOnlyFlag = 0; % ONLY mean center input to PCA (dont divide by std)
centerNormalizeFlag = 1; % normalize cells/variables by their variance in addition to mean-centering (PCA on correlation matrix)    

allCellsPerTrialSpikeCount = data.allCellsPerTrialSpikeCount;
allTrials = data.allTrials;
bhvs2plot = data.bhvs2plot;
taxis = data.ax; 

% Single-Trial PCA 
% Input is nTrials x nCells , 
%    where each entry is # of spikes within a trial
%  Each point in state-space is a TRIAL 
%%
if centerOnlyFlag
     [coeff,score,latent, ~] = pca(allCellsPerTrialSpikeCount); % centering (subtract column/cell means) done automatically by pca()
     neuralInput = allCellsPerTrialSpikeCount - mean(allCellsPerTrialSpikeCount); %for plotting only, later
elseif centerNormalizeFlag
     [coeff,score,latent, ~] = pca(allCellsPerTrialSpikeCount,'VariableWeights', 'variance'); 
     neuralInput = (allCellsPerTrialSpikeCount - mean(allCellsPerTrialSpikeCount)) ./ std(allCellsPerTrialSpikeCount,[],1);
elseif sqrtFlag 
    [coeff,score,latent, ~] = pca(sqrt(allCellsPerTrialSpikeCount)); 
    neuralInput = sqrt(allCellsPerTrialSpikeCount);
else
    [coeff,score,latent, ~] = pca(allCellsPerTrialSpikeCount,'Centered','off');
    neuralInput = allCellsPerTrialSpikeCount;
    disp('performing PCA on uncentered, non-normalized data!')
end


% explained variance 
explained = latent/sum(latent);

% eigenvalues 
figure; 
plot(cumsum(latent/sum(latent)), 'ok','MarkerFaceColor','k') 
xlabel('PC'); ylabel('% variance'); xlim([0 80]);
title('Cumulative Variance Explained');

% plot perTrial data projected onto first 3 PCs 
%%
[~,~,id] = unique(allTrials.trialType);
markers = 'osd*.';
cols = [1, 0 0 ;   
        0, 0, 1;   
        0, 1, 0];      

blendFactor = 0.5;  % 0 = vivid, 1 = pure white
muted = cols * (1 - blendFactor) + blendFactor * ones(3,3);

plotFig = figure; 
plotFig.Position = [476 92 946 774];
hold on

bhvs2plot = unique(allTrials.trialType);  % or specify your desired order
for idx = 1:length(bhvs2plot)
    dat2plot = score(allTrials.trialType == bhvs2plot(idx), :);
    scat(idx) = scatter3(dat2plot(:,1), dat2plot(:,2), dat2plot(:,3), ...
        500, muted(idx,:), markers(idx), 'filled', ...
        'MarkerFaceAlpha', 0.5, 'MarkerEdgeColor', cols(idx,:), 'LineWidth', 2);
end

xlabel('Neural PC1')
ylabel('Neural PC2')
zlabel('Neural PC3')
title([thisRegion ': Single Trials of Facial Gesture']); hold off
legend({'Threat','Lipsmack','Chew'})
view(3); grid on; 

ax = gca;
set(ax, 'CameraPosition',  [-84.9252  76.9082  48.6309], ...
        'CameraTarget',    [2.5 2.5 0], ...
        'CameraUpVector',  [0 0 1], ...
        'CameraViewAngle', 10.18);
drawnow; 
set(ax, 'CameraPositionMode','auto', ...
        'CameraTargetMode','auto', ...
        'CameraUpVectorMode','auto', ...
        'CameraViewAngleMode','auto');

rotate3d(ax, 'on');     
% %plotFig = gcf;
%fileOut = ['/Users/geena/Dropbox/PhD/Manuscript/ScienceSubmission/ScienceRevision/ScienceResubmissionFigures/Fig3A_' thisRegion '.pdf'];
%exportgraphics(plotFig, fileOut, 'ContentType', 'vector');
