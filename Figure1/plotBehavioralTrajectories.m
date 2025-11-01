
%%%%%%%%
%%%%%%%% PLOTTING Figures 1 B-F; 
%%%%%%% Facial gestures are distinguishable during a naturalistic social paradigm
%%%%%%%%%
%%%%%%%%% GRI 09/11/2025 
%%%%%%%  markerTrajectories.mat contains results of behavioral trajectories
%%%%%%%         in face-marker space; contained in bhvTraj structure
%%%%%%%         (nTimepoints x nTrials x nPCs)
%%%%%%%  markerTrajectoriesTSNE.mat contains dimensionality-reduced
%%%%%%%         results (tsneResults.Y is nTimepoints x nDims) 

clear all; close all;
set(0,'defaultAxesFontSize',24)

thisDir = pwd;
workDir = [thisDir '/matfiles/'];
tsneResults = load([workDir '/markerTrajectoriesTSNE.mat']);
load([workDir '/markerTrajectories.mat']);
scalar = 3000; 
tmin = taxis(1); tmax = taxis(end);

%% plot 10% of trials starting position (THIS IS SLOW) 
figure;
colors = 'rbg';
sampleFrac = 0.02;
nTotalTrials = sum([size(markerData2plot(1).data,2),size(markerData2plot(2).data,2),size(markerData2plot(3).data,2)]);
randSample = ceil(sampleFrac*nTotalTrials);

for bhv = 1:length(bhvs2plot)
    
    thisData = markerData2plot(bhv).data;
    
    for tt = 1:size(thisData,2)
        for i = 1:length(markerLabels)
            xpos = (i*2) - 1;
            ypos = (i*2);
            x = markerData2plot(bhv).data(1,tt,xpos);
            y = -markerData2plot(bhv).data(1,tt,ypos);

            s = scatter(x,y,150,colors(bhv),'filled'); 
            s.MarkerFaceAlpha = 0.1; hold on; 

        end
    end
end
title('Starting Position for All Trials');
clear markerData2plot; 

%% tSNE embedding per trial 
allBhv = tsneResults.allBhv;
tsneResults = tsneResults.tsneResults;

figure;
h = gscatter(tsneResults.Y(:,1), tsneResults.Y(:,2),allBhv.trialTypeData,'rbg',[],12);
legend('Threat','Lipsmack','Chew'); legend boxoff
xlabel('t-SNE 1');
ylabel('t-SNE 2');
title('Characterization of Three Facial Expressions');

%% eigenspectrum
fig1 = figure; plot(cumsum(latent) ./sum(latent),'o','linew',8);
xlabel('PC'); ylabel('CumVar, %')
xlim([0 length(latent)])
title(['Marker Positions; Eigenspectrum; tmin = ' num2str(tmin) ' tmax = ' num2str(tmax)]);

%% 2D depiction of PC weights
for PC = 1:5
    fig = figure; 
    x=[]; y = []; radius = [];
    for i = 1:length(markerLabels)
        xpos = (i*2) - 1;
        ypos = (i*2);
        x(i) = meanPos(xpos);
        y(i) = -meanPos(ypos);
        radius(i) = PCweights(ypos,PC);
    end
    s = scatter(x,y,abs(radius)*scalar,radius,'filled'); 
    hold off; colormap(bluewhitered); axes = gca; axes.XTick = []; axes.XTickLabels = '';axes.YTick = []; axes.YTickLabels = '';
    title(['PC' num2str(PC) ' : ' num2str(round(explained(PC),0)) '% var'])
    cb = colorbar; cb.Label.String = 'PC weight';


end

%% trial-average & single-trial trajectories 
fig = figure('Position', get(0, 'Screensize')); 
set(gcf,'color','w');
colors = 'rbg';
markerType = 'osd*.'; markerSize = 300;

for bhv = 1:length(bhvTraj)

    data = (bhvTraj(bhv).data);

    for tt = 1:size(data,2)
        p = plot3(data(:,tt,1), data(:,tt,2),data(:,tt,3),...
            'Color',colors(bhv),'linew',0.5); hold on;
    end

end

for bhv = 1:length(bhvTraj)
    scatter3(bhvTraj(bhv).avg(:,1), bhvTraj(bhv).avg(:,2),bhvTraj(bhv).avg(:,3),...
        markerSize, colormap(customColormap(:,:,bhv)),'filled',markerType(bhv));  hold on     
end

hold off;
xlabel('Dim 1; au'); 
ylabel('Dim 2; au');
zlabel('Dim 3; au');
title(['Bhv Trajectories; tmin = ' num2str(tmin) ' tmax = ' num2str(tmax)])
grid on

view(3)
set(0,'defaultAxesFontSize',20);
set(0,'defaultAxesFontWeight','bold');

%% plot the PCs with mean/SEM
fig = figure('Position', get(0, 'Screensize')); 
for bhv = 1:length(bhvs2plot)
    for PC = 1:5
        subplot(1,5,PC)
        colors = 'rbmgc';

        shadedErrorBar(taxis,bhvTraj(bhv).avg(:,PC),bhvTraj(bhv).sem(:,PC),'lineprops',{'Color',colors(bhvs2plot(bhv)),'linew',8}); hold on;

        title(['PC=' num2str(PC)],'FontSize',20)
    end
end

legend('Thr','LS', 'Chew'); xlabel('Time, s'); ylabel('a.u.')
sgtitle(['FaceMotion PCs; tmin = ' num2str(tmin) ' tmax = ' num2str(tmax)],'FontSize',24)



