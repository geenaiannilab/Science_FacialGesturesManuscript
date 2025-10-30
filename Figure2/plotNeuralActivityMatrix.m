
%%%%%%%%
%%%%%%%% PLOTTING Figure 2B 
%%%%%%%%  Single-Cell Activity and Selectivity in Cortical Face-Motor Regions
%%%%%%%%
%%%%%%%% Sample population activity of all simultaneously recorded cells 
%%%%%%%% (top, sorted by area) and facial marker positions (bottom) 
%%%%%%%%% GRI 09/11/2025 

close all; clear all;
set(0,'defaultAxesFontSize',24);
set(0,'defaultAxesFontWeight','bold');

load('matfiles/Fig2B.mat');

%%%%%
%%%%% PLOTTING
fig = figure('units','normalized','outerposition',[0 0 1 1]);
fig.Position =[ 0.3690    0.0855    0.6138    0.7953];

h1 = axes(fig,'position',[0.13 0.3 0.75, 0.6]);
h2 = axes(fig,'position',[0.13 0.1 0.75, 0.19]);

% raster
axes(h1);
imagesc(newTimeAxis, 1:size(raster,1), raster); colormap(h1,(1-gray)); cb = colorbar; ylabel(cb,'Spikes/second','FontSize',18,'Rotation',270)
ylabel('Cell ID, Sorted (top:S1, bottom, M3')
h1.YLabel; h1.XTick = []; h1.YTick = []; h1.YTickLabel = [];
cb.Position = [0.89 0.3 0.0106    0.6001];
cb.Label.Units = 'normalized';cb.Label.Position = [6.7 0.5 0 ];
caxis([0 80])

% face marker traces
axes(h2)
plot(newTimeAxis,DLCOut(:,50:76),'linew',2); colormap(h2, hsv);
ylabel('Marker Position','FontSize',18); xlabel('Time,s')
h2.YAxisLocation = 'left';
h2.YTick = [];

linkaxes([h1 h2],'x')
xlim([330 360])
axes(h1);title(blockname)
%subtitle([obhv.Info.sessType])
fig; xlim([330 360]);

