
%%%%%% PLOTTING FIGURE 5 (left and middle panels) 
%%%%%% Full cross-temporal generalization (CTG) matrices (one per region)
%%%%%%  Each is the average CTG (over fifty iterations) 
%%%%%%  Each column is a timebin at which the linear classifer was trained,
%%%%%%  Each row is a timebin at which that linear classifier was tested
%%%%%%  Cool colors represent below chance decoding, (defined as 98th
%%%%%%  percentile of one hundred permuted decoding accuracies), hot colors
%%%%%%  indicate above chance decoding accuracy 
%%%%%%  dashed white line is movement onset
%%%%%%
%%%%%% Subplots that appear to the left of each CTG
%%%%%%  represent a horizontal slice through the decoding matrix:
%%%%%% Accuracy for a given training time, is plotted for all test-bins
%%%%%% (mean +/- 2 SEMs across iterations) 
%%%%%% 
%%%%%% written GRI 250923

close all; clear all ;
set(0,'defaultAxesFontSize',24)
set(0,'defaultAxesFontWeight','bold')


%% intial set-up
subject = 'combined';
nIterations = 50; nUnitsPerRegion = 50;
plotIndividualIterations = false; 
binSize = 0.4;

%% load data 
f2load = dir(['matfiles/Fig5_left.mat']);
ctdat = load([f2load.folder '/' f2load.name]);

%% channels
arrayList = ctdat.arrayList;

%% Create timebase
windowsCenter = ctdat.windowsCenter; 

%% For plotting horizontal slices thru CTD matrix 
time2test = -0.8:0.4:0.8; 
skipFactor = 10;

%% Plot
for aa = 1:length(arrayList)

    data = mean(ctdat.validationScoreStruct.(arrayList{aa}),3);
    sem = std(ctdat.validationScoreStruct.(arrayList{aa}),[],3) ./ sqrt(nIterations);
    cmax = max(sem,[],'all');

    diagData = diag(data);
    diagSEM = diag(sem);

    plotFig = figure('Position',[ 476    86   964   780]);
    zeroPos = find(windowsCenter == 0);
    imagesc(data); hold on; 
    line([zeroPos zeroPos], [1 length(windowsCenter)],'Color', 'w', 'LineStyle',':','LineWidth',2.5); 
    line([1 length(windowsCenter)],[zeroPos zeroPos], 'Color', 'w', 'LineStyle',':','LineWidth',2.5); 

    hold off
    set ( gca, 'ydir', 'normal' );
    xlabel('Generalization Time')
    ylabel('Training Time')
    xticks(1:skipFactor:(length(windowsCenter)))
    xticklabels(windowsCenter(xticks))
    yticks(1:skipFactor:(length(windowsCenter)))
    yticklabels(windowsCenter(xticks))
    
    cb = colorbar; cb.Label.String = 'Accuracy'; colorbarpzn(0.3,1,'full',0.6); 
    set ( gca, 'FontSize', 24);
    title([arrayList{aa} ' - Generalization Across Time'], 'FontSize',28)
    %saveas(plotFig, ['/Users/geena/Dropbox/PhD/Manuscript/figures/' arrayList{aa} '_CTD_' subject '_' popType '.fig']);
    %saveas(plotFig, ['/Users/geena/Dropbox/PhD/Manuscript/figures/' arrayList{aa} '_CTD_' subject '_' popType '.png']);
    %%
    plotFig = figure('units','normalized','outerposition',[  0    0.0886    0.2817    1]);
    nSubplots = length(time2test);
    ha = tight_subplot(nSubplots,1, 0.02, [0.05 0.03], [0.12 0.1]);

    for tt = 1:length(time2test)
        t2testIdx = find(windowsCenter >= time2test(tt),1,'first');
        t2testPerform = data(t2testIdx,:);
        t2testErr = sem(t2testIdx,:);
        
        axes(ha((length(time2test)-tt)+1))
        errorbar(windowsCenter, diagData, diagSEM,'k','linew',2); hold on; 
        errorbar(windowsCenter,t2testPerform, t2testErr,'Color', [0 0.4470 0.7410],'linew',2); hold off
        xlim([-0.8 0.8])
        str = {['Train = ' num2str(time2test(tt)) ' s']}; 
        yyaxis right ;ylabel(str); tmp = gca; tmp.YTick = []; tmp.YLabel.Color = 'k'; tmp.FontSize = 18; tmp.FontWeight = 'bold';
       
        if tt > 1
            tmp.XTickLabel = '';
        end
        box off
    end
    
    han=axes(plotFig,'visible','off'); 
    han.Title.Visible='on'; 
    han.XLabel.Visible='on';
    han.YLabel.Visible='on';
    yyaxis right; 
    yyaxis left; ylabel(han,'Accuracy');han.YLabel.Color = 'k';han.YLabel.FontSize = 20;
    box off; 
    han.YLabel.Position = [-0.1 0.5 0 ];

    xlabel(han,'Test Time Bin','FontSize',20);
    han.XLabel.Color = 'k'; 
    han.XLabel.Position = [0.5 -0.1 0];

    han.Title.Position = [0.5 1.06 0.5];    
    title(han,'Decoder Performance','FontSize',20);
    
    plotFig = gcf;
    plotFig.Color = [1 1 1 ];
    %saveas(plotFig, ['/Users/geena/Dropbox/subFigurePDFs/' arrayList{aa} '_CTDScores_' subject '_' extractBefore(popType, '/') '.fig']);
   % exportgraphics(plotFig, ['/Users/geena/Dropbox/subFigurePDFs/' arrayList{aa} '_CTDScores_' subject '_' extractBefore(popType, '/') '.pdf'], 'ContentType', 'vector');
    %saveas(plotFig, ['/Users/geena/Dropbox/PhD/Manuscript/figures/' arrayList{aa} '_CTDScores_' subject '_' popType '.png']);
    
end


%% =======================================
%%  plot individual iterations for interest 
%% =======================================
if plotIndividualIterations
    arrayList = {'PMv'};
    for i = 1:50

        data = ctdat.validationScoreStruct.(arrayList{1})(:,:,i);

        plotFig = figure;

        imagesc(data)
        set ( gca, 'ydir', 'normal' );
        xlabel('Testing Time')
        ylabel('Training Time')
        xticks(1:2:(length(windowsCenter)))
        xticklabels(windowsCenter(xticks))
        yticks(1:2:(length(windowsCenter)))
        yticklabels(windowsCenter(xticks))
        cb = colorbar; cb.Label.String = 'Accuracy'; caxis([0.55 1]); colormap(jet);
        title('Mean CTD')

        title([subject ' - ' arrayList{1} 'N: ' num2str(i)],'FontSize',20)
        pause;

        %saveas(plotFig, ['~/Desktop/' subject '_' arrayList{1} '_N' num2str(i) '_CTD.png']);
    end
end

%%%% helpers 
function [results, sp2keep] = removeNanChannels(arrayList, results)

for aa = 1:length(arrayList)
    
    chlIdx = results.spNames.(arrayList{aa}){:};

    counter = 1;
    for unit = 1:length(results.spNames.(arrayList{aa}){:})
    
        unitData = squeeze(results.betaScoreStruct.(arrayList{aa})(:,unit,:));
    
        if ~any(unitData,'all') %if it's all zeros, get rid of it
            results.betaScoreStruct.(arrayList{aa})(:,unit,:) = NaN;
           
        else
            sp2keep.(arrayList{aa}){counter,:} = chlIdx(unit); % if not, keep it's tag
            counter = counter + 1;
        end
    
    end


% this nasty line removes the NaN-only cells 
results.betaScoreStruct.(arrayList{aa})(:,~any(squeeze(~any(isnan(results.betaScoreStruct.(arrayList{aa})),1)),2),:) = [];
results.spNames.(arrayList{aa}) = sp2keep.(arrayList{aa});
end
end
