
%%%%%%%%
%%%%%%%% PLOTTING Figure 4A, 4B
%%%%%%%%  Kinematic decoding and non-preserved neural correlations
%%%%%%%%
%%%%%%%% Plots per region decoding accuracies 
%%%%%%%% defined as cross-correlation between actual and predicted face PCs
%%%%%%%% on held out individual trials
%%%%%%%%
%%%%%%%% Versus decoding performance of a shuffled, null dataset 
%%%%%%%% and runs stats, FDR-corr (between regions, also each region v. its null) 
%%%%%%%% Input data contains (region).kinematicDecoding, which is nFolds x nBhvsPred x nIter
%%%%%%%  Also contains bhvTest (true) and bhvPred (predicted) for held out
%%%%%%%  trials; in addition to all final model parameters 
%%%%%%%
%%%%%%%  Shuffled results, similar 
%%%%%%%%% GRI 09/11/2025 


clear all; close all; 
set(0,'defaultAxesFontSize', 18); % bc im blind 

% define the dataset you are running
popSize = '50'; % nCells
nIterations = 25; % pseduopopulations
outcomeMetric = 'CC'; % vs R2 vs RSME etc
saveFlag = false; 
statsPrint = false; 

resultsDir = '/Users/geena/Documents/MATLAB/projects/FacialGesturesPaperCode/FaceGesturePaper/Figure4/matfiles';

%% define the dataset you are running
decodingResults = load([resultsDir '/Fig4AResults.mat' ]); %true results 
decodingShuffleResults =  load([resultsDir '/Fig4ANullResults.mat']);

regionList = fieldnames(decodingResults.resultsOut.combined)';

%% parameters on the model search, not relevant here 
outerFolds = 4; 
statesBhvRelv = [2 3 4 8 10 12 20]; %  N1 gridsearch  
additionalStates = [2 4 8 10]; % ** this is NOT Nx (Nx = statesBhvRelv + additional) **
horizonValues = [5 10 15  ];  % i gridsearch 
nz = 5;                       % bhvDims predicted      
taxis = -0.5:0.05:1.5;         % time axis; didnt save in results file, ugh

% region colormap
regionCmap = [0.4940 0.1840 0.5560	
    0.6350 0.0780 0.1840
    0.8500 0.3250 0.0980	
    0.9290 0.6940 0.1250
    0.5 0.5 0.5];

x = 2*[ones(1,nIterations) ones(1,nIterations)*2 ones(1,nIterations)*3 ones(1,nIterations)*4 ones(1,nIterations)*5] ;

allDecodingMean = [];
allDecodingErr = [];
alldistN1 = [];
alldistNX = [];
allIter2plot = [];

%% pool results 
for rr = 1:length(regionList) 
    theseResults = decodingResults.resultsOut.combined.(regionList{rr});
    shuffledResults = decodingShuffleResults.agg.combined.(regionList{rr});

    outcomeMetric = theseResults.outcomeMetric;

    %%%  kinematic decoding is nFolds x nBhvsPred x nIter
    decodingAcc = squeeze(mean(theseResults.kinematicDecoding,1,'omitnan')); % averaged over nOuterFolds
    [tmp, decodingMean] = std(decodingAcc,[],2);
    
    % nBhvPreds x nRegions
    allDecodingMean(:,rr) = decodingMean;
    allDecodingErr(:,rr) = tmp ./ sqrt(nIterations);

    %%% optimal state #s
    alldistN1 = [alldistN1 theseResults.N1];
    alldistNX = [alldistNX theseResults.NX ];

    %%% Neural, full model self-prediction accuracy
    %%%  neuralSelfFullAcc is nCells x nIterations (averaged over folds) 
    %%% neuralSelfFullMean is nCells (averaged over folds, and iterations) 
    neuralSelfFullAcc = round(squeeze(mean(theseResults.neuralSelfPred2stage,1,'omitnan')),2);
    [neuralSelfFullErr, neuralSelfFullMean] = std(neuralSelfFullAcc,[],2);
    
    %%% allNeuralSelfFullMean is nCells x nRegions (averaged over iterations and folds) 
    allNeuralSelfFullMean(:,rr) = round(neuralSelfFullMean,2);
    allNeuralSelfFullErr(:,rr) = round(neuralSelfFullErr,2);

    %%%% analogously, but for neural, 1st stage model self-prediction
    %%%% accuracy 
    neuralSelfFirstAcc = round(squeeze(mean(theseResults.neuralSelfPred1stage,1,'omitnan')),2);
    [neuralSelfFirstErr, neuralSelfFirstMean] = std(neuralSelfFirstAcc,[],2);
    allNeuralSelfFirstMean(:,rr) = round(neuralSelfFirstMean,2);
    allNeuralSelfFirstErr(:,rr) = round(neuralSelfFirstErr,2);

    % shuffled kinematic decoding is nBhvsPred x nIter
    % allShuffledMean is nBhvsPred x nIter x nRegions
    allShuffled(:,:,rr) = shuffledResults.kinematicDecoding;
    allShuffledStd(:,:,rr) = shuffledResults.kinematicDecodingStd;

end

%% mean and std across iterations for shuffles 
allShuffledMean = squeeze(mean(allShuffled,2));
allShuffledErr = squeeze(mean(allShuffledStd,2));

%% run stats 
%% -------- Significance testing: observed vs null, and between regions --------
% We keep per-iteration observed decoding (nBhvs x nIter x nRegions)
% and per-iteration null (nBhvs x nIter x nRegions).

% 1) Recompute/retain per-iteration observed accuracies to use in tests
nRegions   = numel(regionList);
nBhvs      = size(allShuffled, 1);           % from your null array
nIterations= size(allShuffled, 2);
allObsIter = nan(nBhvs, nIterations, nRegions);

for rr = 1:nRegions
    
    theseResults = decodingResults.resultsOut.combined.(regionList{rr});
    
    % Average across folds -> nBhvs x nIter
    decodingAcc = squeeze(mean(theseResults.kinematicDecoding,1, 'omitnan'));
    % in case "decodingAcc" came out 1 x nBhvs x nIter, fix orientation
    if size(decodingAcc,1) ~= nBhvs
        decodingAcc = squeeze(decodingAcc);   % remove stray dims
        if size(decodingAcc,1) ~= nBhvs
            decodingAcc = decodingAcc.';      % transpose if needed
        end
    end
    allObsIter(:,:,rr) = decodingAcc;
end

% Convenience handles
allNullIter = allShuffled;     % nBhvs x nIter x nRegions

%% ------------- (A) Within-region: observed decoding > null -------------
p_within   = nan(nBhvs, nRegions);
z_within   = nan(nBhvs, nRegions);
r_within   = nan(nBhvs, nRegions);  % 
p_permMean = nan(nBhvs, nRegions);  %

for rr = 1:nRegions
    for bb = 1:nBhvs
        obs = squeeze(allObsIter(bb,:,rr));
        nul = squeeze(allNullIter(bb,:,rr));
        
        % keep only columns present in BOTH
        msk = isfinite(obs) & isfinite(nul);
        xo  = obs(msk);
        xn  = nul(msk);

        if isempty(xo) || isempty(xn)
            p_within(bb,rr) = NaN; z_within(bb,rr)=NaN; r_within(bb,rr)=NaN; p_permMean(bb,rr)=NaN;
            continue
        end

        % Wilcoxon signed-rank on paired iterations (one-sided: obs > null)
        d = xo - xn;
        if all(d==0)
            p = 1; z = 0;
        else
            try
                [p,~,stats] = signrank(xo, xn, 'tail','right');
                z = getfield(stats,'zval'); %
            catch
                % fallback if signrank fails (e.g., too many ties): use sign test
                sPos = sum(d>0); N = sum(d~=0);
                p = 1 - binocdf(sPos-1, N, 0.5);  % one-sided
                z = 0;
            end
        end
        p_within(bb,rr) = p;
        z_within(bb,rr) = z;
        r_within(bb,rr) = z ./ sqrt(sum(d~=0));

        % Auxiliary permutation-style p on the MEAN (coarse with N=iterations)
        mObs = mean(xo, 'omitnan');
        p_permMean(bb,rr) = (sum(xn >= mObs) + 1) / (numel(xn) + 1);
    end
end

%% FDR across all within-region tests (bhv x region)
alpha = 0.01;
p_within_fdr = reshape(bh_fdr_vec(p_within(:), alpha), size(p_within));
sig_within   = p_within_fdr < alpha;

%% ------------- (B) Between regions: do accuracies differ? -------------
pairs = nchoosek(1:nRegions, 2);
nPairs = size(pairs,1);

p_between = nan(nBhvs, nPairs);
z_between = nan(nBhvs, nPairs);
r_between = nan(nBhvs, nPairs);
diffMean  = nan(nBhvs, nPairs);   % mean difference (r1 - r2) for direction

for pi = 1:nPairs
    r1 = pairs(pi,1); r2 = pairs(pi,2);
    for bb = 1:nBhvs
        x1 = squeeze(allObsIter(bb,:,r1));
        x2 = squeeze(allObsIter(bb,:,r2));
        msk = isfinite(x1) & isfinite(x2);
        x1 = x1(msk); x2 = x2(msk);
        d  = x1 - x2;
        diffMean(bb,pi) = mean(d, 'omitnan');

        if isempty(d) || all(d==0)
            p_between(bb,pi)=NaN; z_between(bb,pi)=NaN; r_between(bb,pi)=NaN;
            continue
        end

        try
            [p,~,stats] = signrank(x1, x2, 'tail','both');  % two-sided
            z = getfield(stats,'zval');
        catch
            % fallback sign test (two-sided)
            sPos = sum(d>0); N = sum(d~=0);
            p_one = 1 - binocdf(sPos-1, N, 0.5);
            p = min(1, 2*min(p_one, 1-p_one));  % two-sided from one-sided
            z = 0;
        end
        p_between(bb,pi) = p;
        z_between(bb,pi) = z;
        r_between(bb,pi) = abs(z) ./ sqrt(sum(d~=0));
    end
end

%% FDR across ALL between-region tests (bhv x all region-pairs)
p_between_fdr = reshape(bh_fdr_vec(p_between(:), alpha), size(p_between));
sig_between   = p_between_fdr < alpha;

%% -------- Package results for plotting --------
decodingStats = struct();
decodingStats.within.regionList   = regionList(:);
decodingStats.within.p            = p_within;
decodingStats.within.p_fdr        = p_within_fdr;
decodingStats.within.sig_fdr      = sig_within;
decodingStats.within.z            = z_within;
decodingStats.within.r            = r_within;        % effect size
decodingStats.within.p_perm_mean  = p_permMean;      

pairNames = cell(nPairs,1);
for i = 1:nPairs
    pairNames{i} = sprintf('%s_vs_%s', regionList{pairs(i,1)}, regionList{pairs(i,2)});
end
decodingStats.between.pairs       = pairNames;
decodingStats.between.p           = p_between;
decodingStats.between.p_fdr       = p_between_fdr;
decodingStats.between.sig_fdr     = sig_between;
decodingStats.between.z           = z_between;
decodingStats.between.r           = r_between;
decodingStats.between.diffMean    = diffMean;

%% (Optional) quick textual summary
if statsPrint
    fprintf('\nWithin-region (observed > null), FDR α=%.2f:\n', alpha);
    for rr = 1:nRegions
        for bb = 1:nBhvs
            fprintf('  %s | pred %d: p=%.3g (FDR=%.3g), r=%.2f, p_permMean=%.3g\n', ...
                regionList{rr}, bb, ...
                p_within(bb,rr), p_within_fdr(bb,rr), r_within(bb,rr), ...
                p_permMean(bb,rr));
        end
    end

    fprintf('\nBetween-region differences (two-sided), FDR α=%.2f:\n', alpha);
    for pi = 1:nPairs
        for bb = 1:nBhvs
            if sig_between(bb,pi)
                fprintf('  %s | pred %d: p=%.3g (FDR=%.3g), Δ=%.3f\n', ...
                    pairNames{pi}, bb, p_between(bb,pi), p_between_fdr(bb,pi), diffMean(bb,pi));
            end
        end
    end
end

%% plot it all with stats (supplemental figure) 
plot_decoding_bars_with_sig(allDecodingMean, allDecodingErr, ...
    allShuffledMean, allShuffledErr, ...
    regionList, regionCmap, outcomeMetric, decodingStats, ...
    'Alpha', 0.01, 'MaxPairsPerGroup', 3, 'YLim', [0 0.8]);

%% plot all decoding Mean (without stats, Fig 4A) 
figure;
b = bar(allDecodingMean);
for rr = 1:length(regionList)
    b(rr).FaceColor = regionCmap(rr,:);
    b(rr).FaceAlpha = 0.5;
end
hold on;

%% decoding mean error bars 
ngroups = size(allDecodingMean, 1); nbars = size(allDecodingMean, 2); groupwidth = min(0.8, nbars/(nbars + 1.5));
for i = 1:nbars
    tmp = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    er = errorbar(tmp, allDecodingMean(:,i), allDecodingErr(:,i),'HandleVisibility','off');
    er.Color = [0 0 0]; er.LineStyle = 'none'; er.LineWidth = 1.5;
end

%% plot the shuffled decoding accuracies 

b2 = bar(squeeze(allShuffledMean)); hold off
for rr = 1:length(regionList)
    b2(rr).FaceColor = 'k';
    b2(rr).BarWidth = 0.4;
    b2(rr).FaceAlpha = 0.5;
end
hold on;

%% shuffled decoding accuracy error bars 
ngroups = size(allShuffledMean, 1); nbars = size(allShuffledMean, 2); groupwidth = min(0.8, nbars/(nbars + 1.5));

for i = 1:nbars
    tmp = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    er = errorbar(tmp, allShuffledMean(:,i), 2*allShuffledErr(:,i),'HandleVisibility','off');
    er.Color = [0 0 0]; er.LineStyle = 'none'; er.LineWidth = 1.5;
end

ylim([0 0.8])
legend(regionList); xlabel('Face Motion PC'); ylabel('Cross Correlation')
title('Kinematic Decoding Performance');

%% Fig 4B ; plot single iteration result (predicted vs actual face PC) 
[~, iter2plot] = max(decodingAcc,[],2); ii = mode(iter2plot);

for rr = length(regionList) 
    
    theseResults = decodingResults.resultsOut.combined.(regionList{rr});

    outerBhvTestC = cell2mat(theseResults.bhvTest{ii});
    outerBhvPredC = cell2mat(theseResults.bhvPred{ii});
    outerBhvTest = theseResults.bhvTest{ii};
    
    heldOutBhvTest = reshape(outerBhvTestC,length(taxis),length(outerBhvTest),nz);
    heldOutBhvPred = reshape(outerBhvPredC,length(taxis),length(outerBhvTest),nz);
    
    figure('units','normalized','outerposition',[0 0 1 1])
    ha = tight_subplot(1,3,0.125,[0.075 0.1],[0.1 0.05]);

    for bb = 1:nz-2 % each predicted bhvPC
        
        test1Y = mean(heldOutBhvTest(:,:,bb),2);
        test1Err = 2*(std(heldOutBhvTest(:,:,bb),[],2) ./ sqrt(length(outerBhvTest)));
    
        pred1Y = mean(heldOutBhvPred(:,:,bb),2);
        pred1Err = 2*(std(heldOutBhvPred(:,:,bb),[],2) ./ sqrt(length(outerBhvTest)));
        
        axes(ha(bb))
        shadedErrorBar(taxis,test1Y,test1Err,'lineprops',{'Color','k','linew',8});
        hold on;
        shadedErrorBar(taxis, pred1Y,pred1Err,'lineprops',{'Color',[0.3010 0.5 0.9],'linew',8});
        hold off

        title(['Face PC' num2str(bb)],'FontSize',26)
        set(ha(bb), 'YTickLabelMode','auto'); legend ('actual', 'predicted')
        xlim([-0.5 1.5])
        if bb == 1 
            xlabel('Time, s'); ylabel('a.u.'); 
            legend('true','pred')
        end

    end

    set(ha(1:nz-2),'XTickLabel',[-0.5 0 0.5 1 1.5]); 
    sgtitle([regionList{rr} 'Iter = ' num2str(ii)],'FontSize',36,'FontWeight','bold')
end

%% NOT included in paper:
%% for each region, plot the decoding accuracy of the NEURAL data ('allNeuralSelf')
%% for both the full model (two stage, behaviorally relevant + irrelevant dynamics) 
%% and for the limited model (1 stage, behaviorally relevant only) 
%% the ratio of these neural self-prediction accuracies ranges 0-1 , where 1 indicates 
%% ALL behaviorally relevant info (here, face kinematics), and 0 indicates all 'other' info
ratio = mean(allNeuralSelfFirstMean,1) ./ mean(allNeuralSelfFullMean,1);
ratioErr = mean(allNeuralSelfFirstErr,1);

% neural self prediction FULL model
figure;
subplot 131
x = categorical(regionList); x = reordercats(x, {'S1','M1','PMv','M3','All'});
b = bar(x, mean(allNeuralSelfFullMean,1),'facecolor','flat'); b.CData = regionCmap;
hold on; 

er = errorbar(x,mean(allNeuralSelfFullMean,1), mean(allNeuralSelfFullErr,1),'HandleVisibility','off');
er.Color = [0 0 0]; er.LineStyle = 'none'; er.LineWidth = 3;
hold off; 
ylabel(outcomeMetric);ylim([0 1]);
title('Neural Self Predict, Full Model ');

% neural self prediction FIRST stage only 
subplot 132;
x = categorical(regionList); x = reordercats(x, {'S1','M1','PMv','M3','All'});
b = bar(x, mean(allNeuralSelfFirstMean,1),'facecolor','flat'); b.CData = regionCmap;
hold on; 

er = errorbar(x,mean(allNeuralSelfFirstMean,1), mean(allNeuralSelfFirstErr,1),'HandleVisibility','off');
er.Color = [0 0 0]; er.LineStyle = 'none'; er.LineWidth = 3;
hold off; 
ylabel(outcomeMetric); ylim([0 1]);
title('Neural Self Predict, First Stage Only Model ');

% ratio 
subplot 133;
x = categorical(regionList); x = reordercats(x, {'S1','M1','PMv','M3','All'});
b = bar(x, ratio,'facecolor','flat'); b.CData = regionCmap;
hold on; 

er = errorbar(x,ratio, ratioErr,'HandleVisibility','off');
er.Color = [0 0 0]; er.LineStyle = 'none'; er.LineWidth = 3;
hold off; 
ylim([0 1]);
title('Ratio: Kinematic/All Other Information');

disp(['Average Neural, Full Model        ' outcomeMetric ' : ' num2str(mean(allNeuralSelfFullMean,1))]);
disp(['Average Neural First-stage only   ' outcomeMetric ' : ' num2str(mean(allNeuralSelfFirstMean,1))]);
disp(['Kinematic : AllOther Info Ratio = ' num2str(ratio)]);
disp(['                                 +/-' num2str(ratioErr)])




%% plotting helper 
function plot_decoding_bars_with_sig(allDecodingMean, allDecodingErr, ...
    allShuffledMean, allShuffledErr, regionList, regionCmap, ...
    outcomeMetric, decodingStats, varargin)
% Plot decoding (colored bars) vs. shuffle null (grey/black bars) with significance.
%
% Inputs match your arrays:
%   allDecodingMean : [nBhvs x nRegions]
%   allDecodingErr  : [nBhvs x nRegions]  (SEM across iterations)
%   allShuffledMean : [nBhvs x nRegions]
%   allShuffledErr  : [nBhvs x nRegions]  (SEM across iterations; you used 2*SEM in plot)
%   regionList      : 1 x nRegions cellstr
%   regionCmap      : nRegions x 3 RGB in 0–1
%   outcomeMetric   : e.g., 'CC'
%   decodingStats   : struct from the earlier stats block
%
% Name-Value (optional):
%   'Alpha'              : 0.05 (for labeling only)
%   'MaxPairsPerGroup'   : 3    (limit # of between-region connectors per behavior)
%   'YLim'               : []   (set y-limits; auto if empty)

ip = inputParser; ip.FunctionName = mfilename;
addParameter(ip,'Alpha',0.05,@(x) isscalar(x) && x>0 && x<=1);
addParameter(ip,'MaxPairsPerGroup',3,@(x) isscalar(x) && x>=0);
addParameter(ip,'YLim',[],@(x) isempty(x) || (isnumeric(x) && numel(x)==2));
parse(ip,varargin{:});
alpha = ip.Results.Alpha;
maxPairs = ip.Results.MaxPairsPerGroup;
yL = ip.Results.YLim;

[nBhvs, nRegions] = size(allDecodingMean);

%% --- 1) Observed bars (your template) ---
figure; 
b = bar(allDecodingMean); hold on;
for rr = 1:nRegions
    b(rr).FaceColor = regionCmap(rr,:);
    b(rr).FaceAlpha = 0.5;
end

% error bars on observed means
ngroups = size(allDecodingMean, 1);
nbars   = size(allDecodingMean, 2);
groupwidth = min(0.8, nbars/(nbars + 1.5));
xObs = nan(nBhvs, nRegions);  % store x positions of each colored bar

for i = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    xObs(:,i) = x(:);
    er = errorbar(x, allDecodingMean(:,i), allDecodingErr(:,i), 'HandleVisibility','off');
    er.Color = [0 0 0]; er.LineStyle = 'none'; er.LineWidth = 1.5;
end

%% --- 2) Null bars (your template) ---
b2 = bar(squeeze(allShuffledMean)); 
for rr = 1:nRegions
    b2(rr).FaceColor = [0 0 0];
    b2(rr).BarWidth = 0.4;
    b2(rr).FaceAlpha = 0.25; % slightly lighter so colored bars pop
    b2(rr).HandleVisibility = 'off'; % keep legend clean
end

% error bars on null means (you used 2*SEM here)
ngroupsN = size(allShuffledMean, 1); nbarsN = size(allShuffledMean, 2);
groupwidthN = min(0.8, nbarsN/(nbarsN + 1.5));
for i = 1:nbarsN
    x = (1:ngroupsN) - groupwidthN/2 + (2*i-1) * groupwidthN / (2*nbarsN);
    er = errorbar(x, allShuffledMean(:,i), 2*allShuffledErr(:,i),'HandleVisibility','off');
    er.Color = [0 0 0]; er.LineStyle = 'none'; er.LineWidth = 1.5;
end

legend(regionList, 'Location','best'); 
xlabel('Face Motion PC'); ylabel('Cross Correlation, (CC)');
title('Kinematic Decoding Performance');
grid on;

% baseline y-limits
yMaxBars = max([allDecodingMean(:)+allDecodingErr(:); allShuffledMean(:)+2*allShuffledErr(:)], [], 'omitnan');
yl = [0, max(0.8, yMaxBars*1.15)]; % default padding
if ~isempty(yL), yl = yL; end
ylim(yl);

%% --- 3) Significance overlays from decodingStats ---
% A) within-region (obs > null): stars above colored bars
if isfield(decodingStats,'within') && isfield(decodingStats.within,'p_fdr')
    pW = decodingStats.within.p_fdr;  % [nBhvs x nRegions]
    for rr = 1:nRegions
        for bb = 1:nBhvs
            p = pW(bb, rr);
            if ~isfinite(p), continue; end
            stars = p2stars(p);
            if isempty(stars), continue; end

            x = xObs(bb, rr);
            y = allDecodingMean(bb, rr) + allDecodingErr(bb, rr) + 0.02*range(yl);
            text(x, y, stars, 'HorizontalAlignment','center', 'VerticalAlignment','bottom', ...
                'FontWeight','bold', 'FontSize',18, 'Color', regionCmap(rr,:), 'Clipping','off');
        end
    end
end

% B) between-region: pairwise connectors per behavior (limit maxPairs per group)
if isfield(decodingStats,'between') && isfield(decodingStats.between,'p_fdr')
    pB = decodingStats.between.p_fdr;     % [nBhvs x nPairs]
    pairNames = decodingStats.between.pairs(:);
    % Reconstruct index pairs (assumes same ordering used to build stats)
    pairs = nchoosek(1:nRegions, 2);      % size = nPairs x 2

    yPad = 0.03*range(yl);
    for bb = 1:nBhvs
        % pick significant pairs for this behavior
        sigIdx = find(isfinite(pB(bb,:)) & pB(bb,:) < alpha);
        % order by p-value (smallest first) and cap
        [~, ord] = sort(pB(bb, sigIdx), 'ascend');
        sigIdx = sigIdx(ord);
        if ~isempty(maxPairs), sigIdx = sigIdx(1:min(maxPairs, numel(sigIdx))); end

        % base y at top of the tallest observed bar in this behavior group
        yBase = max(allDecodingMean(bb,:)+allDecodingErr(bb,:), [], 'omitnan');
        yH = yBase + 0.05*range(yl);

        for kk = 1:numel(sigIdx)
            pi = sigIdx(kk);
            r1 = pairs(pi,1); r2 = pairs(pi,2);
            x1 = xObs(bb, r1); x2 = xObs(bb, r2);
            if x1 == 0 || x2 == 0, continue; end

            y = yH + (kk-1)*yPad;        % stack connectors to avoid overlap
            % draw connector
            plot([x1 x1 x2 x2], [y-0.004 y y y-0.004], 'k-', 'LineWidth', 1.0, 'HandleVisibility','off');
            % star(s) above midpoint
            pm = pB(bb, pi);
            stars = p2stars(pm);
            text(mean([x1 x2]), y + 0.004, stars, 'HorizontalAlignment','center', ...
                 'VerticalAlignment','bottom', 'FontSize',14, 'FontWeight','bold', 'Color', [0 0 0], 'Clipping','off');
        end
    end
end

% expand ylim if annotations exceed current top
ax = gca;
drawnow;
ylNow = ylim(ax);
allText = findall(ax, 'Type','text');
if ~isempty(allText)
    yText = arrayfun(@(h) h.Position(2), allText);
    yTop = max([ylNow(2); yText(:)]);
    if yTop > ylNow(2)
        ylim([ylNow(1) yTop + 0.02*range(ylNow)]);
    end
end

end

%% --------- helpers ----------
function s = p2stars(p)
% Convert p-value to asterisks (FDR-corrected p recommended)
if ~isfinite(p) || p>=0.05
    s = '';
elseif p<0.001
    s = '***';
elseif p<0.01
    s = '**';
else
    s = '*';
end
end

function p_fdr = bh_fdr_vec(p, alpha)
    if nargin < 2, alpha = 0.05; end
    p = p(:);
    m = numel(p);
    [ps, idx] = sort(p);
    adj = nan(m,1);
    minv = 1;
    for i = m:-1:1
        if ~isfinite(ps(i)), adj(i) = NaN; continue; end
        val = ps(i) * m / i;
        minv = min(minv, val);
        adj(i) = min(minv, 1);
    end
    p_fdr = nan(size(p));
    p_fdr(idx) = adj;
end
