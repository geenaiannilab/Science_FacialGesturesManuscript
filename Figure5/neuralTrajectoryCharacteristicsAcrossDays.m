%%%%%%% PLOTTING FIGURE 5 (right panels) 
%%%%%%% Face gesture neural trajectory velocities differ across regions 
%%%%%%  Velocity (rate of change) of facial gesture trajectories in the neural space defined by the top 12 principle components of neural activity in each region.
%%%%%%% The velocity for each gesture is plotted, as is the all-gesture average (in dashed black line). 
%%%%%%%
%%%%%%% Population trajectory velocity was obtained by 
%%%%%%% 1) obtaining the neural trajectories for each facial gesture in
%%%%%%% each region (as in Figure 5B), averaged over iterations (N=50)
%%%%%%% 2) calculating the distance between two successive time bins in N-dimensional neural space 
%%%%%%  3) to test additional dimensions, edit 'dims2test' below
%%%%%% 
%%%%%   Also runs permutation-based clustering significance testing to
%%%%%%%     compare region-region trajectory velocities; horizontal black
%%%%%%%     bars under region-region velocity comparisons indicate
%%%%%%%     significant differences in regions' speeds 
%%%%%%% Input data -- both velocity and distance metrics are nTimepoints x
%%%%%%%         nPCs x nIterations 
%%%%%%%
%%%%%%% This will also plot the average inter-gesture trajectory distances
%%%%%%
 

clear all;
set(0,'defaultAxesFontSize',24)
set(0,'defaultAxesFontWeight','bold')

%% load data
data = load(['Matfiles/neuralTrajectory_combined.mat']);
subject2analyze = {'combined'};

%% extract parameters
bhvs = {'Thr','LS','Chew'};
dims2test = [12]; % [3 6 12 20];

%% pairwise stats/plotting parameters 
dimUse     = 12;          % which dimensionality to analyze/plot
gest       = 'All';       % 'Thr' | 'LS' | 'Chew' | 'All'
usePaired  = false;        % true = paired sign-flip; false = Welch two-sample
alpha_cf   = 0.01;        % per-bin cluster-forming alpha (Welch only)
alpha_cl   = 0.01;        % cluster-level alpha (paired function uses this)
nPerm      = 5000;
minCluster = 4;
tail       = 'both';      % for Welch version if usePaired=false
colorArray = [0.4940 0.1840 0.5560;0.6350 0.0780 0.1840;0.8500 0.3250 0.0980;0.9290 0.6940 0.1250];

%% define the data 
nIterations = data.nIterations;
arrayList = fieldnames(data.velChewTotal);
win = data.win;
tmin = data.tmin2take;
tmax = data.tmax2take;
taxis = -tmin:win:tmax;  


%% calculate distance between gesture-trajectories, also global average
for ss = 1:length(subject2analyze)

    thisSubj = subject2analyze{ss};

    for aa = 1:length(arrayList)

        distanceThrVLS.(arrayList{aa}).mean = mean(data.distanceThrVLS.(arrayList{aa}).data,3);
        distanceThrVLS.(arrayList{aa}).sem = 2*(std(data.distanceThrVLS.(arrayList{aa}).data,[],3)) ./ sqrt(nIterations);
 
        distanceLSVCh.(arrayList{aa}).mean = mean(data.distanceLSVCh.(arrayList{aa}).data,3);
        distanceLSVCh.(arrayList{aa}).sem = 2*(std(data.distanceLSVCh.(arrayList{aa}).data,[],3)) ./ sqrt(nIterations);

        distanceThrVCh.(arrayList{aa}).mean = mean(data.distanceThrVCh.(arrayList{aa}).data,3);
        distanceThrVCh.(arrayList{aa}).sem = 2*(std(data.distanceThrVCh.(arrayList{aa}).data,[],3)) ./ sqrt(nIterations);

        thrVelocity.(arrayList{aa}).mean = mean(data.velThrTotal.(arrayList{aa}).data,3);
        thrVelocity.(arrayList{aa}).sem = 2*(std(data.velThrTotal.(arrayList{aa}).data,[],3)) ./ sqrt(nIterations);

        lsVelocity.(arrayList{aa}).mean = mean(data.velLSTotal.(arrayList{aa}).data,3);
        lsVelocity.(arrayList{aa}).sem = 2*(std(data.velLSTotal.(arrayList{aa}).data,[],3)) ./ sqrt(nIterations);

        chewVelocity.(arrayList{aa}).mean = mean(data.velChewTotal.(arrayList{aa}).data,3);
        chewVelocity.(arrayList{aa}).sem = 2*(std(data.velChewTotal.(arrayList{aa}).data,[],3)) ./ sqrt(nIterations);

    end % end all arrays 

    %% plot distance between gesture-trajectories, also global average, using all dimensions
    %% plot velocity of gesture-trajectories, also global average, using all dimensions

    fig1 = figure('Position', [476    83   975   783]);
    for aa = 1:length(arrayList)
        subplot(2,2,aa)
        errorbar(taxis(1:end-1), thrVelocity.(arrayList{aa}).mean(:,end), thrVelocity.(arrayList{aa}).sem(:,end) ,'Color','r','linew',2); hold on
        errorbar(taxis(1:end-1), lsVelocity.(arrayList{aa}).mean(:,end), lsVelocity.(arrayList{aa}).sem(:,end),'Color','b','linew',2);
        errorbar(taxis(1:end-1), chewVelocity.(arrayList{aa}).mean(:,end), chewVelocity.(arrayList{aa}).sem(:,end) ,'Color','g','linew',2); 
        avgTraj = mean([thrVelocity.(arrayList{aa}).mean(:,end), lsVelocity.(arrayList{aa}).mean(:,end), chewVelocity.(arrayList{aa}).mean(:,end)],2);

        plot(taxis(1:end-1), avgTraj,'--k','linew',5);
        ylabel('Speed, au'); xlabel('Time, s'); ylim([0 3]); xlim([-0.8 0.8]);
        title([arrayList{aa}])
    end
    sgtitle('Neural Trajectory Velocities, All Dimensions','FontSize',28,'FontWeight','bold');
 

    fig3 = figure('Position', [476    83   975   783]);
    for aa = 1:length(arrayList)
        avgDistance = mean([distanceThrVLS.(arrayList{aa}).mean(:,end) distanceThrVCh.(arrayList{aa}).mean(:,end) distanceLSVCh.(arrayList{aa}).mean(:,end)],2);
        subplot(2,2,aa)
        errorbar(taxis, distanceThrVLS.(arrayList{aa}).mean(:,end), distanceThrVLS.(arrayList{aa}).sem(:,end) ,'Color',[0.3010 0.7450 0.9330],'linew',2); hold on
        errorbar(taxis, distanceThrVCh.(arrayList{aa}).mean(:,end), distanceThrVCh.(arrayList{aa}).sem(:,end) ,'Color',[0.4660 0.6740 0.1880],'linew',2);
        errorbar(taxis, distanceLSVCh.(arrayList{aa}).mean(:,end), distanceLSVCh.(arrayList{aa}).sem(:,end) ,'Color',[0.6350 0.0780 0.1840]	,'linew',2);
        plot(taxis, avgDistance,'k','linew',2); hold off
        ylabel('Distance, au'); xlabel('Time, s'); ylim([5 28 ])
        title([arrayList{aa} '; All Dims'])
        legend('Thr-LS','Thr-Ch','LS-Ch');
    end
    sgtitle([subject2analyze{:} ' All Neural Distances'],'FontSize',20);

    fig4 = figure('Position', [476    83   975   783]);
    for aa = 1:length(arrayList)
        avgDistance = mean([distanceThrVLS.(arrayList{aa}).mean(:,end) distanceThrVCh.(arrayList{aa}).mean(:,end) distanceLSVCh.(arrayList{aa}).mean(:,end)],2);
        avgSEM = sum([distanceThrVLS.(arrayList{aa}).sem(:,end) distanceThrVCh.(arrayList{aa}).sem(:,end) distanceLSVCh.(arrayList{aa}).sem(:,end)],2);
        errorbar(taxis,avgDistance, avgSEM,'Color', colorArray(aa,:),'linew',3.5); hold on
        ylim([0 25]); ylabel('Distance, au'); xlabel('Time, s')

        legend('S1','M1','PMv','M3');
    end
    sgtitle(['Average Intra-Trajectory Distances; All Dimensions'],'FontSize',28);


    %% repeat for a number of different dimensions
    for dd = 1:length(dims2test)

    fig5 = figure('Position', [476    83   975   783]);
        for aa = 1:length(arrayList)
            avgDistance = mean([distanceThrVLS.(arrayList{aa}).mean(:,dd) distanceThrVCh.(arrayList{aa}).mean(:,dd) distanceLSVCh.(arrayList{aa}).mean(:,dd)],2);
            avgSEM = sum([distanceThrVLS.(arrayList{aa}).sem(:,dd) distanceThrVCh.(arrayList{aa}).sem(:,dd) distanceLSVCh.(arrayList{aa}).sem(:,dd)],2);
            errorbar(taxis,avgDistance, avgSEM,'Color', colorArray(aa,:),'linew',3.5); hold on
            ylim([0 25]); ylabel('Distance, au'); xlabel('Time, s')
            
            legend('S1','M1','PMv','M3');
        end
       sgtitle(['Average Intra-Trajectory Distances; Dims = ' num2str(dims2test(dd)) ],'FontSize',28);

       fig6 = figure('Position', [476    83   975   783]);
       for aa = 1:length(arrayList)
           avgDistance = mean([distanceThrVLS.(arrayList{aa}).mean(:,dd) distanceThrVCh.(arrayList{aa}).mean(:,dd) distanceLSVCh.(arrayList{aa}).mean(:,dd)],2);
           subplot(2,2,aa)
            errorbar(taxis, distanceThrVLS.(arrayList{aa}).mean(:,dd), distanceThrVLS.(arrayList{aa}).sem(:,dd) ,'Color',[0.3010 0.7450 0.9330],'linew',2); hold on
            errorbar(taxis, distanceThrVCh.(arrayList{aa}).mean(:,dd), distanceThrVCh.(arrayList{aa}).sem(:,dd) ,'Color',[0.4660 0.6740 0.1880],'linew',2);
            errorbar(taxis, distanceLSVCh.(arrayList{aa}).mean(:,dd), distanceLSVCh.(arrayList{aa}).sem(:,dd) ,'Color',[0.6350 0.0780 0.1840]	,'linew',2); 
            plot(taxis, avgDistance,'k','linew',2); hold off
            ylim([0 25]); ylabel('Distance, au'); xlabel('Time, s')
            title([arrayList{aa} '; Dims = ' num2str(dims2test(dd))])
            legend('Thr-LS','Thr-Ch','LS-Ch');
        end
        sgtitle([subject2analyze{:} ' All Neural Distances'],'FontSize',20);

        fig7 = figure('Position', [476    83   975   783]);
        for aa = length(arrayList):-1:1
            avgSEM = mean([thrVelocity.(arrayList{aa}).sem(:,dims2test(dd)), lsVelocity.(arrayList{aa}).sem(:,dims2test(dd)), chewVelocity.(arrayList{aa}).sem(:,dims2test(dd))],2);
            avgTraj2 = mean([thrVelocity.(arrayList{aa}).mean(:,dims2test(dd)), lsVelocity.(arrayList{aa}).mean(:,dims2test(dd)), chewVelocity.(arrayList{aa}).mean(:,dims2test(dd))],2);
            errorbar(taxis(1:end-1), avgTraj2,avgSEM, 'Color',colorArray(aa,:),'linew',5); hold on
            ylabel('Speed, au'); xlabel('Time, s'); ylim([0 2]); xlim([-0.8 0.8]); 
            l = legend('M3','PMv','M1','S1');
        end
        hold off
        sgtitle(['Neural Trajectory Velocities, Dim = ' num2str(dims2test(dd))],'FontSize',28,'FontWeight','bold');


        fig8 = figure('Position', [476    83   975   783]);
        for aa = 1:length(arrayList)
            subplot(2,2,aa)

            errorbar(taxis(1:end-1), thrVelocity.(arrayList{aa}).mean(:,dims2test(dd)), thrVelocity.(arrayList{aa}).sem(:,dims2test(dd)) ,'Color','r','linew',2); hold on
            errorbar(taxis(1:end-1), lsVelocity.(arrayList{aa}).mean(:,dims2test(dd)), lsVelocity.(arrayList{aa}).sem(:,dims2test(dd)),'Color','b','linew',2);
            errorbar(taxis(1:end-1), chewVelocity.(arrayList{aa}).mean(:,dims2test(dd)), chewVelocity.(arrayList{aa}).sem(:,dims2test(dd)) ,'Color','g','linew',2); 
            avgTraj3 = mean([thrVelocity.(arrayList{aa}).mean(:,dims2test(dd)), lsVelocity.(arrayList{aa}).mean(:,dims2test(dd)), chewVelocity.(arrayList{aa}).mean(:,dims2test(dd))],2);
            plot(taxis(1:end-1), avgTraj3,'--k','linew',5);       
            ylabel('Speed, au'); xlabel('Time, s'); ylim([0 2.2]); xlim([-0.8 0.8]);
            title([arrayList{aa}])

        end
       sgtitle(['Neural Trajectory Speeds, ' num2str(dims2test(dd)) '-D'],'FontSize',28,'FontWeight','bold');

    end % end dimensionsList

    %% --- Pairwise comparison of velocities (one subplot per region pair) ---
    figPairs = plot_pairwise_velocities_with_stats( ...
        data, arrayList, taxis, colorArray, dimUse, gest, ...
        usePaired, alpha_cf, alpha_cl, nPerm, minCluster, tail);


end % end subjectList


%%%%%%%%% ======== CLUSTER PERMUTATION HELPERS ======== %%%%%%%%%%%

function out = cluster_perm_2group(A, B, varargin)
% Two-group cluster permutation test across time (Welch's t)
%
% A, B : [nSamples x T] matrices (rows = samples, cols = time)
% Optional parameters:
%   'tail'      : 'both' (default) | 'right' | 'left'
%   'alpha_cf'  : cluster-forming alpha (default 0.05)
%   'nPerm'     : number of permutations (default 5000)
%
% Returns struct with fields:
%   t, p_unc, sig_mask, clusters, maxMass0, info

    p = inputParser;
    addParameter(p,'tail','both',@(s)ischar(s)||isstring(s));
    addParameter(p,'alpha_cf',0.05,@(x)isnumeric(x)&&isscalar(x)&&x>0&&x<1);
    addParameter(p,'nPerm',5000,@(x)isnumeric(x)&&isscalar(x)&&x>=100);
    parse(p,varargin{:});
    tail     = lower(string(p.Results.tail));
    alpha_cf = p.Results.alpha_cf;
    nPerm    = p.Results.nPerm;

    % bookkeeping
    nA_t = sum(isfinite(A),1);
    nB_t = sum(isfinite(B),1);
    keep = (nA_t >= 2) & (nB_t >= 2);
    A = A(:,keep);  B = B(:,keep);
    [nA,T] = size(A);  nB = size(B,1);

    % Observed t, df, and two-sided p (per time)
    [tObs, p2Obs, dfObs] = welch_t(A,B);

    % cluster-forming mask from per-time p-values
    switch tail
        case "both"
            supra = (p2Obs < alpha_cf);
            massfun = @(t) sum(abs(t));
            tForMass = abs(tObs);
        case "right"
            pRight = 1 - tcdf(tObs, dfObs);
            supra  = (pRight < alpha_cf);
            massfun = @(t) sum(t);
            tForMass = tObs;
        case "left"
            pLeft = tcdf(tObs, dfObs);
            supra = (pLeft < alpha_cf);
            massfun = @(t) sum(t);
            tForMass = -tObs;
        otherwise
            error('tail must be both/right/left');
    end

    obsClusters = label_clusters(supra);
    for c = 1:numel(obsClusters)
        seg = obsClusters(c).start:obsClusters(c).end;
        obsClusters(c).mass = massfun(tForMass(seg));
    end

    % permutation: shuffle whole trajectories
    X = [A; B]; nTot = size(X,1);
    maxMass0 = zeros(nPerm,1);
    for pidx = 1:nPerm
        idx = randperm(nTot);
        Ai = X(idx(1:nA),:);
        Bi = X(idx(nA+1:end),:);
        [t0, p20, df0] = welch_t(Ai,Bi);
        switch tail
            case "both"
                supra0 = (p20 < alpha_cf);
                mass0  = max_cluster_mass(abs(t0), supra0, @(u) sum(u));
            case "right"
                pR0    = 1 - tcdf(t0, df0);
                supra0 = (pR0 < alpha_cf);
                mass0  = max_cluster_mass(t0, supra0, @(u) sum(u));
            case "left"
                pL0    = tcdf(t0, df0);
                supra0 = (pL0 < alpha_cf);
                mass0  = max_cluster_mass(-t0, supra0, @(u) sum(u));
        end
        maxMass0(pidx) = mass0;
    end

    sig_mask = false(1,T);
    for c = 1:numel(obsClusters)
        cmass = obsClusters(c).mass;
        pcl = (sum(maxMass0 >= cmass) + 1) / (nPerm + 1);
        obsClusters(c).p = pcl;
        if pcl < 0.05
            sig_mask(obsClusters(c).start:obsClusters(c).end) = true;
        end
    end

    % pack output
    fullT = false(1, numel(keep)); fullT(keep) = true;
    out.t        = nan(1, numel(keep)); out.t(keep) = tObs;
    out.p_unc    = nan(1, numel(keep)); out.p_unc(keep) = p2Obs;
    out.sig_mask = false(1, numel(keep)); out.sig_mask(keep) = sig_mask;
    out.clusters = obsClusters;
    out.maxMass0 = maxMass0;
    out.info     = struct('tail',tail,'alpha_cf',alpha_cf,'nPerm',nPerm,'nA',nA,'nB',nB,'T',T,'keptMask',fullT);

end

function [clusters, sigMask, pvals, tvals] = paired_cluster_perm_1d(VA, VB, alpha, nPerm, minCluster)
% Paired cluster permutation test across 1-D time series
% ------------------------------------------------------
% VA, VB : [T x N] matrices (T = time bins, N = paired samples)
% alpha  : cluster-level alpha (e.g. 0.05)
% nPerm  : number of permutations (e.g. 5000)
% minCluster : minimum contiguous bins to define a cluster
%
% Returns:
%   clusters : struct array with fields idx, sumT, p
%   sigMask  : 1 x T logical mask of significant time bins
%   pvals    : cluster p-values
%   tvals    : 1 x T t-statistics per bin

    [T, N] = size(VA);
    if any(size(VB) ~= [T N])
        error('VA and VB must be T x N with the same number of samples (paired).');
    end

    % Paired differences
    D = VA - VB; 
    mD = mean(D, 2);
    sD = std(D, 0, 2) + eps;
    tvals = (mD ./ (sD ./ sqrt(N)))'; % 1 x T
    tcrit = tinv(1 - alpha/2, N - 1);

    supra = abs(tvals) > tcrit;
    supra = prune_short_segments(supra, minCluster);

    % --- Observed clusters ---
    obsSegs = mask2segments(supra);
    obsScores = zeros(size(obsSegs,1),1);
    for k = 1:size(obsSegs,1)
        idx = obsSegs(k,1):obsSegs(k,2);
        obsScores(k) = sum(abs(tvals(idx)));
    end

    % --- Permutation null ---
    maxNull = zeros(nPerm,1);
    for p = 1:nPerm
        flips = (rand(1,N) > 0.5)*2 - 1;  % random ±1 sign flips
        Dp = D .* flips;
        mDp = mean(Dp,2);
        sDp = std(Dp,0,2) + eps;
        tp = (mDp ./ (sDp ./ sqrt(N)))';
        suprap = abs(tp) > tcrit;
        suprap = prune_short_segments(suprap, minCluster);
        segp = mask2segments(suprap);
        if isempty(segp)
            maxNull(p) = 0;
        else
            scores = zeros(size(segp,1),1);
            for kk = 1:size(segp,1)
                idx = segp(kk,1):segp(kk,2);
                scores(kk) = sum(abs(tp(idx)));
            end
            maxNull(p) = max(scores);
        end
    end

    % --- Cluster p-values ---
    clusters = struct('idx', {}, 'sumT', {}, 'p', {});
    sigMask = false(1,T);
    pvals = [];

    if ~isempty(obsSegs)
        for k = 1:size(obsSegs,1)
            sc = obsScores(k);
            p = (sum(maxNull >= sc) + 1) / (nPerm + 1);
            clusters(k).idx = [obsSegs(k,1) obsSegs(k,2)];
            clusters(k).sumT = sc;
            clusters(k).p = p;
            pvals(k,1) = p;
            if p < alpha
                sigMask(1, obsSegs(k,1):obsSegs(k,2)) = true;
            end
        end
    end
end

function [t, p2, df] = welch_t(A,B)
    mA = mean(A,1,'omitnan'); vA = var(A,0,1,'omitnan'); nA = sum(isfinite(A),1);
    mB = mean(B,1,'omitnan'); vB = var(B,0,1,'omitnan'); nB = sum(isfinite(B),1);
    se = sqrt(vA./nA + vB./nB);
    t  = (mA - mB) ./ se;
    num = (vA./nA + vB./nB).^2;
    den = ((vA./nA).^2) ./ max(nA-1,1) + ((vB./nB).^2) ./ max(nB-1,1);
    df  = num ./ max(den, eps);
    p2  = 2 * (1 - tcdf(abs(t), df));
end

function clusters = label_clusters(mask)
    idx = find(mask(:)');
    clusters = struct('start',{},'end',{},'mass',{});
    if isempty(idx), return; end
    splits = [1, find(diff(idx) > 1) + 1, numel(idx)+1];
    for s = 1:numel(splits)-1
        seg = idx(splits(s):splits(s+1)-1);
        clusters(end+1).start = seg(1);
        clusters(end).end     = seg(end);
    end
end

function m = max_cluster_mass(t, supra, massfun)
    clusters = label_clusters(supra);
    if isempty(clusters)
        m = 0; return;
    end
    masses = zeros(1,numel(clusters));
    for k = 1:numel(clusters)
        seg = clusters(k).start:clusters(k).end;
        masses(k) = massfun(t(seg));
    end
    m = max(masses);
end

%%%%%%%%% ======== PLOTTING HELPERS ======== %%%%%%%%%%%

function figPairs = plot_pairwise_velocities_with_stats( ...
    data, arrayList, taxis, colorArray, dimUse, gest, ...
    usePaired, alpha_cf, alpha_cl, nPerm, minCluster, tail)

% Build T x Niter matrices per region for the chosen gesture
getVelMat = @(reg, g) ...
    (strcmp(g,'Thr')  .* squeeze(data.velThrTotal.(reg).data(:,dimUse,:)) + ...
     strcmp(g,'LS')   .* squeeze(data.velLSTotal.(reg).data(:,dimUse,:))  + ...
     strcmp(g,'Chew') .* squeeze(data.velChewTotal.(reg).data(:,dimUse,:)));

getAll = @(reg) ( squeeze(data.velThrTotal.(reg).data(:,dimUse,:)) + ...
                  squeeze(data.velLSTotal.(reg).data(:,dimUse,:)) + ...
                  squeeze(data.velChewTotal.(reg).data(:,dimUse,:)) )/3;

velMats = struct();
for ai = 1:numel(arrayList)
    reg = arrayList{ai};
    if strcmpi(gest,'All')
        velMats.(reg) = getAll(reg);           % T x N
    else
        velMats.(reg) = getVelMat(reg, gest);  % T x N
    end
end

% Region pairs
pairs = nchoosek(1:numel(arrayList), 2);

% Precompute means to set y-lims and align time base (velocities use taxis(1:end-1))

Tmin = inf;
for ai = 1:numel(arrayList)
    Tmin = min(Tmin, size(velMats.(arrayList{ai}),1));
end
for ai = 1:numel(arrayList)
    velMats.(arrayList{ai}) = velMats.(arrayList{ai})(1:Tmin, :);
end
tForVel = taxis(1:Tmin);

% Means across iterations for y-lims
allMeans = [];
for k = 1:size(pairs,1)
    A = velMats.(arrayList{pairs(k,1)});
    B = velMats.(arrayList{pairs(k,2)});
    allMeans = [allMeans; mean(A,2)'; mean(B,2)']; 
end
yMin = min(allMeans, [], 'all'); if ~isfinite(yMin), yMin = 0; end
yMax = max(allMeans, [], 'all'); if ~isfinite(yMax), yMax = 1; end
yPad   = 0.05*(yMax - yMin + eps);
yLimLo = yMin - 0.25*yPad;  
yLimHi = yMax + 0.10*yPad;

nSub = size(pairs,1);
nCols = 2; nRows = ceil(nSub / nCols);
figPairs = figure('Position',[300 100 1200 800]);

for k = 1:nSub
    subplot(nRows, nCols, k); hold on
    i = pairs(k,1); j = pairs(k,2);
    regA = arrayList{i}; regB = arrayList{j};

    A = velMats.(regA);  % T x N_A
    B = velMats.(regB);  % T x N_B

    % Mean and SEM across iterations
    muA = mean(A, 2);       seA = 2*(std(A, 0, 2) ./ sqrt(size(A,2)));
    muB = mean(B, 2);       seB = 2*(std(B, 0, 2) ./ sqrt(size(B,2)));

    % Error bars 
    hA = errorbar(tForVel, muA, 2*seA, 'LineWidth', 2, 'Color', colorArray(i,:), 'CapSize', 0);
    hB = errorbar(tForVel, muB, 2*seB, 'LineWidth', 2, 'Color', colorArray(j,:), 'CapSize', 0);

    % ---------- Significance via cluster permutation ----------
    if usePaired
        % paired iterations: same N and matched resample indices across regions
        [clusters, sigMask, ~, ~] = paired_cluster_perm_1d(A, B, alpha_cl, nPerm, minCluster);
    else
        % Welch two-sample cluster permutation (unpaired, heteroscedastic)
        % cluster_perm_2group expects nSamples x Time
        outWelch = cluster_perm_2group(A', B', 'tail', tail, 'alpha_cf', alpha_cf, 'nPerm', nPerm);
        sigMask  = outWelch.sig_mask;
        T_AB = summarize_cluster_perm(taxis(1:size(A,1)), regA, regB, outWelch, false, A, B, 0.05);
    end

    % Draw a horizontal black bar at the bottom for significant bins
    if any(sigMask)
        segs = local_mask2segments(sigMask);
        yBar = yLimLo + 0.02*(yLimHi - yLimLo);  % near the bottom
        for s = 1:size(segs,1)
            idx1 = segs(s,1); idx2 = segs(s,2);
            plot(tForVel(idx1:idx2), yBar*ones(1, idx2-idx1+1), 'k-', 'LineWidth', 6);
        end
    end

    % Axes / labels
    xlim([-0.8 0.8]);                    
    ylim([yLimLo yLimHi]);
    xlabel('Time (s)'); ylabel('Speed (a.u.)');
    title(sprintf('%s vs %s (%s, %dD)', regA, regB, gest, dimUse));
    legend([hA hB], {regA, regB}, 'Location','best'); box off
end

sgtitle(sprintf('Pairwise Neural Trajectory Speeds — %s, Dim %d', gest, dimUse), ...
        'FontSize', 20, 'FontWeight', 'bold');

end % subfunction

function segs = mask2segments(mask)
    mask = mask(:)'; 
    if ~any(mask), segs = zeros(0,2); return; end
    d = diff([false mask false]);
    segs = [find(d==1)' find(d==-1)' - 1];
end

function maskOut = prune_short_segments(mask, minLen)
    segs = mask2segments(mask);
    maskOut = false(size(mask));
    for k = 1:size(segs,1)
        if (segs(k,2) - segs(k,1) + 1) >= minLen
            maskOut(segs(k,1):segs(k,2)) = true;
        end
    end
end

function segs = local_mask2segments(mask)
    mask = logical(mask(:))';
    if ~any(mask), segs = zeros(0,2); return; end
    d = diff([false mask false]);
    starts = find(d==1);
    ends   = find(d==-1) - 1;
    segs   = [starts(:) ends(:)];
end

%%%%%%%%%% ========= STATS REPORTING HELPERS ====== %%%%%%
function T = summarize_cluster_perm(taxis, regA, regB, result, isPaired, VA, VB, pmax)
% taxis: 1xT time vector (s)
% regA, regB: region names (strings)
% result: either out (Welch) or [clusters,sigMask,pvals,tvals] packed into a struct
% isPaired: true if 'result' came from paired test
% VA, VB: T x N matrices used in the test (for direction/mean diff)
% pmax: report clusters with p <= pmax (default 0.05)

    if nargin < 7 || isempty(pmax), pmax = 0.05; end

    % ---- Normalize to common fields: t (1xT), clusters with start/end/p, sig_mask (1xT) ----
    if ~isPaired
        % Welch path: you already have these fields
        tvec      = result.t(:)';                   % 1 x Tkept (NaN elsewhere)
        sig_mask  = logical(result.sig_mask(:)');
        keptMask  = logical(result.info.keptMask(:)');
        clustersW = result.clusters;                % fields: start, end, p
        % Expand to full time base
        tFull = nan(1, numel(keptMask));
        tFull(keptMask) = tvec;
        sigFull = false(1, numel(keptMask));
        sigFull(keptMask) = sig_mask;
        clusters = arrayfun(@(c) struct('start', c.start, 'end', c.end, 'p', c.p), clustersW);
        tvals = tFull;  % use Welch t’s in full time base
        kept  = keptMask;

    else
        % Paired path: result is {clusters, sigMask, pvals, tvals}; pack to common
        % Expect fields: result.clusters(k).idx ( [i1 i2] ), and pvals
        % Build a keptMask = true(1,T)
        kept = true(1, numel(result.tvals));
        tvals = result.tvals(:)';    % 1 x T
        sigFull = logical(result.sigMask(:))';
        clustersP = result.clusters; % fields: idx=[start end], p in pvals(k)
        clusters = repmat(struct('start',[], 'end',[], 'p',[]), numel(clustersP),1);
        for k=1:numel(clustersP)
            clusters(k).start = clustersP(k).idx(1);
            clusters(k).end   = clustersP(k).idx(2);
            clusters(k).p     = result.pvals(k);
        end
        tFull = tvals; % already full length
    end

    % Final time vector aligned to what the test used
    t = taxis(kept);

    % ---- Build table, compute direction and effect descriptors ----
    rows = {};
    for k = 1:numel(clusters)
        if isempty(clusters(k).start), continue; end
        if clusters(k).p > pmax, continue; end
        idx1 = clusters(k).start; idx2 = clusters(k).end;
        t0 = t(idx1); t1 = t(idx2);
        % Winner by mean(VA-VB) within the cluster
        d  = mean(VA(idx1:idx2,:),2) - mean(VB(idx1:idx2,:),2);
        dbar = mean(d,'omitnan');
        winner = regA; if dbar < 0, winner = regB; end
        % Peak t in cluster
        [~,mx] = max(abs(tvals(idx1:idx2)));
        peak_t = tvals(idx1+mx-1);

        rows(end+1,:) = {regA, regB, round(1000*t0), round(1000*t1), clusters(k).p, string(winner), peak_t, dbar}; %#ok<AGROW>
        fprintf('%s vs %s — %d–%d ms, p_cluster=%.3g; winner=%s; peak t=%.2f; Δ=%.3g a.u.\n', ...
            regA, regB, round(1000*t0), round(1000*t1), clusters(k).p, winner, peak_t, dbar);
    end

    if isempty(rows)
        T = table([],[],[],[],[],[],[],[], 'VariableNames', ...
            {'RegionA','RegionB','tStart_ms','tEnd_ms','p_cluster','Winner','Peak_t','MeanDiff'});
        return;
    end

    T = cell2table(rows, 'VariableNames', ...
        {'RegionA','RegionB','tStart_ms','tEnd_ms','p_cluster','Winner','Peak_t','MeanDiff'});
    T = sortrows(T, 'p_cluster');
end
