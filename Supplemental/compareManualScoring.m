
%%%%% PLOTS Fig. S1 
%%%%% Handscored gesture onsets compared to automatic detection by continuous facial marker tracking
%%%%%     'markerTrajectories.mat' been smoothed with smoothdata(DLCmarkersIn, 'gaussian') which is truncated 50 ms gaussian (SD is 1/5th of window), then data rebinned at 10 ms.
%%%%%%     This is the same input as Figures 1 B-F
%%%%%% 

clear all; close all;

load(['../Figure1/matfiles/markerTrajectories.mat']);

tmin = taxis(1); tmax = taxis(end); %this is in seconds
tbin=0.01;                          %bin width is 10 ms
t0=find(taxis==0);
tspon=[1 80];                       %bins of spontaenous activity period
ntspon_bins=tspon(2)-tspon(1)+1;
pre_mov_period=taxis(1:t0);

n_bhvs = length(bhvs2plot);

tc_times_std=3; %number of standard deviations for threshold crossing
sw_ROC=1; %step width for ROC analysis

cols = [1, 0 0 ;   
        0, 0, 1;   
        0, 1, 0];      

blendFactor = 0.5;  % 0 = vivid, 1 = pure white
colors = cols * (1 - blendFactor) + blendFactor * ones(3,3);

SZ=get(0, 'Screensize');
SZhw=SZ;SZhw(3)=SZhw(3)/2;

fig1 = figure('Name','Peak and Peak Slope','Position', SZhw); 
fig1_nrow=5; fig1_ncol=3;

for bhv = 1:n_bhvs
  
    %trajectories by time, trials, and PCs
    trajs = squeeze(bhvTraj(bhv).data); %time, trial, 
    [time, trials, PCnum] = size(trajs);
    
    %trajectory length by time, trials
    trajs_norm=zeros(time, trials);
    for i=1:trials
      trajs_norm(:,i)=vecnorm(squeeze(trajs(:,i,:)),2,2);
    end
    
    %(1) Basic Characterization 
    % mean and standard deviation over trials, i.e. function of time
    avg_traj=mean(trajs_norm,2);
    std_traj=zeros(time,1);
    for i=1:time
      std_traj(i)=std(trajs_norm(i,:));
    end
    
    %mean during spont activity,function of trial
    avg_spon_traj=mean(trajs_norm(tspon(1):tspon(2),:),1);
    
    %mean and standard deviation during spontaneous activity period
    avg_spon=mean(avg_traj(tspon(1):tspon(2)));
    std_spon=std(reshape(trajs_norm(tspon(1):tspon(2),:),ntspon_bins*trials,1));
    
    %plot in first row
    subplot(fig1_nrow,fig1_ncol,bhv)
    hold off ; ax = gca; ax.YAxis.Visible = 'off';
    plot(taxis, trajs_norm,'Color',colors(bhv,:));hold on;
    xlabel('Time [sec]','FontSize',12); ylabel('Amplitude [a.u.]','FontSize',12);
    plot(taxis, avg_traj  ,'LineWidth',3,'Color','black');hold on;
    line([tmin tmax],avg_spon+0*[std_spon std_spon],'Color','black');
    line([tmin tmax],avg_spon+tc_times_std*[std_spon std_spon],'Color','black');
    line([tmin tmax],avg_spon-tc_times_std*[std_spon std_spon],'Color','black');
    ax.XAxis.TickDirection = 'out'; ax.YAxis.TickDirection = 'out';
    fontsize(ax,12,'point');
    title([gestureNames{bhv}],'FontSize',14);
    
    %(3) Threshold crossing (simply the first transition)
    thr=avg_spon+tc_times_std*std_spon;
    fa_thr_cr=zeros(trials,1)-2;
    for i=1:trials
       %thr_cr=min(find(trajs_norm(:,i) > thr));
       thr_cr=find(trajs_norm(:,i) > thr,1);
       if ~isempty(thr_cr)
         fa_thr_cr(i)=thr_cr;
       else
         fa_thr_cr(i)=NaN;
       end
    end
    disp(newline);
    disp(gestureNames{bhv});
    disp(['Number of Trials: ' num2str(trials)]);
    disp(['Threshold set at ' num2str(tc_times_std) ' Standard Deviations']);
    disp(['Number of Trials without Threshold crossing: ',num2str(sum(isnan(fa_thr_cr)))]);

    fa_thr_cr=tbin*(fa_thr_cr-t0);
    subplot(fig1_nrow,fig1_ncol,2*fig1_ncol+bhv);
    h_TCR=histogram(fa_thr_cr,tmin-0.05:0.1:tmax+0.05);
    std_TCR = std(fa_thr_cr(~isnan(fa_thr_cr)));
    set(h_TCR,'FaceColor',colors(bhv,:));
    set(gca,'TickDir','out');
    xlabel('Time [sec]','FontSize',12); ylabel('# Trials','FontSize',12);
    title('Response Latency','FontSize',12);
    fontsize(gca,12,'point');

    %(5) ROC analysis
    % change threshold and determine fraction of trials with crossing prior
    % to hand scored onset (false negative) and fraction of trials with
    % crossing after hand scored (true positive)
    % True Positive Ratio (TPR) (true positive over (true positive plus false negative)
    % False Positive Ratio (FPR) (falso positive over (false positive plus true negative)
    low_thr=floor(min(trajs_norm,[],"all")-0.5*sw_ROC);
    lar_thr=ceil(max(trajs_norm,[],"all")+0.5*sw_ROC); 
    fa_thr=low_thr:sw_ROC:lar_thr; 

    n_thr=length(fa_thr);
    fa_TP=zeros(n_thr,trials);
    fa_FP=zeros(n_thr,trials);
    fa_thr_cr=zeros(trials,1)-2;
    for i=1:n_thr
      thr=fa_thr(i);
      for j=1:trials
        fp=find(trajs_norm(   1:t0-1,j) > thr,1); %false positive 
        tp=find(trajs_norm(t0+1:end ,j) > thr,1); %true positive
        fa_FP(i,j)=~isempty(fp);
        fa_TP(i,j)=~isempty(tp);
      end
    end
    fa_TPR=mean(fa_TP,2);
    fa_FPR=mean(fa_FP,2);
    
    AUC = trapz(flip(fa_FPR), flip(fa_TPR));
    disp([gestureNames{bhv} ':  AUC = ' num2str(AUC)]);
    fprintf('std_spon=%5.3f\n',std_spon);

    subplot(fig1_nrow,fig1_ncol,4*fig1_ncol+bhv);
    hold off ; ax = gca; ax.YAxis.Visible = 'off';
    plot(fa_FPR,fa_TPR,'-o','Color',colors(bhv,:)); axis square
    set(gca,'TickDir','out');
    xlabel('False Positive Ratio','FontSize',12); ylabel('True Positive Ratio','FontSize',12);
    title('ROC','FontSize',12);
    fontsize(gca,12,'point');

    %(4) Compute peak and peak slope paramters
    
    %get maximum positions:
    [mr, mi]=max(trajs_norm);
    tnmax_p=mi;
   
    tnmax_pp=tbin*(tnmax_p-t0);
    subplot(fig1_nrow,fig1_ncol,1*fig1_ncol+bhv);
    h_MAXP=histogram(tnmax_pp,tmin-0.05:0.1:tmax+0.05);
    std_MAXP = std(tnmax_pp(~isnan(tnmax_pp)));
    disp([gestureNames{bhv} ' - Max Response Peak - Std = ' num2str(std_MAXP)]);

    set(h_MAXP,'FaceColor',colors(bhv,:));
    set(gca,'TickDir','out');
    xlabel('Time [sec]','FontSize',12); ylabel('# Trials','FontSize',12);
    title('Max Response Peak','FontSize',12);
    fontsize(gca,12,'point');

    %compute differential
    tnt=cat(1,trajs_norm(1,:),trajs_norm);
    trajs_norm_diff=tnt(2:end,:)-tnt(1:end-1,:);
    avg_traj_diff=mean(trajs_norm_diff,2);

    [mrd, mi]=max(trajs_norm_diff);
    tdmax_p=mi;
   
    tdmax_pp=tbin*(tdmax_p-t0);
    subplot(fig1_nrow,fig1_ncol,3*fig1_ncol+bhv);
    h_MAXDP=histogram(tdmax_pp,tmin-0.05:0.1:tmax+0.05);
    std_MAXDP = std(tdmax_pp(~isnan(tdmax_pp)));
    disp([gestureNames{bhv} ' - Max Slope - Std = ' num2str(std_MAXDP)]);

    set(h_MAXDP,'FaceColor',colors(bhv,:));
    set(gca,'TickDir','out');
    xlabel('Time [sec]','FontSize',12); ylabel('# Trials','FontSize',12);
    title('Max Slope','FontSize',12);
    fontsize(gca,12,'point');

    %compute sliding average comparing distributions across time
    %we need multiple comparison correction
    
    ww=1;           %window width in number of bins 
    wc=ceil(ww/2.); %window center
    sw=9;           %step width from window to window
    sb=ceil(sw/2.); %start bin
    n_steps=floor(time/sw);
    tnsc=mean(trajs_norm(1:ww,:),1); %comparison distribution
    sw_taxis=taxis(sb:sw:end);
    fa_tns=zeros(n_steps,trials);
    for i=1:n_steps 
      tns=mean(trajs_norm(sb+(wc-ww:ww-wc)+(i-1)*sw,:),1);
      fa_tns(i,:)=tns;
    end
    [p_kw,tbl_kw,stats_kw] = kruskalwallis(fa_tns',[],'off');
    disp(['Kruskal Wallis Test for ' gestureNames{bhv} ': p=' num2str(p_kw,'%e')]);

    cc_mc_kw=multcompare(stats_kw,"CriticalValueType","dunnett","Approximate",1,"Display","off");
    p_mc_kw=cc_mc_kw(:,6); %now limited to the comparison to control group

    subplot(fig1_nrow,fig1_ncol,1*fig1_ncol+bhv);
    hold off ; ax = gca; ax.YAxis.Visible = 'off';
    semilogy(sw_taxis,[1 p_mc_kw'],'Color',colors(bhv,:),'LineStyle','-','LineWidth',2);
    ylim([10^-15 2]);
    line([tmin tmax],0.01*[1 1],'Color','black','LineStyle',':');
    title('Kruskal Wallis Test','FontSize',12); 
    xlabel('Time [sec]','FontSize',12); ylabel('Probability','FontSize',12);
    ax.XAxis.TickDirection = 'out'; ax.YAxis.TickDirection = 'out';
    fontsize(ax,12,'point');
end