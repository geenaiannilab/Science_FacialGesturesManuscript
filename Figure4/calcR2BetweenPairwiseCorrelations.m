% written GI 
%%%     
%%% PUPRPOSE: for each neural signal correlation matrix (ie Fig 4C), and
%%% each day, compute the R^2 between pairwise correlations
%%% (threat-lipsmack, threat-chew, lipsmack-chew); loop over all days and
%%% average 

close all; clear all

dates = {'171005','171010','171027','171128'};
subject ='Thor';
regions = {'S1','M1','PMv','M3','All'};
bhvs2plot = [1 2 4];
bhvStrings = {'Threat','Lipsmack','Chew'};

for dd = 1:length(dates)
    
    workdir = (['/Users/geena/Dropbox/PhD/SUAinfo/' subject '_' dates{dd} '/CorrMatricesV2']);
    disp(['Date is ' dates{dd}]);

    for rr = 1:length(regions)
        data = load([workdir '/' regions{rr} '.mat']);
        
        disp(['Region is ' regions{rr}])

        taxis = data.taxis;

        % plot pairwise correlations
        tmp1 = fitlm(data.prCorr{:,1}', data.prCorr{:,2}');
        r2_1v2(dd,rr) = tmp1.Rsquared.Ordinary;
        p_1v2(dd,rr) = tmp1.ModelFitVsNullModel.Pvalue;

        tmp2 = fitlm(data.prCorr{:,1}', data.prCorr{:,3}');
        r2_1v3(dd,rr) = tmp2.Rsquared.Ordinary;
        p_1v3(dd,rr) = tmp2.ModelFitVsNullModel.Pvalue;

        tmp3 = fitlm(data.prCorr{:,2}', data.prCorr{:,3}');
        r2_2v3(dd,rr) = tmp3.Rsquared.Ordinary;
        p_2v3(dd,rr) = tmp3.ModelFitVsNullModel.Pvalue;

        clear tmp1 tmp2 tmp3;
    end
end

save(['/Users/geena/Dropbox/PhD/SUAinfo/' subject '_R2ValsAcrossDays.mat'], 'dates','subject', 'regions', 'r2_1v3', 'r2_1v2', 'r2_2v3','p_1v2','p_1v3','p_2v3');

