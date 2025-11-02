% for Fig 4E, extract STAs and normalized PETH
%

clear all;
%data = load('~/Dropbox/geena/PhD/SUAInfo/Fig4E.mat');

cellPairs = {[10 9] [18 11] [17 9]};
bhvs2plot = [1 2];

newPixR = data.newPixR;
newPixC = data.newPixC;
taxis = data.corrMatrices.taxis;

%% for each cell pair, for each behavior, calculate the 2D spatial
% correlation between whole frame STAs
counter = 1;
for pp = 1:length(cellPairs)

    cc = cellPairs{pp}(1);
    dd = cellPairs{pp}(2);
    
    thisCellLS = data.lsSTAs.allCellsSTMAacrossTime.pos(1:newPixR,1:newPixC,cc); 
    thatCellLS = data.lsSTAs.allCellsSTMAacrossTime.pos(1:newPixR,1:newPixC,dd);
 
    thisCellThr = data.thrSTAs.allCellsSTMAacrossTime.pos(1:newPixR,1:newPixC,cc ); 
    thatCellThr = data.thrSTAs.allCellsSTMAacrossTime.pos(1:newPixR,1:newPixC,dd); 
   

    dataLSOut(:,:,counter) = thisCellLS;
    dataLSOut(:,:,counter + 1) = thatCellLS;

    dataThrOut(:,:,counter) = thisCellThr;
    dataThrOut(:,:,counter+1) = thatCellThr;

    counter = counter +2;

end

%% extract their signalCorrelations (as in Fig 4D)
%% extrac their PETH for plotting 
counter = 1;
for  pp = 1:length(cellPairs)

    cc = cellPairs{pp}(1);
    dd = cellPairs{pp}(2);
    
    for bhv = bhvs2plot
        normTrialAvg(bhv).data(:,counter) = zscore(data.corrMatrices.trialAvg(bhv).data(:,cc));
        normTrialAvg(bhv).data(:,counter+1) = zscore(data.corrMatrices.trialAvg(bhv).data(:,dd));
        
    end

    % extract the signalCorr all cellpairs in each behavior
    for bhv = 1:length(bhvs2plot)

        allCorr(pp,bhv) = data.corrMatrices.signalCorrOverTime(bhv).corr(cc,dd);
        
    end
    counter = counter + 2;
end

save('matfiles/Fig4E_STAsOnly.mat', 'dataThrOut', 'dataLSOut', 'allCorr', 'normTrialAvg', 'cellPairs', 'bhvs2plot','taxis');