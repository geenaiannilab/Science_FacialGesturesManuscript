% Preserved = 
% 38-1 (10) and 29-1 (9)

% Non-presreved = 
% 29-1 (9) and 56 -1  (18)

% Random = 38-2 (11) and 56-2 (19) 

clear all;
data = load('~/Desktop/Fig4E.mat');

cellPairs = {[10 9] [11 19] [9 18] };
bhvs2plot = [1 2];

newPixR = data.newPixR;
newPixC = data.newPixC;
taxis = data.corrMatrices.taxis;

%% for comparsion, another cell pair correlation metric -- 
% for each cell pair, for each behavior, calculate the 2D spatial
% correlation between WHOLE FRAME, BINARIZED STAs (ignoring magnitude of STA)

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
    
    % extract the signalCorr & normalized PETHs of all cellpairs in each behavior
    for bhv = 1:length(bhvs2plot)

        allCorr(pp,bhv) = data.corrMatrices.signalCorrOverTime(bhv).corr(cc,dd);
        
    end

end

counter = 1;
for  pp = 1:length(cellPairs)

    cc = cellPairs{pp}(1);
    dd = cellPairs{pp}(2);
    
    for bhv = bhvs2plot
        normTrialAvg(bhv).data(:,counter) = zscore(data.corrMatrices.trialAvg(bhv).data(:,cc));
        normTrialAvg(bhv).data(:,counter+1) = zscore(data.corrMatrices.trialAvg(bhv).data(:,dd));
        
    end
    counter = counter + 2;
end

save('matfiles/Fig4E_STAsOnly.mat', 'dataThrOut', 'dataLSOut', 'allCorr', 'normTrialAvg', 'cellPairs', 'bhvs2plot','taxis');