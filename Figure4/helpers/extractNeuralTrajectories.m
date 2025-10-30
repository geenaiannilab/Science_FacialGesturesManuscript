function [neuralTraj,evLatentStInd] = extractNeuralTrajectories(obhv, bhvs2plot, allTestInds,LatentStTestPred,subsessionInds,tmin, tmax, bhvSR)

%%% PURPOSE: given you have run PSID_continuous_DLC.m, and 
%%%     and now want to extract neural trajectories from 
%%%     the predicted latent states
%%%%    
%%%%   This will 1) for each subsession, identify the facial expression trials 
%%%%              that were left out in the test indices
%%%%             2) extract their "new" indices, in the latentSt
%%%%             predictions (this is more complicated than it may seem, 
%%%%              since the testing set, is made up of the last 1 - trainProp% indices 
%%%%               *in each subsession* 
%%%%             3) extract the latent states, from -tmin prior to onset,
%%%%             to +tmax after onset, for each held-out trial 
%%%%            
%%%%    INPUTS: obhv = 1 x nSubs structure of blended bhvdat; *should be concatenated across
%%%%                    subsessions* 
%%%%            allTestInds = nTimepoints x 1 vector of test indices 
%%%%            LatentStTestPred = nTimepoints x nDims (output from PSID)
%%%%            subsessionInds = 1 x nSubs vector of # samples/frames in
%%%%                            each subsession 
%%%%            tmin, tmax = prior/after movement onset to examine in trajs
%%%%            bhvSR = camera rate; 70 or 120 hz (also bin window size) 

%%%%    OUTPUTS: neuralTraj = 1 x nBhvs structure; 
%%%%                         neuralTraj(bhv).data = nTimepoints x
%%%%                         nLatentDims x nTrials  

% these should be same size
latentStIndices = (1:length(LatentStTestPred))';
allTestInds;

for bhv = 1:length(bhvs2plot)

    thisBhv = bhvs2plot(bhv);

    for ss = 1:length(obhv)
        theseTrials = [obhv(ss).evScore{:,5}] == thisBhv; % trials of desired bhv
        evHandScoredInd = cell2mat(obhv(ss).evScore(theseTrials,1)); % handscored index, in subsession
        evHandScoredIndAllSubs = evHandScoredInd + sum(subsessionInds(1:ss)); % handscored index, over all subsessions
        [~,evTestingSetInd] = ismember(evHandScoredIndAllSubs,allTestInds); % index in testing set 
        evLatentStInd{ss,bhv,:} = latentStIndices(evTestingSetInd(evTestingSetInd >0));

    end
end

tmin = tmin*bhvSR;
tmax = tmax*bhvSR;

for bhv = 1:length(bhvs2plot)

    thisBhv = bhvs2plot(bhv);
    onsets = cell2mat(evLatentStInd(:,bhv));

    for tt = 1:length(onsets)
        onsetMin = onsets(tt) - tmin;
        onsetMax = onsets(tt) + tmax;
        if onsetMin <= 0 
            continue
        end
        neuralTraj(bhv).data(:,:,tt) = LatentStTestPred(onsetMin:onsetMax,:);
    end

end

%figure; plot(latentStIndices,'o'); hold on; plot( allTestInds,'o');


% example 1
% evHandScoredInd = obhv(1).evScore{end,1}; % handscored index, in 1st subsession
% evHandScoredIndAllSubs = evHandScoredInd + length(allPreviousSubs); % handscored index, over all subsessions
% evTestingSetInd = find(allTestInds == evHandScoredIndAllSubs); % index in testing set 
% evLatentStInd = latentStIndices(evTestingSetInd);
% 
% figure; 
% plot(LatentStTestPred(evLatentStInd-100:evLatentStInd+100,1),'linew',6); hold on; ...
% plot(LatentStTestPred(evLatentStInd-100:evLatentStInd+100,2),'linew',6);
% plot(LatentStTestPred(evLatentStInd-100:evLatentStInd+100,3),'linew',6); hold off
% 
% figure; 
% plot3(LatentStTestPred(evLatentStInd-100:evLatentStInd+100,1),LatentStTestPred(evLatentStInd-100:evLatentStInd+100,2),LatentStTestPred(evLatentStInd-100:evLatentStInd+100,3),'linew',8);
% view(3)
% 
% % example 2
% evHandScoredInd = obhv(2).evScore{end,1}; % handscored index, in 1st subsession
% evHandScoredIndAllSubs = evHandScoredInd + subsessionInds(1); % handscored index, over all subsessions
% evTestingSetInd = find(allTestInds == evHandScoredIndAllSubs); % index in testing set 
% evLatentStInd = latentStIndices(evTestingSetInd);
% 
% figure; 
% plot(LatentStTestPred(evLatentStInd-100:evLatentStInd+100,1),'linew',6); hold on; ...
% plot(LatentStTestPred(evLatentStInd-100:evLatentStInd+100,2),'linew',6);
% plot(LatentStTestPred(evLatentStInd-100:evLatentStInd+100,3),'linew',6); hold off
% 
% figure; 
% plot3(LatentStTestPred(evLatentStInd-100:evLatentStInd+100,1),LatentStTestPred(evLatentStInd-100:evLatentStInd+100,2),LatentStTestPred(evLatentStInd-100:evLatentStInd+100,3),'linew',8);
% view(3)

end

    