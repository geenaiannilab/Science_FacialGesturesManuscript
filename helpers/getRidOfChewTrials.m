
function [allSpikesAllSubsOut, allTrialsOut] = getRidOfChewTrials(allSpikesAllSubsIn, allTrialsIn)

%get rid of chew trials happening inside faceExp sets (for Thor's data
%specifically) 
for i = 1:length(allTrialsIn.condition)
    if strcmp(class(allTrialsIn.condition{i}), 'char')
        allTrialsIn.condition{i} = str2double(allTrialsIn.condition{i}); 
    end
end

getRidOf = (allTrialsIn.trialType == 4 & cell2mat(allTrialsIn.condition) ~= 2);
allSpikesAllSubsIn.binnedSpikes = allSpikesAllSubsIn.binnedSpikes(~getRidOf,:,:);
allTrialsIn.frameON = allTrialsIn.frameON(~getRidOf);
allTrialsIn.trialType = allTrialsIn.trialType(~getRidOf);
allTrialsIn.condition = allTrialsIn.condition(~getRidOf);
allTrialsIn.comments = allTrialsIn.comments(~getRidOf);

allSpikesAllSubsOut = allSpikesAllSubsIn;
allTrialsOut = allTrialsIn;
end