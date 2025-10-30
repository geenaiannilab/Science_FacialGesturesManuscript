function [MI_corr, MI_raw, MI_shuff] = computeMI(spikeCounts, labels, nShuffles)
% computeMI
% Compute mutual information (bits) between categorical labels and neural responses.
% 
% INPUTS:
%   spikeCounts : trials x features (can be [trials x timeBins x neurons], reshape to 2D)
%   labels      : vector of category labels (length = #trials)
%   nShuffles   : number of label permutations for shuffle baseline
%
% OUTPUTS:
%   MI_corr  : bias-corrected MI (observed - mean shuffle)
%   MI_raw   : observed MI
%   MI_shuff : [nShuffles x 1] shuffled MI values

if ndims(spikeCounts) > 2
    % Flatten timebins × neurons into "features"
    nTrials = size(spikeCounts,1);
    spikeCounts = reshape(spikeCounts, nTrials, []);
end

% Discretize spike counts (integers assumed)
counts = spikeCounts;
labels = labels(:);

% Compute observed MI
MI_raw = computeMI_discrete(counts, labels);

% Compute shuffle distribution
MI_shuff = zeros(nShuffles,1);
for s = 1:nShuffles
    shuf_labels = labels(randperm(numel(labels)));
    MI_shuff(s) = computeMI_discrete(counts, shuf_labels);
end

% Bias correction: subtract mean shuffle MI
MI_corr = MI_raw - mean(MI_shuff);
end

% -------------------------
function mi = computeMI_discrete(counts, labels)
% counts: trials x features (integer values)
% labels: categorical labels
% For multidimensional counts, treat each unique row as a response symbol.

% Convert counts to strings to treat each row as a unique "symbol"
rowStr = cellstr(num2str(counts));
[respVals,~,respIdx] = unique(rowStr);
[labelVals,~,labIdx] = unique(labels);

% Joint distribution
joint = accumarray([respIdx,labIdx],1,[numel(respVals),numel(labelVals)]);
joint = joint / sum(joint(:));

p_resp =_
