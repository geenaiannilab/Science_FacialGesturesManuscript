function [MI_corr, MI_raw, MI_shuff, pvals] = computeMI(spikeCounts, labels, nShuffles)
% spikeCounts: nTrials × nCells
% labels:      nTrials × 1 (categorical or numeric)
% nShuffles:   number of shuffles for bias correction
%
% Returns:
%   MI_corr   : 1 × nCells bias-corrected MI
%   MI_raw    : 1 × nCells observed MI
%   MI_shuff  : nShuffles × nCells shuffle distribution
%   pvals     : 1 × nCells empirical p-values

    nCells = size(spikeCounts,2);
    MI_raw   = zeros(1,nCells);
    MI_corr  = zeros(1,nCells);
    MI_shuff = zeros(nShuffles,nCells);
    pvals    = zeros(1,nCells);

    % --- observed MI for each cell ---
    for c = 1:nCells
        MI_raw(c) = computeMI_single(spikeCounts(:,c), labels);
    end

    % --- shuffle distribution ---
    for s = 1:nShuffles
        shuffled_labels = labels(randperm(length(labels)));
        for c = 1:nCells
            MI_shuff(s,c) = computeMI_single(spikeCounts(:,c), shuffled_labels);
        end
    end

    % --- bias-corrected MI & p-values ---
    MI_corr = MI_raw - mean(MI_shuff,1);

    for c = 1:nCells
        pvals(c) = (1 + sum(MI_shuff(:,c) >= MI_raw(c))) / (1 + nShuffles);
    end
end


function mi = computeMI_single(resp, labels)
% Computes MI between one neuron's responses and category labels

    % Make a joint histogram
    [~,~,labelBins] = unique(labels); 
    nLabels = max(labelBins);

    % response histogram (discretize spike counts as states)
    respVals = unique(resp);
    nResp = numel(respVals);
    respMap = zeros(size(resp));
    for i = 1:nResp
        respMap(resp==respVals(i)) = i;
    end

    joint = accumarray([respMap, labelBins], 1, [nResp, nLabels]);
    joint = joint / sum(joint(:));

    p_resp  = sum(joint,2);
    p_label = sum(joint,1);

    % mutual information
    mi = 0;
    for i = 1:nResp
        for j = 1:nLabels
            if joint(i,j) > 0
                mi = mi + joint(i,j) * log2(joint(i,j)/(p_resp(i)*p_label(j)));
            end
        end
    end
end

