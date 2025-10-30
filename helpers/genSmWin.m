function smWin = genSmWin(smParams)

% smParams.filterType = 1,2,3 (see below)
% smParams.smoothSize = std of gaussian window in MILLISECONDS

switch smParams.filterType
    case 1 % gaussian
        gw = gausswin(round(smParams.smoothSize*6),3);
        smWin = gw./sum(gw);
        %fprintf(1, 'filter is gaussian with stdev %.2f ms\n', smParams.smoothSize);
    case 2 % half gaussian, causal
        gw = gausswin(round(smParams.smoothSize*6*3),3);
        gw(1:round(numel(gw)/2)) = 0;
        smWin = gw./sum(gw);
        %fprintf(1, 'filter is causal half-gaussian with stdev %.2f ms\n', smParams.smoothSize*3);
    case 3 % box
        smWin = ones(round(smParams.smoothSize*3),1);
        %fprintf(1, 'filter is box with width %.2f ms\n', smParams.smoothSize*3);
end
end