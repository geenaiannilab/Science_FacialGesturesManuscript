function [gaussKern, tminLong, tmaxLong, taxisLong, taxis2take] = fixTimeAxisEdgeArtifacts(gaussKernStd, win, tmin, tmax)

smParams.smoothSize = gaussKernStd; %std in milisec
smParams.filterType = 1; % 1= gauss, 2= caus half-gauss, 3=box
smWin = genSmWin(smParams);
gaussKern = dwnsampleKernel(smWin,win);
extraBins = length(gaussKern); % extra bins to discard due to conv edges 

% original taxis
taxis = -tmin:win:tmax;

% these nasty bits remove edge artifacts from gaussian convolv
tminLong = tmin + (extraBins*win);
tmaxLong = tmax + (extraBins*win);
taxisLong = -tminLong:win:tmaxLong;

taxis2take = zeros(1,length(taxisLong));
[~, firstOverlappingIndex] = (min(abs(taxisLong - taxis(1))));
[~, lastOverlappingIndex] = (min(abs(taxisLong - taxis(end))));
taxis2take(firstOverlappingIndex:lastOverlappingIndex) = 1;
taxis2take = logical(taxis2take);

end