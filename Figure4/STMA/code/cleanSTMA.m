function[cleanedSTMA] = cleanSTMA(faceDataZscoredFiltThresh, varargin)

% PURPOSE: will remove pixels that are *not* above zscore threshold in >=2
%           temporally consecutive frames in STMA 

if length(varargin) ~= 0
    plotFlag = varargin{1,1};
else
    plotFlag = 0;
end

% pattern of temporally sequential pixels to find (ie., not above thres, ABOVE, not above) 
str2find = [0 1 0];

pixW = size(faceDataZscoredFiltThresh,1);
pixH = size(faceDataZscoredFiltThresh,2);
nFrames = size(faceDataZscoredFiltThresh,3);

reshapedSTMA = (reshape(faceDataZscoredFiltThresh, [pixW*pixH], nFrames));

reshapedCleanedSTMA = zeros(size(reshapedSTMA));
for i = 1:size(reshapedSTMA,1)
    currRow = reshapedSTMA(i,:) ~= 0; %frames where pixel > threshold
    ind2Remove = strfind(currRow, str2find); %frames -1; where it goes 0-1-0
    ind2Remove = ind2Remove+1;
    currRow(ind2Remove) = 0; %cleaned logical vector
    newRow = reshapedSTMA(i,:) .* currRow;
    reshapedCleanedSTMA(i,:) = newRow;
end

cleanedSTMA = reshape(reshapedCleanedSTMA, pixW, pixH, nFrames); %reshape back to 3D
cmax = max(cleanedSTMA,[],'all');
cmin = min(cleanedSTMA,[],'all');

if plotFlag
    figure;
    for j = 1:nFrames
       subplot(121);
       imagesc(faceDataZscoredFiltThresh(:,:,j)); colorbar; caxis([cmin cmax])
       subplot(122)
       imagesc(cleanedSTMA(:,:,j)); colorbar; caxis([cmin cmax])
       pause(0.2); 
    end
end
end