function correlationMatrix = correlatefROIs(varargin)

%% Input parser
p = inputParser;

% Optional params
p.addParameter('videoFileName','D:\Dropbox\FaceExpressionAnalysis\STMA2021\Barney_210704\Barney-210704-110137\cropped_Barney_FT1_-07042021110105-0000.mat',@ischar);
p.addParameter('pathToROIs','D:\Dropbox\FaceExpressionAnalysis\STMA2021\Barney_210704\Barney-210704-110137\',@ischar);
p.addParameter('statistic','mean',@ischar);
p.addParameter('spatialFilterSize',3,@isscalar);
p.addParameter('subject','Barney', @ischar);

% parse
p.parse(varargin{:})

close all; 

%% Load up some variables
subject = p.Results.subject;

% Load up the video
videoFile = load(p.Results.videoFileName);
spatialFilterSize = p.Results.spatialFilterSize;

% Load up ROIs
rois.eyes = load(fullfile(p.Results.pathToROIs, 'faceROI_eyes.mat'));
rois.LCheek = load(fullfile(p.Results.pathToROIs, 'faceROI_Lcheek.mat'));
rois.LEar = load(fullfile(p.Results.pathToROIs, 'faceROI_Lear.mat'));
rois.mouth = load(fullfile(p.Results.pathToROIs, 'faceROI_mouth.mat'));
rois.nose = load(fullfile(p.Results.pathToROIs, 'faceROI_nose.mat'));
rois.RCheek = load(fullfile(p.Results.pathToROIs, 'faceROI_Rcheek.mat'));
rois.REar = load(fullfile(p.Results.pathToROIs, 'faceROI_Rear.mat'));
rois.background.position2save = [415, 315, 581-415, 421-315];  %rectangle that doesnt overlap w/ any face ROIs

roiNames = fieldnames(rois);

%% process pixels
if strcmp(subject,'Thor')
    videoFile.catFrames = rot90(videoFile.catFrames);
    rois.background.position2save = [243, 213, 315-243, 293-13]; 
end

diffFrames = ProcessPixels(videoFile.catFrames);

%% Loop through video frames
nFrames = size(diffFrames, 3);

% Loop over frames
for ii = 1:nFrames
    
    % Smooth the frame
    diffFrames(:,:,ii) = imgaussfilt(diffFrames(:,:,ii), spatialFilterSize, 'FilterDomain', 'spatial');
    
    % Within each frame, loop over ROIs
    for rr = 1:length(roiNames)
        
        % Extract the relevant pixels
        roiPixelsDiff = imcrop(diffFrames(:,:,ii),rois.(roiNames{rr}).position2save);
        
        % Compute and stash the statistic, whether mean or variance
        if strcmp(p.Results.statistic, 'mean')
            
            result.(roiNames{rr})(ii) = mean(roiPixelsDiff(:));
            
        elseif strcmp(p.Results.statistic, 'variance')
            
            result.(roiNames{rr})(ii) = var(double(roiPixelsDiff(:)));
            
        end
        
    end
    
    
end

%% Do some plotting
c = distinguishable_colors(length(roiNames));

% Plot each time series
plotFig = figure;
nROIs = length(roiNames);

for ii = 1:nROIs
    
    subplot(nROIs/2,2, ii)
    plot(result.(roiNames{ii}),'k')
    xlabel('Frame Number')
    ylabel(p.Results.statistic)
    title(roiNames{ii},'Color',c(ii,:))
end
sgtitle('Face ROIs TimeSeries')

% Create the correlation matrix
correlationMatrix = nan(nROIs, nROIs);
plotFig = figure;
for xx = 1:nROIs
    for yy = 1:nROIs
        
        correlationMatrix(xx,yy) = corr2(result.(roiNames{xx}), result.(roiNames{yy}));
        
    end
end

imagesc(correlationMatrix); colorbar; colormap(jet)
xticklabels(roiNames);
yticklabels(roiNames);
set(gca,'xaxisLocation','top');
title('fROI vs. fROI correlation')

% Plot the ROIs 
plotFig = figure;
nROIs = length(roiNames);

imagesc(videoFile.catFrames(:,:,1)); colormap(gray); hold on
for ii = 1:nROIs
    r = rectangle('Position',rois.(roiNames{ii}).position2save,'LineWidth',4);
    r.EdgeColor = c(ii,:);
end
colororder(c);
hold off; 

end