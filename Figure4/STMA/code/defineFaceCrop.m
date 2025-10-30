function[rect] = defineFaceCrop(videoInObj,downsampleScale)

% First do the cropping
% read the source video frame into memory
cropFrame = readFrame(videoInObj);
cropFrameGray = rgb2gray (cropFrame);
cropFrameGraySmall = imresize(cropFrameGray, (1/downsampleScale));    

hFig = figure( 'Visible', 'on');
hAxes = gca();
hCrop = imshow(cropFrameGraySmall,'Border', 'tight','Parent',hAxes);


% Get the frame
drawnow;
thisFrame=getframe(hAxes);
frame2crop = rgb2gray(thisFrame.cdata); %converts RGB to a greyscale intensity image 

% do the cropping
r1 = imrect;
position = wait(r1);
[~, rect] = imcrop(frame2crop);

close all; 
end