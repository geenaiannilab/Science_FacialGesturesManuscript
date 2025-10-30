
function applyFaceCrop(videoInObj, outDir, rect, downsampleScale, nFrames)

hFig = figure( 'Visible', 'on');
hAxes = gca();

for ii = 1:nFrames
    
    
    disp(ii)
    
    % read the source video frame into memory
    sourceFrame = readFrame(videoInObj);
    sourceFrame_gray = rgb2gray (sourceFrame);
    sourceFrame_graySmall = imresize(sourceFrame_gray, (1/downsampleScale));    
    
    % show the initial frame
    if ii==1
        hImage = imshow(sourceFrame_graySmall,'Border', 'tight','Parent',hAxes);
        
    else
        hold off
        set(hImage,'CData',sourceFrame_graySmall);
    end
    hold on
    
    hPlot = gobjects(0);
        
    % Get the frame
    drawnow;
    thisFrame=getframe(hAxes);
    frame2write = rgb2gray(thisFrame.cdata); % converts RGB to a greyscale intensity image 
                                               % by eliminating hue/saturation while retaining luminance
    frame2write = imcrop(frame2write,rect) ;

    
    % Clear the plot objects
    if exist('hRender', 'var')
        delete(hRender);
    end
    delete(hPlot);
    
   % concatenate frames over time
    catFrames(:,:,ii) = frame2write; %greyscale uint8 type 
   
   
end

% Save the output 
save([outDir 'catFrames.mat'], 'catFrames', '-v7.3');
save([outDir 'cropRect.mat'],'rect');

close all;
end