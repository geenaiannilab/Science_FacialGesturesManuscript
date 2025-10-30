% PURPOSE: load a single .avi file (one 10 minute session)
% PURPOSE: downsample each frame to fewer pixels, and then concatenate
            % these matrices of downsampled pixels 
% OUTPUT: catFrames = downsampled pixels Width X Height X Frames, where
%            values are uint8 intensity values 
% OUTPUT: note: this is a large, large .mat file (~600MB) that you generated once (in testing 1 take) & saved    
%%% 'Barney-210704-110137'
%%% 'Barney_FT2_-07042021110110-0000.avi'
%%% 'Barney_FT1_-07042021110105-0000.avi'
% written by GI 180821 

%%% eg -- ExtractMovieFrames2021('G:/Barney_210704/FT1/Barney_FT1_-07042021110105-0000.avi','dirOut/testOut.mat','downsampleScale','8');

function ExtractMovieFrames2021(videoInFileName, outDir, cropOptions, )
 
close all;
workDir = extractBefore(outDir,'cropped');
workDir = extractBefore(outDir,'Take');

% Prepare the video object
%downsampleScale = 1;%str2double(p.Unmatched.downsampleScale);
%videoInObj = VideoReader([videoInFileName '.avi']);
videoInObj = VideoReader([videoInFileName '.mp4']);

% get number of frames
nFrames = videoInObj.NumFrames;

% get video dimensions
videoSizeX = videoInObj.Height;
videoSizeY = videoInObj.Width;

if strcmpi(cropOptions,'defineCrop') 
    [rect] = defineFaceCrop(videoInObj,downsampleScale);
    applyFaceCrop(videoInObj, outDir, rect, downsampleScale, nFrames);
    
elseif  strcmpi(cropOptions,'applyCrop') 
    load([workDir 'cropRect.mat']);
    applyFaceCrop(videoInObj, outDir, rect, downsampleScale, nFrames);
end

% % Save the output 
% outDir = p.Results.videoOutFileName;
% save([outDir 'catFrames.mat'], 'catFrames', '-v7.3');
% save([outDir 'cropRect.mat'],'rect');
% 
% close all;

end