function diffFrames = ProcessPixels(catFrames)

diffFrames = diff(catFrames,1,3); 
diffFrames = abs(diffFrames);

end
