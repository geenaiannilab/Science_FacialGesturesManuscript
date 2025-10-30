function [custom_colormap] = continuousColorMap(dtpts)

% Define the colormap
pad = 10;
nDatapts = dtpts+pad;
nBhvs = 3;
custom_colormap = zeros(nDatapts, 3,nBhvs);

% Define the color values 
whiteColor = [1, 1, 1];
redColor = [1,0,0]; blueColor = [0, 0, 1];greenColor = [0,1,0];
bhvColors = cat(1,redColor, blueColor, greenColor);

% Interpolate to create the colormap
for b = 1:nBhvs
    for i = 1:nDatapts
        custom_colormap(i, :,b) = (1 - (i-1) / (nDatapts-1)) * whiteColor + ((i-1) / (nDatapts-1)) * bhvColors(b,:);
    end
end

custom_colormap = custom_colormap(pad+1:end,:,:);

end
