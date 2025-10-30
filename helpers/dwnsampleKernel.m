function dwnsmpKernel = dwnsampleKernel(fullKernel,win)

% fullKernel = L-point filter you generated (output of genSmWin.m or
%               gausswin() or something similar 
% win = bin size in SECONDS, you used to bin spikes (dwnsmple factor
%           needed) 
%
% this will downsample your smoothing kernel (normally a gaussian or causal
% half gaussian) to same resolution as binned spikes, in order to do proper
% convolution

% note: alternatively, bin your spikes at 1 ms & then no need to downsample
% kernel 

originalTimebase = 1:length(fullKernel);
newtimebase=1:(win*1e3):length(fullKernel);
dwnsmpKernel=interp1(originalTimebase,fullKernel,newtimebase);

end