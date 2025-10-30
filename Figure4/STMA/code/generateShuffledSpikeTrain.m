function [randSpikeTrain] = generateShuffledSpikeTrain(spikes, varargin)
 
% Generates shuffled spike train within specified window for comparing to
 % real STMA 
 if size(varargin,2) == 1
    shuffleWindow = varargin{1};
    plotFlag = 0;
 elseif size(varargin,2) == 2
    shuffleWindow = varargin{1};
    plotFlag = varargin{2};
 else
    disp('indicate shuffle window size!')
end
% count number of real spikes 
nSpikes = length(spikes);

% generate random timestamps of spike train, same length as actual data; this neuron is ~30
% hz

randSpikeTrain = zeros(1,nSpikes);

for i = 1:nSpikes
    lowerBound = spikes(i) - shuffleWindow/2;
    upperBound = spikes(i) + shuffleWindow/2;
    randSpikeTrain(i) = (upperBound-lowerBound).*rand(1,1) +lowerBound;
end



if plotFlag
    ax1 = ones(1,nSpikes);
    figure; 
    plot(spikes, ax1, 'o');
    hold on 
    plot(randSpikeTrain, ax1, 'o');
    title('spikeTimes'); xlabel('s')
    legend('realSpikes', 'shuffledSpikes')

    bin = 0.04;
    ax = 0:bin:spikes(end)+1;
    fr_spikes = histcounts(spikes,ax);
    fr_shuffled = histcounts(randSpikeTrain, ax);

    figure; 
    plot(ax(1:end-1), fr_spikes);
    hold on 
    plot(ax(1:end-1), fr_shuffled);
    title('PSTH, real vs shuffled Spikes')
    legend('realSpikes', 'shuffledSpikes')
    hold off

    taxis = -1000:1000;  % make a time axis of 2000 ms
    t_wid = 5;  % width of kernel (in ms)
    gauss_kernel = normpdf(taxis, 0, t_wid);
    gauss_kernel = gauss_kernel ./ sum(gauss_kernel); %

    fr_spikes_sm = conv(fr_spikes, gauss_kernel, 'same')/bin; %
    fr_shuffled_sm = conv(fr_shuffled, gauss_kernel, 'same')/bin;

    figure; 
    plot(ax(1:end-1), fr_spikes_sm); hold on 
    plot(ax(1:end-1), fr_shuffled_sm); hold off
    title('SDF; real vs shuffled spikes')
    legend('realSpikes', 'shuffledSpikes')

end

randSpikeTrain = randSpikeTrain';

end