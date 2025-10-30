function dat = simulate_timeResolved_MI_dataset()
% Simulate dataset for time-resolved MI validation
% Format matches your pipeline:
%   dat.allSpikesAllSubs.binnedSpikes : [nTrials x nTime x nCells]
%   dat.allTrials.trialType           : [nTrials x 1]  (gesture labels; non-contiguous: 2,4,7)
%   dat.cortRegion                    : {1 x nCells}   (region labels per cell)
%   dat.taxis2take                    : [1 x nTime]    (seconds)

% Rate offsets (R1_offsets, R2_peakOffsets) control how separable the gestures are; increase/decrease to make MI larger/smaller.
% 
% R2_sigma sets the temporal width of the Region-2 transient around t0 = 1.0 s.
% 
% Heterogeneity (baseRate_sd, gain_sd) keeps the population realistic (not all cells identical).
% 
% binSize = 10 ms keeps window math simple for your 400 ms window (40 bins) and 10 ms step (1 bin).
% rng(7);

%% ===================== PARAMETERS =====================
nRegions      = 3;
cellsPerReg   = 50;                 % 50 neurons / region
nCells        = nRegions * cellsPerReg;

binSize       = 0.010;              % 10 ms bins (sec)
T             = 2.0;                % total duration (sec)
nTime         = round(T / binSize); % number of time bins
t             = (0:nTime-1) * binSize;  % time axis (sec)
t0            = 1.0;                % "time zero" at the middle (sec) for Region 2

nGest         = 3;                  % three gestures
gestLabels    = [2 4 7];            % non-contiguous labels
trialsPerGest = 40;                 % trials per gesture (adjust freely)
nTrials       = nGest * trialsPerGest;

% Firing-rate baseline & modulation (Hz)
baseRate_mean = 5;                  % baseline firing (all regions)
baseRate_sd   = 1;                  % per-neuron baseline heterogeneity

% Region 1: constant separation over entire time (in Hz)
R1_offsets = [0 6 12];              % gesture-specific constant separation

% Region 2: transient separation only around t0
R2_peakOffsets = [0 8 16];          % peak separation between gestures (Hz)
R2_sigma       = 0.15;              % temporal width (sec) of the transient bump

% Region 3: no separation (same rates across gestures)
R3_offsets = [0 0 0];

% Additional per-neuron gain variability
gain_mean  = 1.0;
gain_sd    = 0.15;

% Optional per-trial multiplicative noise (kept modest)
trialGain_mean = 1.0;
trialGain_sd   = 0.05;

%% ===================== ALLOCATE OUTPUTS =====================
allSpikes = zeros(nTrials, nTime, nCells, 'uint16');  % store counts safely
allLabels = zeros(nTrials, 1);
cortRegion = cell(1, nCells);

% Build labels (balanced per gesture), then shuffle trial order
labels_vec = repelem(gestLabels, trialsPerGest)';
perm = randperm(numel(labels_vec));
labels_shuf = labels_vec(perm);

%% ===================== REGION RATE TEMPLATES =====================
% Region 1 templates (time x gesture): constant offsets
% base(t) ~ baseRate_mean (later individualized per neuron), plus gesture offsets
R1_template = zeros(nTime, nGest);
for g = 1:nGest
    R1_template(:, g) = baseRate_mean + R1_offsets(g);
end

% Region 2 templates: Gaussian bump centered at t0, only there gestures differ
gauss = exp(-0.5 * ((t - t0) / R2_sigma).^2);
R2_template = zeros(nTime, nGest);
for g = 1:nGest
    R2_template(:, g) = baseRate_mean + R2_peakOffsets(g) * gauss;
end

% Region 3 templates: equal rates (no separation)
R3_template = zeros(nTime, nGest);
for g = 1:nGest
    R3_template(:, g) = baseRate_mean + R3_offsets(g);
end

% Put into a cell for easy indexing by region
regionTemplates = {R1_template, R2_template, R3_template};
regionNames     = {'RegionEarly','RegionMid','RegionNone'};

%% ===================== SIMULATE ALL NEURONS =====================
cellIdx = 0;
for r = 1:nRegions
    tmpl = regionTemplates{r};
    for c = 1:cellsPerReg
        cellIdx = cellIdx + 1;

        % Per-neuron heterogeneity
        base_jitter = max(0.5, baseRate_mean + baseRate_sd * randn());  % keep positive
        gain        = max(0.2, gain_mean + gain_sd * randn());          % positive multiplicative
        % Build neuron-specific rate templates by shifting baseline and scaling
        rate_by_gesture = tmpl - baseRate_mean;       % remove global mean baseline
        rate_by_gesture = (rate_by_gesture .* gain) + base_jitter; % add neuron baseline back

        % Simulate all trials for this neuron
        for tr = 1:nTrials
            gLabel = labels_shuf(tr);
            gIdx   = find(gestLabels == gLabel, 1, 'first');

            % Per-trial multiplicative noise
            trialGain = max(0.2, trialGain_mean + trialGain_sd * randn());

            % Time-varying rate (Hz)
            lambda_t = trialGain * rate_by_gesture(:, gIdx)';  % 1 x nTime

            % Convert to counts per bin (Poisson with mean = rate*binSize)
            counts_t = poissrnd(max(lambda_t,0) * binSize);
            allSpikes(tr, :, cellIdx) = uint16(counts_t);
        end

        % Region label
        cortRegion{cellIdx} = regionNames{r};
    end
end

%% ===================== PACK INTO 'dat' STRUCT =====================
dat = struct();
dat.allSpikesAllSubs = struct();
dat.allSpikesAllSubs.binnedSpikes = allSpikes;               % [nTrials x nTime x nCells]
dat.allTrials = struct('trialType', labels_shuf);             % [nTrials x 1] (2,4,7, shuffled)
dat.cortRegion = cortRegion;                                  % {1 x nCells}
dat.taxis2take = t;                                           % [1 x nTime] seconds

% (Optional convenience fields)
dat.info = struct( ...
    'binSize_sec', binSize, ...
    'regions', {regionNames}, ...
    'cellsPerRegion', cellsPerReg, ...
    'gestLabels', gestLabels, ...
    'trialsPerGest', trialsPerGest, ...
    'notes', 'RegionEarly separates always; RegionMid separates only near t0; RegionNone never separates.' ...
);

fprintf('Simulated %d trials, %d time bins (%.0f ms), %d cells (%d per region).\n', ...
    nTrials, nTime, binSize*1000, nCells, cellsPerReg);

end
