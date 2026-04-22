% Author: Gabriel Nicholson
% Date: 11/11/2025
% Contributors: ChatGPT assisted
%IMPORTANT
%The order of running files is follows:
% 1. InitializationofVariables
% 2. RadarSimulinkProject (simulation is run)
% 3. This script

%% Reformat Results
signalsRange = {'RangeMAG1', 'RangeINDX1', 'RangeMAG2','RangeINDX2'};
signalsVelocity =  {'VelocityMAG1', 'VelocityINDX1', 'VelocityMAG2', 'VelocityINDX2'};
RangeSize = 128; % Chirps
VelocitySize = 2048; % Number of FFTs

for s = 1:length(signalsRange)
    name = signalsRange{s};
    data = out.(name);
    numFrames = length(data)/RangeSize;
    results.(name) = reshape(data, RangeSize, numFrames);

end

for s = 1:length(signalsVelocity)
    name = signalsVelocity{s};
    data = out.(name);
    numFrames = length(data)/VelocitySize;
    results.(name) = reshape(data, VelocitySize, numFrames);
end

if (numFrames >= 12)
    numFrames = 11;
end

%% Interpret data part 1
Nrange = 16384;   % FIXED (full FFT size)
Nvel   = 128;

fieldsRange = {'RangeINDX1','RangeINDX2'};
fieldsVel   = {'VelocityINDX1','VelocityINDX2'};

%% --- Proper index unwrapping (CRITICAL FIX) ---

% Range: convert to signed BEFORE averaging
for f = 1:numel(fieldsRange)
    idx = results.(fieldsRange{f});
    
    idx(idx >= Nrange/2) = idx(idx >= Nrange/2) - Nrange;
    
    results.(fieldsRange{f}) = idx; % KEEP signed
end

% Velocity: same idea
for f = 1:numel(fieldsVel)
    idx = results.(fieldsVel{f});
    
    idx(idx >= Nvel/2) = idx(idx >= Nvel/2) - Nvel;
    
    results.(fieldsVel{f}) = idx;
end



%% Reinterpreting as object data

object.RangeMAG = sqrt(results.RangeMAG1.^2 + results.RangeMAG2.^2); 
object.VelocityMAG = sqrt(results.VelocityMAG1.^2 + results.VelocityMAG2.^2);

% Denominator (avoid divide-by-zero properly)
denomRange = results.RangeMAG1 + results.RangeMAG2;
denomVel   = results.VelocityMAG1 + results.VelocityMAG2;

denomRange(denomRange == 0) = NaN;
denomVel(denomVel == 0) = NaN;

% Weighted averaging (SIGNED indices)
object.RangeINDX = zeros(size(results.RangeINDX1));
object.VelocityINDX = zeros(size(results.VelocityINDX1));

validRange = denomRange > 0;
validVel   = denomVel > 0;

object.RangeINDX(validRange) = ...
    (results.RangeINDX1(validRange).*results.RangeMAG1(validRange) + ...
     results.RangeINDX2(validRange).*results.RangeMAG2(validRange)) ./ denomRange(validRange);

object.VelocityINDX(validVel) = ...
    (results.VelocityINDX1(validVel).*results.VelocityMAG1(validVel) + ...
     results.VelocityINDX2(validVel).*results.VelocityMAG2(validVel)) ./ denomVel(validVel);

% Convert range to absolute AFTER averaging
object.RangeINDX = abs(object.RangeINDX);

%%Reduce values for spillover of large objects.


%% --- FALSE TARGET REJECTION ---

% Magnitude threshold (tune this)
magThresh = 0.01 * max(object.RangeMAG(:));

invalid = object.RangeMAG < magThresh;

object.RangeINDX(invalid) = NaN;
object.VelocityINDX(invalid) = 0;
object.VelocityMAG(invalid) = 0;

%% --- REMOVE VELOCITY-ONLY GARBAGE (FRAME-BASED FIX) ---

rangeThresh = 1e-6;

for j = 1:numFrames
    
    % Check if this frame has any valid range
    hasRange = any(object.RangeMAG(:,j) > rangeThresh);
    
    if ~hasRange
        % Kill entire velocity column for this frame
        object.VelocityINDX(:,j) = 0;
        object.VelocityMAG(:,j) = 0;
    end
end

% Optional: clean range NaNs → 0
object.RangeINDX(isnan(object.RangeINDX)) = 0;

%% --- SUDDEN JUMP FILTER (FIXED) ---

maxJump = 20;

for j = 2:size(object.RangeINDX,2)
    prev = object.RangeINDX(:,j-1);
    curr = object.RangeINDX(:,j);
    
    % Only consider bins with actual signal
    valid = (prev > 0) & (curr > 0);
    
    jump = abs(curr - prev);
    
    bad = (jump > maxJump) & valid;
    
    if any(bad)
        tmp = curr;
        tmp(bad) = prev(bad);
        object.RangeINDX(:,j) = tmp;
    end
end

%% Reduce data based on centeroid magnitude
thresholdFactor = 0.4; % keep strong bins only

Range = zeros(1,numFrames);
Velocity = zeros(1,numFrames);

for j = 1:numFrames
    mags = object.RangeMAG(:,j);
    vel = object.VelocityMAG(:,j);
    rangeidx  = object.RangeINDX(:,j);
    velidx  = object.VelocityINDX(:,j);

    if all(mags == 0)
        Range(j) = 0;
        continue;
    end
     if all(vel == 0)
        Velocity(j) = 0;
        continue;
    end

    Rngmask = mags > (max(mags)*thresholdFactor);
    Velmask = mags > (max(mags)*thresholdFactor);

    Range(j) = sum(rangeidx(Rngmask).*mags(Rngmask)) / sum(mags(Rngmask));
    Velocity(j) = sum(velidx(Velmask).*vel(Velmask)) / sum(vel(Velmask));
end


%% --- CONVERT TO PHYSICAL UNITS ---

Bw = ChirpMaxFreq - ChirpMinFreq;
Tchirp = ChirpSizeLength / SystemFs;
f_d = SystemFs / Nrange;

object.RangeMeters = (Range * Tchirp * f_d * c) / (2 * Bw);
object.Velocity = ( Velocity / Nvel / Tchirp) * lambda / 2;


%% Perform data compaction and prepare for plot
numFrames = size(object.RangeMeters, 2);
clustered.RangeMeters = zeros(1, numFrames);
clustered.RangeMAG = zeros(1, numFrames);

thresholdFactor = 0.4; % keep bins > 40% of the frame max

for j = 1:numFrames
    mags = object.RangeMAG(:, j);
    ranges = object.RangeMeters(:, j);

    if all(mags == 0)
        % No signal this frame
        clustered.RangeMeters(j) = 0;
        clustered.RangeMAG(j) = 0;
        continue;
    end

    % --- Step 1: threshold
    magThresh = max(mags) * thresholdFactor;
    mask = mags > magThresh;

    % --- Step 2: find contiguous nonzero clusters manually
    edges = find(diff([0; mask; 0])); % 1 where mask turns on/off
    starts = edges(1:2:end-1);
    ends = edges(2:2:end)-1;

    if isempty(starts)
        % No cluster found
        clustered.RangeMeters(j) = 0;
        clustered.RangeMAG(j) = 0;
        continue;
    end

    % --- Step 3: compute total magnitude per cluster
    numClusters = length(starts);
    clusterSum = zeros(1, numClusters);

    for k = 1:numClusters
        bins = starts(k):ends(k);
        clusterSum(k) = sum(mags(bins));
    end

    % --- Step 4: pick the strongest cluster
    [~, bestIdx] = max(clusterSum);
    bins = starts(bestIdx):ends(bestIdx);

    % --- Step 5: compute weighted centroid of range
    clustered.RangeMeters(j) = sum(ranges(bins) .* mags(bins)) / sum(mags(bins));
    clustered.RangeMAG(j) = sum(mags(bins));
end

%% Clustering
rangeRes = 16.7;

numFrames = numel(clustered.RangeMeters);
rangeVals = clustered.RangeMeters;
magVals   = clustered.RangeMAG;

% Preallocate
mergingRange = zeros(size(rangeVals));
M = zeros(size(rangeVals));       % cluster sizes (temporary)
M_full = zeros(size(rangeVals));  % same length as rangeVals

% Initialization
previous = false;
bgnidx = 1;
incr = 1;

for k = 1:numFrames-1
    diffVal = abs(rangeVals(k+1) - rangeVals(k));
    
    if diffVal >= rangeRes
        if previous
            [~, idx] = max(magVals(bgnidx:k));
            clusterRange = rangeVals(bgnidx:k);
            repIdx = bgnidx + idx - 1;

            mergingRange(repIdx) = mean(clusterRange);
            M(incr) = numel(clusterRange);
            M_full(repIdx) = M(incr);  % put value in correct position

            incr = incr + 1;
            previous = false;
        end
        bgnidx = k + 1;
        continue;
    end

    if ~previous
        previous = true;
        bgnidx = k;
    end

    if k == numFrames-1
        [~, idx] = max(magVals(bgnidx:k+1));
        clusterRange = rangeVals(bgnidx:k+1);
        repIdx = bgnidx + idx - 1;

        mergingRange(repIdx) = mean(clusterRange);
        M(incr) = numel(clusterRange);
        M_full(repIdx) = M(incr);
    end
end

% Trim M to actual clusters only
M = M(1:incr-1);


clustered.RangeMeters = mergingRange;

%% Begin visualization
% Convert your 11-frame angular positions to radians
angles_deg = linspace(-50, 50, 11);
angles_rad = deg2rad(angles_deg);

% Example data (use your processed data)
ranges = clustered.RangeMeters./2;    % 1x11 array
velocities = object.Velocity;   % 1x11 array (m/s)
magnitudes = clustered.RangeMAG;   % optional for intensity

% Create polar axes
pax = polaraxes;
hold(pax, 'on');
pax.ThetaLim = [-50 50];      % limit to semi-circle
pax.ThetaZeroLocation = 'top'; % 0° = front center
pax.ThetaDir = 'clockwise';   % right side = positive angles
pax.Color = [0.05 0.05 0.1];      % dark navy background
pax.GridColor = [0.3 0.3 0.5];    % softer gridlines
pax.GridAlpha = 0.5;              % transparency
pax.ThetaColor = 'black';         % make theta labels white
pax.RColor = 'black';             % make range labels white

% Create scaled point for size
scaleFactor = 300; 
markerSizes = scaleFactor * (M_full / max(M_full));  
markerSizes(isnan(markerSizes)) = 0;  % handle any NaN
markerSizes(markerSizes == 0) = 1;    % minimum visible point size
% Optional: color scale for velocity
scatter(pax, angles_rad, ranges, markerSizes, velocities, 'filled');

% Labeling and appearance
colormap(turbo);
cb = colorbar;
cb.Label.String = 'Velocity (m/s)';
title('Radar Polar Plot (-50° to +50°)');
