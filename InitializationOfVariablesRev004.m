% Author: Gabriel Nicholson
% Date: -/-/2025
% Contributors: Chat GPT
%
%% Parameters
theta_min   = -50;          % min angle (degrees)
theta_max   =  50;          % max angle (degrees)
N_angles    = 10;           % number of steps
f_c         = 24e9;         % carrier freq (Hz)
c           = 3e8;          % speed of light
d           = 0.5;          % element spacing in lambda (λ/2)
PHASE_WIDTH = 24;           % DDS phase input width (bits)
NFilterOrder           = 24;  
SystemFs = 192e6;
% Bits

%% Compute
lambda = c / f_c;
theta = linspace(theta_min, theta_max, N_angles+1)*pi/180; % radians
phi = -2*pi*d*sin(theta);     % radians (relative phase shift)

%% Normalize to [0,1)
phi_norm = mod(phi, 2*pi)/(2*pi);

disp('Phase shifts normalized [0,1):')
disp(phi_norm)

%% Adjustable PARAMETERS
%/////////////////////////////////////////////////////////////////////////
%200,000 = 1ms, 20,000 = 100us, 2,000 = 10us, 200 = 1us
%/////////////////////////////////////////////////////////////////////////
TimeValue = 16384;
ChirpNc = 128; %How many chirps per frame?
AngleTime = TimeValue*ChirpNc; %10 chirps for each frame
SweepTime= TimeValue*ChirpNc*11; %2.2ms
P1 = 1/SystemFs;
P2 = 1/(SystemFs/16);

ChirpMaxFreq = 9.5e6;
ChirpMinFreq = 0.5e6;
ChirpFreqStepSize = 0.05e6;
DopplerShiftVmin = 0;
DopplerShiftVmax = 36;
DopplerStep = 3;
RXFFTLength= 16384/16;
WindowLength = 16384;
%ChirpSizeLength = FreqDeltaPeriod*361;
%Adjust OffsetDelayR based on delay between beginning, and the ASR
OffsetDelayR = 1;

%Range param for button sim
R_min = 0;
R_max = 300;
N_pos = 11;
Btn_Debounce = 20000;

%/////////////////////////////////////////////////////////////////////////
%/////////////////////////////////////////////////////////////////////////
%% Generate frequencies
% desired chirp samples
desiredLen = 16384;

% choose Nfreq as a power of two
Nfreq = 4096;

% build frequency vector with Nfreq uniformly spaced points
FreqArray = linspace(ChirpMinFreq, ChirpMaxFreq, Nfreq);

% forward sweep
FCW_p1 = FreqArray ./ SystemFs;

% backward sweep (exclude first & last to avoid duplicates)
FCW_p2 = fliplr(FCW_p1(2:end-1));

% triangular base waveform
FCW_base = [FCW_p1 FCW_p2];

% exact-length interpolation
FCW = interp1(1:length(FCW_base), FCW_base, ...
    linspace(1, length(FCW_base), desiredLen), 'linear');

% Now chirp length is exactly:
ChirpSizeLength = length(FCW);   % = 8192 exactly
FreqDeltaPeriod = 1;             % because each sample is one FCW update


%Generates 0 at either end, unlike hamming
%LengthofFCW*FreqDeltaPeriod is how many samples per chirp
WindowCoefficients = hann(WindowLength); 
WindowCoefficients = transpose(WindowLength);

% Generate velocity vector
v_mps = DopplerShiftVmin:DopplerStep:DopplerShiftVmax;

% Doppler shift in Hz
f_doppler = (2 .* v_mps) ./ lambda;    

% Normalize to sampling rate (dimensionless phase increment per sample)
DPShift = f_doppler./SystemFs;  


%% Range delays
R = linspace(R_min, R_max, N_pos);

ASR_LUT = round((R*2*SystemFs)/c);
k_test = ASR_LUT; 
R_check = k_test*c/(2*SystemFs);

%% Filter creation
%% Parameters
fc1 = 0.3e6;           % Cutoff frequency for filter 1
fc2 = 10e6;          % Cutoff frequency for filter 2
NFilterOrder = 255;             % Filter order (change if needed)
intBits = 16;        % Number of bits for integer coefficients

%% Design low-pass filters using fir1
h1 = fir1(NFilterOrder, fc1/(SystemFs/2), 'low');   % 1 MHz filter
h2 = fir1(NFilterOrder, fc2/(SystemFs/2), 'low');   % 10 MHz filter

%% Normalize to use full integer range
maxInt = 2^(intBits-1)-1;

h1_int = round(h1 / max(abs(h1)) * maxInt);
h2_int = round(h2 / max(abs(h2)) * maxInt);

filt1divide = sum(h1_int);
filt2divide = sum(h2_int);

movavgcoe = ones(1, 10)/10;
