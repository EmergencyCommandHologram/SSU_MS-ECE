%% Parameters
theta_min   = -50;          % min angle (degrees)
theta_max   =  50;          % max angle (degrees)
N_angles    = 10;           % number of steps
f_c         = 24e9;         % carrier freq (Hz)
c           = 3e8;          % speed of light
d           = 0.5;          % element spacing in lambda (λ/2)
PHASE_WIDTH = 24;           % DDS phase input width (bits)
N           = 24;           % Bits

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
FreqDeltaPeriod = 45;   %45*361 = 16245~ close to 2^14 for FFT length
TimeValue = 16384;
ChirpNc = 128; %How many chirps per frame?
AngleTime = TimeValue*ChirpNc; %10 chirps for each frame
SweepTime= TimeValue*ChirpNc*11; %2.2ms

ChirpMaxFreq = 9.5e6;
ChirpMinFreq = 0.5e6;
ChirpFreqStepSize = 0.05e6;
DopplerShiftVmin = 0;
DopplerShiftVmax = 40;
DopplerStep = 5;
SystemFs = 200e6;
RXFFTLength= 16384;
ChirpSizeLength = FreqDeltaPeriod*361;
%Adjust OffsetDelayR based on delay between beginning, and the ASR
OffsetDelayR = 1;

%Range param for button sim
R_min = 0;
R_max = 300;
N_pos = 10;
Btn_Debounce = 20000;

%/////////////////////////////////////////////////////////////////////////
%/////////////////////////////////////////////////////////////////////////
%% Generate frequencies
FreqArray = ChirpMinFreq:ChirpFreqStepSize:ChirpMaxFreq;
v_mps = DopplerShiftVmin:DopplerStep:DopplerShiftVmax;
% Normalize and stack the chrip back to back for a triangle format
FCW_p1 = (FreqArray./SystemFs);
FCW_p2 = fliplr(FCW_p1(1:end-1));
FCW = [FCW_p1, FCW_p2];

%Generates 0 at either end, unlike hamming
%LengthofFCW*FreqDeltaPeriod is how many samples per chirp
WindowCoefficients = hann(ChirpSizeLength); 
WindowCoefficients = transpose(WindowCoefficients);

%DopplerShift
DPShift = 2.*v_mps./lambda;
DPShift = (DPShift/SystemFs);
DPShift_Rom = min(max(DPShift, -2^(N-1)), 2^N-1);

% !! Temporary !! DELETE, only for testing
DPShift_Rom = fliplr(DPShift_Rom);

%% Range delays
R = linspace(R_min, R_max, N_pos);
tau = 2*R./c;                 % seconds
Rsamp = round(tau*SystemFs); % integer samples
ASR_LUT = Rsamp-OffsetDelayR;
ASR_LUT([1])=0;

%Expect freq = 
TestDelay = ((2*(60)/c)*SystemFs)

H_calc_Freq = (2*(600)*((ChirpMaxFreq-ChirpMinFreq)/100e-6))/c
H_calc_period = 1/H_calc_Freq
H_calc_Range = (c*(90)*1000)/(2*((ChirpMaxFreq-ChirpMinFreq)/81e-6))

