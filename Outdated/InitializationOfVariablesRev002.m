%% Parameters
theta_min   = -50;          % min angle (degrees)
theta_max   =  50;          % max angle (degrees)
N_angles    = 10;           % number of steps
f_c         = 24e9;         % carrier freq (Hz)
c           = 3e8;          % speed of light
d           = 0.5;          % element spacing in lambda (λ/2)
PHASE_WIDTH = 24;           % DDS phase input width (bits)

%% Compute
lambda = c / f_c;
theta  = linspace(theta_min, theta_max, N_angles) * pi/180; % radians
phi    = -2*pi*d*sin(theta);     % radians (relative phase shift)

%% Normalize to [0,1)
phi_norm = mod(phi, 2*pi) / (2*pi);

disp('Phase shifts normalized [0,1):')
disp(phi_norm)

%% Adjustable PARAMETERS
%/////////////////////////////////////////////////////////////////////////
%/////////////////////////////////////////////////////////////////////////
SweepTime= 20000;
AngleTime=SweepTime/10;
ChirpSamp = 51;
ChirpMaxFreq = 19e6;
ChirpMinFreq = 1e6;
ChirpFreqStepSize = 0.5e6;
SystemFs = 200e6;

%/////////////////////////////////////////////////////////////////////////
%/////////////////////////////////////////////////////////////////////////
%% Generate frequencies
FreqArray = ChirpMinFreq:ChirpFreqStepSize:ChirpMaxFreq;
% Normalize and stack the chrip back to back for a triangle format
FCW_p1 = (FreqArray./SystemFs);
FCW_p2 = fliplr(FCW_p1(1:end-1));
FCW = [FCW_p1, FCW_p2];


