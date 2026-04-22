% Author: Gabriel Nicholson
% Date: 4/20/2026
% Contributors: Chat GPT
%
%% Parameters
theta_min   = -50;          % min angle (degrees)
theta_max   =  50;          % max angle (degrees)
N_angles    = 10;           % number of steps
f_c         = 24e9;         % carrier freq (Hz)
c           = 3e8;          % speed of light
d           = 0.5;          % element spacing in lambda (λ/2)
SystemFs = 192e6;

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
AngleTime = TimeValue*ChirpNc; %Number of samples in each angle
SweepTime= TimeValue*ChirpNc*11; %Number of samples in whole sweep
DeciFactor = 8; %how much to decimate
P1 = 1/SystemFs; %Common period for most blocks
P2 = 1/(SystemFs/DeciFactor); %Decimated period (not used?)

ChirpMaxFreq = 11e6; 
ChirpMinFreq = 1e6;
ChirpFreqStepSize = 0.00055e6;
DopplerShiftVmin = 0;
DopplerShiftVmax = 36;
DopplerStep = 3;
RXFFTLength= TimeValue/DeciFactor;
WindowLength = TimeValue;

%Range param for button sim
R_min = 0;
R_max = 300;
N_pos = 9;
Btn_Debounce = 20000;

%/////////////////////////////////////////////////////////////////////////
%/////////////////////////////////////////////////////////////////////////
%% Generate frequencies
%% Parameters
desiredLen = 16384;  % samples per chirp
t = (0:desiredLen-1)/SystemFs;
T = desiredLen / SystemFs;

k = (ChirpMaxFreq - ChirpMinFreq) / (T/2); % slope

% Triangular chirp
f = zeros(size(t));

for i = 1:length(t)
    if t(i) <= T/2
        f(i) = ChirpMinFreq + k*t(i);
    else
        f(i) = ChirpMaxFreq - k*(t(i)-T/2);
    end
end

FCW = f / SystemFs;
%FCW = (ChirpMinFreq:ChirpFreqStepSize:ChirpMaxFreq)./SystemFs;

% Now chirp length is exactly:
ChirpSizeLength = TimeValue;   % = 8192 exactly
%Generates 0 at either end, unlike hamming
%LengthofFCW*FreqDeltaPeriod is how many samples per chirp
WindowCoefficients = hamming(ChirpSizeLength); 
WindowCoefficients = transpose(WindowCoefficients);

DopplerWindow = hamming(ChirpNc);
DopplerWindow = transpose(DopplerWindow);
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

%% Design low-pass filters using fir2
f1 = [0 0.4 0.5 1];
m1 = [1 1 0 0];
b1 = fir2(30,f1,m1);

f2 = [0 0.2 0.3 1];
m2 = [1 1 0 0];
b2 = fir2(30,f2,m2);

f3 = [0 0.15 0.15 1];
m3 = [1 1 0 0];
b3 = fir2(30,f3,m3);

figure(1);
[h1,w1] = freqz(b1,1);
[h2,w2] = freqz(b1,1);
[h3,w3] = freqz(b1,1);
subplot(1,3,1);
plot(f1,m1,w1/pi,abs(h1))
subplot(1,3,2);
plot(f2,m2,w2/pi,abs(h2))
subplot(1,3,3);
plot(f3,m3,w3/pi,abs(h3))

