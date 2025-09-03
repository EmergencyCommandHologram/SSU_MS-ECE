%% Parameters
theta_min   = -50;          % min angle (degrees)
theta_max   =  50;          % max angle (degrees)
N_angles    = 10;           % number of steps
f_c         = 24e9;         % carrier freq (Hz)
c           = 3e8;          % speed of light
d           = 0.5;          % element spacing in lambda (λ/2)
PHASE_WIDTH = 16;           % DDS phase input width (bits)

%% Compute
lambda = c / f_c;
theta  = linspace(theta_min, theta_max, N_angles) * pi/180; % radians
phi    = -2*pi*d*sin(theta);     % radians (relative phase shift)

%% Normalize to [0,1)
phi_norm = mod(phi, 2*pi) / (2*pi);

%% Optional: quantize to DDS word (if needed later)
phi_quant = round(phi_norm * (2^PHASE_WIDTH - 1));

%% Display
disp('Phase shifts (radians):')
disp(phi)

disp('Phase shifts normalized [0,1):')
disp(phi_norm)

%% Frequency Calculator
function [fcw, binstr] = make_fcw(f_out, f_clk, N)
    % make_fcw: compute DDS phase increment word and formatted binary string
    %   f_out: desired output freq (Hz)
    %   f_clk: DDS clock freq (Hz)
    %   N:     phase accumulator width (bits)
    %
    % returns:
    %   fcw    : integer phase increment word
    %   binstr : N-bit binary string grouped every 4 bits

    fcw = round((f_out / f_clk) * 2^N);
    % Saturate just in case
    fcw = min(max(fcw,0), 2^N-1);

    b = dec2bin(fcw, N);
    binstr = b;

    % compute actual produced freq
    f_prod = fcw / (2^N) * f_clk;
    fprintf('Binary (%d-bit): %s\n', N, binstr);
    fprintf('Produced f_out = %.9f Hz  (error = %.6f Hz)\n', f_prod, f_prod - f_out);
end

%% sweep times for coutnere
%For 100ms, 10,000,000, for 10ms, 1,000,000
%for 1ms 100,000, for 100us, 10,000
SweepTime= 100000;
AngleTime=SweepTime/10;


%% get frequencies
%make_fcw(freq, clock, bitsize)
[fcw24, b24] = make_fcw(1e6, 100e6, 24);
disp(b24)
% 2.1, 1.63, 1.13, 0.56, -0.063, -0.69, -1.25, -1.88 ,-2.38, -3.14
%330ns, 260ns, 180ns, 90ns, -10ns, -110ns, -200ns, -300ns, -380ns, -500ns
