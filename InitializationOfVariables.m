%% Parameters
theta_min   = -50;          % min angle (degrees)
theta_max   =  50;          % max angle (degrees)
N_angles    = 10;           % number of steps
f_c         = 24e9;         % carrier freq (Hz)
c           = 3e8;          % speed of light
d           = 0.5;          % element spacing in lambda (λ/2)
PHASE_WIDTH = 16;           % DDS phase input width (bits)

%% Custom Fuction for Phase Calculation
function [phi_ROM_A, phi_ROM_B, theta_rad, Phi_rad] = make_2el_phase_words(theta_min, theta_max, N_angles, ...
                                                f_c, c, d_lambda, PHASE_WIDTH)
% make_2el_phase_words  compute per-antenna phase words for a 2-element array
%   d_lambda = element spacing in wavelengths (e.g. 0.5 for lambda/2)
%   returns:
%     phi_words_A, phi_words_B : uint arrays (0..2^PHASE_WIDTH-1)
%     theta_rad : angles (radians)
%     Phi_rad   : inter-element phase difference (radians)
%
% Example:
% [A,B,th,P] = make_2el_phase_words(-50,50,10,24e9,3e8,0.5,16);

lambda = c / f_c;
theta_rad = linspace(theta_min, theta_max, N_angles) * pi/180;

% inter-element phase difference (radians)
% positive Phi corresponds to lag on element 2 relative to element 1
Phi_rad = (2*pi * d_lambda) .* sin(theta_rad);   % = 2*pi*d/λ * sin(theta)

% choose per-element offsets as +Phi/2 and -Phi/2 (centred)
phaseA =  +0.5 * Phi_rad;   % element A
phaseB =  -0.5 * Phi_rad;   % element B

% normalize to [0,1) where 1 = 2*pi
phi_ROM_A = mod(phaseA, 2*pi) / (2*pi);
phi_ROM_B = mod(phaseB, 2*pi) / (2*pi);

end

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
    % group every 4 bits for readability
    groups = regexp(b, '.{1,4}', 'match');
    binstr = strjoin(groups, ' ');

    % compute actual produced freq
    f_prod = fcw / (2^N) * f_clk;
    fprintf('N = %d, FCW = %d\n', N, fcw);
    fprintf('Binary (%d-bit): %s\n', N, binstr);
    fprintf('Produced f_out = %.9f Hz  (error = %.6f Hz)\n', f_prod, f_prod - f_out);
end

%% Get Phase values
%Params: Angle1, Angle2, Angle #, RF, c, 0.5(half space), word size
%angles,
[ROMA, ROMB, theta, Phi] = make_2el_phase_words(-50,50,10,24e9,3e8,0.5,16);
disp(table((1:10)', theta'*180/pi, Phi', ROMA', ROMB', 'VariableNames', ...
    {'Idx','Theta_deg','Phi_rad','ROMA','ROMB'}))

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