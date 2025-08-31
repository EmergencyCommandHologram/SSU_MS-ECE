%% Parameters
N       = 2;             % number of antennas
f_c     = 24e9;          % carrier (Hz)
c       = 3e8;           
lambda  = c / f_c;
d       = 0.5 * lambda;  % spacing = λ/2
theta0  = 30*pi/180;     % steering angle (rad)

%% Array factor computation
theta_scan = linspace(-90,90,1000)*pi/180;  % angles to plot
AF = zeros(size(theta_scan));

for k = 1:length(theta_scan)
    % phase shift relative to steering angle
    for n = 0:N-1
        AF(k) = AF(k) + exp(1j * 2*pi*d*n/lambda * ...
                            (sin(theta_scan(k)) - sin(theta0)));
    end
end

%% Normalize & plot
AF_mag = abs(AF)/max(abs(AF));
plot(theta_scan*180/pi, 20*log10(AF_mag));
xlabel('Angle (degrees)');
ylabel('Normalized gain (dB)');
grid on;
title(sprintf('Array Factor (%d elements, steer=%d°)', N, round(theta0*180/pi)));
ylim([-40 0]);
