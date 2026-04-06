clc; clear; close all;

%% Parameters
Fs = 192e6;          % Sampling rate (192 MHz)
f_sig = 300e3;       % Signal frequency (300 kHz)
duration = 200e-6;   % 200 microseconds
t = 0:1/Fs:duration-1/Fs;

%% Generate Signal
signal = sin(2*pi*f_sig*t);

% Add noise
noise_power = 0.5;
noisy_signal = signal + noise_power*randn(size(signal));

%% -------- FIR Lowpass Filter (Direct at 192 MHz) --------
fc = 1e6; % 1 MHz cutoff

fir_filtered = lowpass(noisy_signal, fc, Fs);

%% -------- CIC Decimator --------
R = 32;   % Decimation factor
N = 6;    % Number of stages
M = 1;    % Differential delay

% Create CIC decimator object
cic = dsp.CICDecimator(R, M, N);

cic_out = cic(noisy_signal.');

% Convert back to row vector
cic_out = cic_out.';

% New sampling rate after CIC
Fs_cic = Fs / R;

% Normalize CIC gain
cic_gain = (R*M)^N;
cic_out = cic_out / cic_gain;

%% -------- FIR AFTER CIC --------
fc_post = 1e6; % still 1 MHz cutoff (now relative to lower Fs)

fir_post = lowpass(cic_out, fc_post, Fs_cic);

%% -------- FFT Function --------
function fft_plot(x, Fs, titleStr)
    figure,
    Nfft = length(x);
    f = linspace(-Fs/2, Fs/2, Nfft); 
    X = fftshift(abs(fft(x))); 
    plot(f/1e6, X), grid on, 
    xlabel('Frequency (MHz)'), ylabel('Magnitude'),
    title(titleStr);
end
%% -------- Plot Time Domain --------
figure;
subplot(3,1,1)
plot(t*1e6, noisy_signal)
title('Original Noisy Signal')
xlabel('Time (us)')

subplot(3,1,2)
plot(t*1e6, fir_filtered)
title('Direct FIR @ 192 MHz')

subplot(3,1,3)
t_cic = (0:length(fir_post)-1)/Fs_cic;
plot(t_cic*1e6, fir_post)
title('CIC + FIR (Decimated)')
xlabel('Time (us)')

%% -------- Plot Frequency Domain --------
fft_plot(noisy_signal, Fs, 'Original Spectrum')
fft_plot(fir_filtered, Fs, 'Direct FIR Spectrum')
fft_plot(fir_post, Fs_cic, 'CIC + FIR Spectrum')
