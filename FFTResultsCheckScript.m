FFTINsize = 16384*128;
re_IN  = double(out.REValues1(1:FFTINsize));
im_IN  = double(out.IMGValues1(1:FFTINsize));
x = complex(re_IN, im_IN);
x_reshaped = reshape(x, 16384,128);
X = fft(x_reshaped, 16384, 1);


re_raw = real(X);
im_raw = imag(X);




%% Code begin
N_range = 2048;
N = N_range * 128;
N_chirp = 128;

%re_raw  = out.REValues1(1:N);
%im_raw  = out.IMGValues1(1:N);



%re_ord  = out.REValuesReordered(1:N);
%im_ord  = out.IMGValuesReordered(1:N);

fft_raw = re_raw + 1j*im_raw;

fft_mat_m = zeros(262144,1);
%fft_mat_m2 = zeros(262144,1);

for i = 1:N_range
    for j = 1:128
        src_index  = i + (j-1)*N_range;     % your original source logic
        dst_index  = (i-1)*128 + j;      % correct linear mapping
        fft_mat_m(dst_index) = fft_raw(src_index);
       % fft_mat_m2(dst_index) = fft_raw(src_index);
    end
end


FFT_ord_mat = reshape(fft_mat_m,128,N_range);
%FFT_ord_simulink = reshape(fft_mat_m2,128,16384);


figure(1);
%subplot(1,2,1)
hold on;
plot(abs(fft(FFT_ord_mat(:,400))));
plot(abs(fft(FFT_ord_mat(:,800))));
plot(abs(fft(FFT_ord_mat(:,1200))));
hold off;


figure(2);
energy = sum(abs(FFT_ord_mat).^2,2);
plot(energy)

N_range = size(FFT_ord_mat,2);

peak_idx_mat = zeros(1,N_range);
peak_idx_sim = zeros(1,N_range);

for k = 1:N_range
    spec_m = fftshift(fft(FFT_ord_mat(:,k)));
    %spec_s = fftshift(fft(FFT_ord_simulink(:,k)));

    [~, peak_idx_mat(k)] = max(abs(spec_m));
   % [~, peak_idx_sim(k)] = max(abs(spec_s));
end


figure(3)

subplot(2,1,1)
plot(peak_idx_mat)
title('MATLAB Doppler peak index vs Range bin')
xlabel('Range bin')
ylabel('Doppler bin')

% subplot(2,1,2)
% plot(peak_idx_sim)
% title('Simulink Doppler peak index vs Range bin')
% xlabel('Range bin')
% ylabel('Doppler bin')
