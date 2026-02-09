N = 2048 * 128;
N_range = 2048;
N_chirp = 128;

re_raw  = out.REValues1(1:N);
im_raw  = out.IMGValues1(1:N);

re_ord  = out.REValuesReordered(1:N);
im_ord  = out.IMGValuesReordered(1:N);

fft_raw = re_raw + 1j*im_raw;
fft_ord = re_ord + 1j*im_ord;

figure (1);
subplot(3,1,1)
plot(abs(fft_raw));
title('Simulink Raw FFT')

subplot(3,1,2)
plot(abs(fft_ord));
title('Simulink Reordered FFT')

fft_mat_m = zeros(262144,1);
for i = 1:2048
    for j = 1:128
        src_index  = i + (j-1)*2048;     % your original source logic
        dst_index  = (i-1)*128 + j;      % correct linear mapping
        fft_mat_m(dst_index) = fft_raw(src_index);
    end
end

subplot(3,1,3)
plot(abs(fft_mat_m));
title('Matlab Reordered FFT')
fft_ord_shift = [fft_ord(2:end); 0];

FFT_ord_mat = reshape(fft_mat_m,128,2048);
FFT_ord_simulink = reshape(fft_ord_shift,128,2048);
doppler_out_matlab = fftshift(fft(FFT_ord_mat, [], 2), 2);
doppler_out_simulink = zeros(N_range, N_chirp);
doppler_out_simulink = fftshift(fft(FFT_ord_simulink, [], 2), 2);

figure(2);
subplot(3,1,1); 
plot(real(FFT_ord_simulink(:,10)));
title('Real');
subplot(3,1,2); 
plot(imag(FFT_ord_simulink(:,10)));
title('Imag');
subplot(3,1,3); 
plot(abs(FFT_ord_simulink(:,10))); 
title('Magnitude');


figure(3);
subplot(1,2,1)
hold on;
plot(abs(fft(FFT_ord_mat(:,300))));
plot(abs(fft(FFT_ord_mat(:,600))));
plot(abs(fft(FFT_ord_mat(:,1300))));
hold off;
subplot(1,2,2)
hold on;
plot(abs(fft(FFT_ord_simulink(:,300))));
plot(abs(fft(FFT_ord_simulink(:,600))));
plot(abs(fft(FFT_ord_simulink(:,1300))));
hold off;
figure(5);
energy = sum(abs(FFT_ord_mat).^2,2);
plot(energy)
