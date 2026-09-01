%% Part 1
clear; clc; close all;
image_t1 = imread("t1.jpg");
figure('Name', 'T1 Visualization', 'NumberTitle', 'off');
imshow(image_t1);
title('Base Medical Image (T1)');

row_128_data = image_t1(128, :, 1);
fft_1d_result = fft(row_128_data);

figure('Name', '1D FFT Characteristics', 'NumberTitle', 'off');
subplot(2, 1, 1);
plot(20 * log(abs(fft_1d_result)), 'LineWidth', 1.2, 'Color', 'b'); 
title('Magnitude Spectrum');
ylabel('Magnitude(log)');
xlabel('Frequency(Hz)');
grid on; axis tight;

subplot(2, 1, 2);
plot(angle(fft_1d_result), 'Color', 'r'); 
title('Phase Angle');
ylabel('Phase(radians)');
xlabel('Frequency(Hz)');
grid on; axis tight;

t1_double_precision = double(image_t1(:, :, 1));
two_dimensional_fft = fft2(t1_double_precision);
shifted_fft2 = fftshift(two_dimensional_fft);

figure('Name', '2D Frequency Domain', 'NumberTitle', 'off');
subplot(1, 2, 1);
imshow(log(abs(shifted_fft2)), []);
title('2D FFT (Log Magnitude)');
subplot(1, 2, 2);
imshow(image_t1);
title('Spatial Domain Input');

%% Part 2
mask_size = 256;

circle_mask = zeros(mask_size);
[X_coords, Y_coords] = meshgrid(-128:127, -128:127);
radius_condition = sqrt(X_coords.^2 + Y_coords.^2) < 15;
circle_mask(radius_condition) = 1;

impulse_matrix = zeros(mask_size);
impulse_matrix(100, 50) = 1;
impulse_matrix(120, 48) = 2;

freq_kernel = fft2(circle_mask);
freq_impulses = fft2(impulse_matrix);
conv_freq_domain = freq_kernel .* freq_impulses;
spatial_convolution = ifft2(conv_freq_domain);

figure('Name', 'Impulse Convolution Process', 'NumberTitle', 'off');
subplot(1, 3, 1);
imshow(circle_mask, []);
title('Circular Kernel');
subplot(1, 3, 2);
imshow(impulse_matrix, []);
title('Impulse Signals');
subplot(1, 3, 3);
imshow(spatial_convolution, []);
title('Kernel * Impulses Result');

pd_scan = imread('pd.jpg');
pd_channel1 = double(pd_scan(:, :, 1));

freq_pd_scan = fft2(pd_channel1);
blurred_freq_result = freq_pd_scan .* freq_kernel;
blurred_spatial = ifft2(blurred_freq_result);

figure('Name', 'PD Image Spatial Blurring', 'NumberTitle', 'off');
subplot(1, 2, 1);
imshow(pd_channel1, []);
title('Proton Density (PD) Scan');
subplot(1, 2, 2);
imshow(fftshift(abs(blurred_spatial)), []);
title('Filtered / Blurred PD Scan');

%% Part 3
ct_scan_input = imread("ct.jpg");
ct_base = double(ct_scan_input(:, :, 1));

freq_domain_ct = fftshift(fft2(ct_base));
[r_count, c_count] = size(ct_base);
scale_multiplier = 2;

padded_rows = scale_multiplier * r_count;
padded_cols = scale_multiplier * c_count;
padded_spectrum = zeros(padded_rows, padded_cols);

idx_r_start = floor((padded_rows - r_count) / 2) + 1;
idx_c_start = floor((padded_cols - c_count) / 2) + 1;
idx_r_end = idx_r_start + r_count - 1;
idx_c_end = idx_c_start + c_count - 1;

padded_spectrum(idx_r_start:idx_r_end, idx_c_start:idx_c_end) = freq_domain_ct;

unshifted_padded = ifftshift(padded_spectrum);
interpolated_spatial = real(ifft2(unshifted_padded));
normalized_zoomed = mat2gray(interpolated_spatial);

center_cropped = normalized_zoomed(idx_r_start:idx_r_end, idx_c_start:idx_c_end);

figure('Name', 'Frequency Interpolation via Zero Padding', 'NumberTitle', 'off', 'Position', [100, 100, 1200, 400]);
subplot(1, 3, 1);
imshow(ct_base, []);
title('Input CT Scan (256x256)');

subplot(1, 3, 2);
imshow(normalized_zoomed, []);
title(sprintf('Interpolated Result (%dx%d)', padded_rows, padded_cols));

subplot(1, 3, 3);
imshow(center_cropped, []);
title('Cropped High-Res Center Region');