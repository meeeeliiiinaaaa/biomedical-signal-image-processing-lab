%% Q1: 
clearvars; close all; clc;

raw_data = double(imread("t2.jpg"));
base_img = raw_data(:, :, 1);
noise_level = 15;
img_corrupted = base_img + (noise_level * randn(size(base_img)));

[grid_x, grid_y] = meshgrid(-127:128, -127:128);
spatial_kernel = zeros(256, 256);
spatial_kernel((-1 <= grid_x) & (grid_x <= 2) & (-1 <= grid_y) & (grid_y <= 2)) = 1;

spatial_kernel = spatial_kernel / sum(spatial_kernel(:)); 

H_freq = fftshift(fft2(ifftshift(spatial_kernel)));
Img_freq = fftshift(fft2(ifftshift(img_corrupted)));

Filtered_freq = Img_freq .* H_freq;
restored_custom = fftshift(ifft2(ifftshift(Filtered_freq)));

figure('Name', 'Custom Kernel Smoothing');
subplot(1, 3, 1);
imshow(abs(base_img), []);
title('Source Image (Clean)');

subplot(1, 3, 2);
imshow(abs(img_corrupted), []);
title('Corrupted via Gaussian Noise');

subplot(1, 3, 3);
imshow(abs(restored_custom), []);
title('Recovered (4x4 Matrix Filter)');

gauss_smoothed_img = imgaussfilt(base_img, 1);

figure('Name', 'Filtering Methods Comparison');
subplot(1, 2, 1);
imshow(abs(restored_custom), []);
title('Result of Custom FFT Filter');

subplot(1, 2, 2);
imshow(abs(gauss_smoothed_img), []);
title('Result of imgaussfilt');

%% Q2: 

img_data2 = imread("t2.jpg");
true_signal = double(img_data2(:, :, 1));

blur_std = 0.8;
psf_kernel = Gaussian(blur_std, [256, 256]);
blurred_observation = conv2(true_signal, psf_kernel, 'same');

H_psf = fftshift(fft2(ifftshift(psf_kernel)));
G_blurred = fftshift(fft2(ifftshift(blurred_observation)));

F_estimate = G_blurred ./ H_psf;
recovered_clean_signal = fftshift(ifft2(ifftshift(F_estimate)));

figure('Name', 'Deconvolution (No Noise)');
subplot(1, 3, 1);
imshow(abs(true_signal), []);
title('Reference Ground Truth');

subplot(1, 3, 2);
imshow(abs(blurred_observation), []);
title('Degraded (Blurred) Image');

subplot(1, 3, 3);
imshow(abs(recovered_clean_signal), []);
title('Inverse Filter Output');

noisy_blurred_obs = blurred_observation + (randn(size(blurred_observation)) * 0.001);
G_noisy_obs = fftshift(fft2(ifftshift(noisy_blurred_obs)));

F_noisy_estimate = G_noisy_obs ./ H_psf;
recovered_noisy_signal = fftshift(ifft2(ifftshift(F_noisy_estimate)));

figure('Name', 'Deconvolution (With Noise Effect)');
subplot(1, 3, 1);
imshow(abs(true_signal), []);
title('Reference Ground Truth');

subplot(1, 3, 2);
imshow(abs(blurred_observation), []);
title('Degraded (Blurred) Image');

subplot(1, 3, 3);
imshow(abs(recovered_noisy_signal), []);
title('Recovery Failure (Ill-Conditioned)');

%% Q3: 

image = imread("t2.jpg");

N = 64;
image = image(:, :, 1);
f = imresize(double(image), [N N]);

K = zeros(N, N);
h = [0 1 0; 1 2 1; 0 1 0];
K(1:3, 1:3) = h;

D = zeros(N^2, N^2);

for r = 0:N-1
    for c = 0:N-1

        K_shift = circshift(K, [r, c]);

        D(r*N + c + 1, :) = K_shift(:)';
    end
end

% h * f = Df
g = D * f(:);

% add noise
sigma = 0.05;
g_noisy_vec = g + sigma * randn(size(g));
g_noisy = reshape(g_noisy_vec, N, N);

f_recovered = pinv(D) * g_noisy_vec;  
f_recovered = reshape(f_recovered, N, N);

figure('Position', [100 100 1200 400]);
subplot(1,3,1);
imshow(f, []);
title('Original image');

subplot(1,3,2);
imshow(g_noisy, []);
title('Noisy filtered image');

subplot(1,3,3);
imshow(f_recovered, []);
title('Recovered image');

%% Q4: 

image = imread("t2.jpg");

N = 64;
image = image(:, :, 1);
f = imresize(double(image), [N N]);

% create D
K = zeros(N, N);
h = [0 1 0; 1 2 1; 0 1 0];
K(1:3, 1:3) = h;

D = zeros(N^2, N^2);

for r = 0:N-1
    for c = 0:N-1

        K_shift = circshift(K, [r, c]);

        D(r*N + c + 1, :) = K_shift(:)';
    end
end

% h * f = Df
g = D * f(:);

sigma = 0.05;
g_noisy_vec = g + sigma * randn(size(g));
g_noisy = reshape(g_noisy_vec, N, N);

f_recovered = pinv(D) * g_noisy_vec;  
f_recovered = reshape(f_recovered, N, N);

betha = 0.01;
max_iter = 100;
toleration = 1e-6;

[f_recovered_grd, history] = gradient_descent(max_iter, D, g_noisy_vec, betha, toleration);
f_recovered_grd = reshape(f_recovered_grd, N, N);


figure('Position', [100 100 1200 500]);

subplot(1,4,1);
imshow(f, []);
title('Original image');

subplot(1,4,2);
imshow(g_noisy, []);
title('Noisy filtered image');

subplot(1,4,3);
imshow(f_recovered, []);
title('Recovered image using inverse');

subplot(1,4,4);
imshow(f_recovered_grd, []);
title(sprintf('Recovered image using gradient descent', length(history)));


figure;
plot(1:length(history), history, 'LineWidth', 2);
xlabel('Num iteration');
ylabel('Cost');
title('Convergence');
grid on;



%%

function [f, history] = gradient_descent(max_iter, D, g, betha, toleration)

    history = [];
    f = zeros(size(g)); 

    for k = 1:max_iter
        
        f = f + betha * D' * (g - D * f); 

        cost = (g - D * f)' * (g - D * f);

        history = [history; cost];
        if k > 1 && abs(history(k) - history(k-1)) < toleration
            break;
        end
    end
end



