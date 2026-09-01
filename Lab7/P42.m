%% Part 4.2

orig_image = imread("ct.jpg");

angle = 30;

% Consider first slice
orig_image = orig_image(:,:,1); 

orig_image = im2double(orig_image);
[M, N] = size(orig_image);

% Kernel
rotated_image = imrotate(orig_image, -angle, 'bilinear', 'crop');

F_orig = fftshift(fft2(orig_image));
F_rotated = fftshift(fft2(rotated_image));

% Display images
figure;
subplot(1,2,1);
imshow(orig_image, []);
title('Original Image');
subplot(1,2,2);
imshow(rotated_image, []);
title(sprintf('Rotated by %d°', angle));

% Dispaly in frequency domain
mag_orig = log(1 + abs(F_orig));
mag_rot = log(1 + abs(F_rotated));

figure;
subplot(1,2,1);
imshow(mag_orig, []);
title('Fourier Magnitude (Original)');
colormap('jet'); colorbar;
subplot(1,2,2);
imshow(mag_rot, []);
title('Fourier Magnitude (Rotated)');
colormap('jet'); colorbar;

%%
% Image rotation using frequency-domain processing

theta = -deg2rad(angle);

F = fft2(fftshift(orig_image));

[U, V] = meshgrid(1:N, 1:M);
U_centered = U - ceil(N/2) - 1;
V_centered = V - ceil(M/2) - 1;

% Rotate the grid by theta 
U_rot = U_centered * cos(theta) - V_centered * sin(theta);
V_rot = U_centered * sin(theta) + V_centered * cos(theta);

U_rot = U_rot + ceil(N/2) + 1;
V_rot = V_rot + ceil(M/2) + 1;

F_rot = interp2(U, V, F, U_rot, V_rot, 'linear', 0);

image_rotated = real(ifftshift(ifft2(F_rot)));

% Display 
figure;
subplot(1,2,1);
imshow(orig_image, []);
title('Original Image');

subplot(1,2,2);
imshow(image_rotated, []);
title('Rotated via Frequency Domain');


