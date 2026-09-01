%% Part 4.1

orig_image = imread("ct.jpg");

dx = 20;
dy = 40;

% Consider first slice
orig_image = orig_image(:,:,1); 

orig_image = im2double(orig_image);
[M, N] = size(orig_image);

% Fourier transform of image
F_orig_image = fft2(orig_image);


[u, v] = meshgrid(0:N-1, 0:M-1);

% Shift kernel
H = exp(-2j * pi * (u * dx / N + v * dy / M));

F_shifted = F_orig_image .* H;

% Inverse FFT to get shifted spatial image
I_shifted = real(ifft2(F_shifted));

% Display images
subplot(1,2,1);
imshow(orig_image, []);
title('Original Image');

subplot(1,2,2);
imshow(I_shifted, []);
title(sprintf('Shifted Image (right %d, down %d)', dx, dy));

% Plot kernel
figure;
surf(1:N, 1:M, abs(H), 'EdgeColor', 'none');
xlabel('u');
ylabel('v');
zlabel('|H(u,v)|');
title('Magnitude of the Shifting Kernel H in the Frequency Domain');
view(2);
colorbar;

%%

% Zero padded to avoid circular shift

orig_image = imread("ct.jpg");

orig_image = orig_image(:,:,1);  
orig_image = im2double(orig_image);
[M, N] = size(orig_image);

dx = 20;   % shift right (columns)
dy = 40;   % shift down (rows)

% Zero-pad the image
M_pad = M + dy;
N_pad = N + dx;
padded_image = zeros(M_pad, N_pad);
padded_image(1:M, 1:N) = orig_image;


F_pad = fft2(padded_image);

[u, v] = meshgrid(0:N_pad-1, 0:M_pad-1);  % u: columns, v: rows


H = exp(-2j * pi * (u * dx / N_pad + v * dy / M_pad));


F_shifted_pad = F_pad .* H;

shifted_padded = real(ifft2(F_shifted_pad));

I_shifted_linear = shifted_padded(1:M, 1:N);

% Display
figure;
subplot(1,2,1);
imshow(orig_image, []);
title('Original Image');

subplot(1,2,2);
imshow(I_shifted_linear, []);
title(sprintf('Linearly Shifted Image (right %d, down %d)', dx, dy));

