%% Part 5

orig_image = imread("t1.jpg");

orig_image = orig_image(:,:,1); 

orig_image = im2double(orig_image);


% Horizontal derivative
I_left  = circshift(orig_image, [0,  1]);   % x+1
I_right = circshift(orig_image, [0, -1]);   % x-1
I_x = (I_left - I_right) / 2;

% vertical derivative 
I_up    = circshift(orig_image, [ 1, 0]);   % y+1
I_down  = circshift(orig_image, [-1, 0]);   % y-1
I_y = (I_up - I_down) / 2;

% Gradient magnitude
grad_mag = sqrt(I_x.^2 + I_y.^2);

% Display
figure;
subplot(2,2,1);
imshow(orig_image, []);
title('Original Image');

subplot(2,2,2);
imshow(I_x, []);
title('Horizontal Derivative');

subplot(2,2,3);
imshow(I_y, []);
title('Vertical Derivative');

subplot(2,2,4);
imshow(grad_mag, []);
title('Gradient Magnitude');

colormap(gray);