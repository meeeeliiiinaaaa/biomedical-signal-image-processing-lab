%% Part 6

orig_image = imread("t1.jpg");

orig_image = orig_image(:,:,1); 

orig_image = im2double(orig_image);

% Sobel
sobel_x = [-1 0 1; -2 0 2; -1 0 1];
sobel_y = sobel_x';
I_sobel_x = imfilter(orig_image, sobel_x, 'replicate');
I_sobel_y = imfilter(orig_image, sobel_y, 'replicate');
sobel_mag = sqrt(I_sobel_x.^2 + I_sobel_y.^2);

% Canny
canny_edges = edge(orig_image, 'canny');

% Display
figure;

subplot(1,3,1);
imshow(orig_image, []);
title('Original Image');

subplot(1,3,2);
imshow(sobel_mag, []);
title('Sobel Magnitude');

subplot(1,3,3);
imshow(canny_edges, []);
title('Canny Edges');

colormap(gray);