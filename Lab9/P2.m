%% Part2

pd = imread("pd.jpg");
t1 = imread("t1.jpg");
t2 = imread("t2.jpg");

pd = im2double(pd(:,:,1));
t1 = im2double(t1(:,:,1));
t2 = im2double(t2(:,:,1));

k = 6;
[m,n] = size(pd);

feature = [pd(:), t1(:), t2(:)];
labels = kmeans(feature, k);

clustered_image = reshape(labels,m,n);

% Display
figure;

subplot(3,3,2)
imshow(t1,[])
title('Original T1')

for i = 1:k

    subplot(3,3,i+3)

    cluster = (clustered_image == i);

    imshow(cluster)

    title(['Cluster ',num2str(i)])

end
