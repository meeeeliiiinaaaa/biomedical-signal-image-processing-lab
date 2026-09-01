%% Part 3

scan_pd = imread("pd.jpg");
scan_t1 = imread("t1.jpg");
scan_t2 = imread("t2.jpg");

[rowSize, colSize, ~] = size(scan_pd);
totalPixels = rowSize * colSize;

feature1 = double(reshape(scan_pd(:,:,1), totalPixels, 1));
feature2 = double(reshape(scan_t1(:,:,1), totalPixels, 1));
feature3 = double(reshape(scan_t2(:,:,1), totalPixels, 1));
featureMatrix = [feature1, feature2, feature3];


numSegments = 6;

segmentLabels = custom_kmeans_clustering(featureMatrix, numSegments);
segmentedImage = reshape(segmentLabels, rowSize, colSize);


fig = figure('Name', 'Image Segmentation Results', 'Color', 'white');

for idx = 1:numSegments
    currentMask = (segmentedImage == idx);
    subplot(2, 3, idx);
    imshow(currentMask, []);
    title(sprintf('Segment ID: %d', idx), 'FontWeight', 'bold', 'FontSize', 10);
end


function labels_out = custom_kmeans_clustering(data, K)
    [numPoints, numFeatures] = size(data);
    random_idx = randperm(numPoints, K);
    centroids = data(random_idx, :);
    
    labels_out = zeros(numPoints, 1);
    maxIterations = 300;
    iterCount = 0;
    
    while iterCount < maxIterations
        iterCount = iterCount + 1;
        distMatrix = pdist2(data, centroids, 'euclidean');
        [~, labels_out] = min(distMatrix, [], 2);

        updated_centroids = zeros(K, numFeatures);
        for g = 1:K
            clusterData = data(labels_out == g, :);
            if ~isempty(clusterData)
                updated_centroids(g, :) = mean(clusterData, 1);
            else
                updated_centroids(g, :) = data(randi(numPoints), :); 
            end
        end
        
        diff_error = sum(abs(centroids(:) - updated_centroids(:)));
        if diff_error < 1e-5
            break; 
        end
        centroids = updated_centroids;
    end
end