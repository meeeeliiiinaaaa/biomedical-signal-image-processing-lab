%% Part 4

fileNames = {'pd.jpg', 't1.jpg', 't2.jpg'};

refImage = imread(fileNames{1});
[imgRows, imgCols, ~] = size(refImage);
totalPixels = imgRows * imgCols;

featureDataset = zeros(totalPixels, length(fileNames), 'double');

for i = 1:length(fileNames)
    tempImg = imread(fileNames{i});
    firstChannel = double(tempImg(:, :, 1));
    featureDataset(:, i) = firstChannel(:);
end

numClasses = 6;
fcmOptions = [2.0, 100, 1e-5, 0]; 

[~, membershipMat] = fcm(featureDataset, numClasses, fcmOptions);
[~, pixelLabels] = max(membershipMat, [], 1);

segmentationMap = reshape(pixelLabels, imgRows, imgCols);

fig = figure('Name', 'Tissue Segmentation via FCM', 'Color', 'white');

for classId = 1:numClasses
    binaryMask = (segmentationMap == classId);
    subplot(2, 3, classId);
    imshow(binaryMask, []);
    customTitle = sprintf('FCM Region: %d', classId);
    title(customTitle, 'FontWeight', 'bold', 'Color', [0.1 0.2 0.5]);
end