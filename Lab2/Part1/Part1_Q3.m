clearvars; clc; close all;

load('X_org.mat', 'X_org');
load('X_noise.mat', 'X_noise');

powerClean = sum(X_org(:).^2);
powerNoise = sum(X_noise(:).^2);

coeff1 = sqrt((powerClean / powerNoise) * (10^0.5));
coeff2 = sqrt((powerClean / powerNoise) * (10^1.5));


mixedSignal_1 = X_org + (coeff1 * X_noise);
mixedSignal_2 = X_org + (coeff2 * X_noise);

visualizeEEG(mixedSignal_2);

function visualizeEEG(signalMat) 
    load('Electrodes.mat', 'Electrodes');
    channelSpacing = max(abs(signalMat(:)));
    sampRate = 250;
    chNames = Electrodes.labels;
    renderMultichannelSignal(signalMat, channelSpacing, sampRate, chNames, 'Noisy EEG Data - Config 2');
end

function timeVector = renderMultichannelSignal(dataBlock, spacing, sampRate, labels, figTitle)
    [numChans, numSamples] = size(dataBlock);

    if ~exist('spacing', 'var') || isempty(spacing)
        spacing = max(abs(dataBlock(:)));
    end
    
    if ~exist('sampRate', 'var') || isempty(sampRate)
        sampRate = 1;
    end
    
    if ~exist('labels', 'var') || isempty(labels)
        labels = arrayfun(@(x) sprintf('S%d', x), 1:numChans, 'UniformOutput', false);
    end
    
    if ~exist('figTitle', 'var')
        figTitle = '';
    end

    timeVector = (1:numSamples) / sampRate;

    shiftValues = spacing * (0:-1:-(numChans-1))';
    dataBlock = dataBlock + (shiftValues * ones(1, numSamples));

    tickPositions = shiftValues;
    reverseIndices = numChans:-1:1;
    yMax = max(dataBlock(1,:)) + spacing;
    yMin = min(dataBlock(end,:)) - spacing;
    
    fig = figure('Color', 'w', 'Name', 'EEG Viewer');
    ax = axes('Parent', fig);
    
    plot(ax, timeVector, dataBlock', 'LineWidth', 0.85);
    
    axis(ax, 'tight');
    ylim(ax, [yMin yMax]);
    
    set(ax, 'YTick', tickPositions(reverseIndices), ...
            'YTickLabel', labels(reverseIndices), ...
            'FontSize', 8, 'FontWeight', 'bold', 'TickDir', 'out');
            
    grid(ax, 'on');
    set(ax, 'GridLineStyle', '--', 'GridAlpha', 0.5);
    
    xlabel(ax, 'Time [s]', 'FontSize', 11, 'FontWeight', 'bold');
    ylabel(ax, 'Electrode Channels', 'FontSize', 11, 'FontWeight', 'bold');
    
    if ~isempty(figTitle)
        title(ax, figTitle, 'FontSize', 12, 'FontWeight', 'bold', 'Color', [0.15 0.15 0.55]);
    end
end
