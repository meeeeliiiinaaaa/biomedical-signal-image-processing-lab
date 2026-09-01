clearvars; clc; close all;

load('X_noise.mat');
visualizeEEGData(X_noise);

function visualizeEEGData(signalMat)
    load('Electrodes.mat'); 

    fs = 250; 
    channelSpacing = max(abs(signalMat(:)));
    chNames = Electrodes.labels;
    renderMultichannelSignal(signalMat, channelSpacing, fs, chNames, 'Noisy EEG Signals');
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
        labels = arrayfun(@(x) sprintf('Ch%d', x), 1:numChans, 'UniformOutput', false);
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
   
    fig = figure('Color', 'w', 'Name', 'EEG Data Viewer');
    ax = axes('Parent', fig);
    
    plot(ax, timeVector, dataBlock', 'LineWidth', 0.8);
    
    
    axis(ax, 'tight');
    ylim(ax, [yMin yMax]);
    
    set(ax, 'YTick', tickPositions(reverseIndices), ...
            'YTickLabel', labels(reverseIndices), ...
            'FontSize', 8, 'FontWeight', 'bold', 'TickDir', 'out');
            
    grid(ax, 'on');
    set(ax, 'GridLineStyle', '--', 'GridAlpha', 0.6);
    
    xlabel(ax, 'Time [s]', 'FontSize', 11, 'FontWeight', 'bold');
    ylabel(ax, 'Electrode Channels', 'FontSize', 11, 'FontWeight', 'bold');
    
    if ~isempty(figTitle)
        title(ax, figTitle, 'FontSize', 13, 'FontWeight', 'bold', 'Color', [0.1 0.1 0.4]);
    end
end
