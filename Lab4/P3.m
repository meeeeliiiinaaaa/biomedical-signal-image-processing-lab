%% Phase 1: 
clearvars; close all; clc;
load('FiveClass_EEG.mat'); 

sRate = 256; 
filtOrder = 4;
numChans = 30;
signalLength = size(X, 1);

sig_delta = zeros(signalLength, numChans);
sig_theta = zeros(signalLength, numChans);
sig_alpha = zeros(signalLength, numChans);
sig_beta  = zeros(signalLength, numChans);

hd_delta = design(fdesign.bandpass('N,Fp1,Fp2,Ap', filtOrder, 1,  4,  1, sRate), 'cheby1');
hd_theta = design(fdesign.bandpass('N,Fp1,Fp2,Ap', filtOrder, 4,  8,  1, sRate), 'cheby1');
hd_alpha = design(fdesign.bandpass('N,Fp1,Fp2,Ap', filtOrder, 8,  13, 1, sRate), 'cheby1');
hd_beta  = design(fdesign.bandpass('N,Fp1,Fp2,Ap', filtOrder, 13, 30, 1, sRate), 'cheby1');

for ch = 1:numChans
    sig_delta(:, ch) = filter(hd_delta, X(:, ch));
    sig_theta(:, ch) = filter(hd_theta, X(:, ch));
    sig_alpha(:, ch) = filter(hd_alpha, X(:, ch));
    sig_beta(:, ch)  = filter(hd_beta, X(:, ch));
end


t_axis = 0 : 1/sRate : 5;
idx_5s = 1 : (5 * sRate + 1);

figure('Name', 'Channel 1 Sub-bands', 'Color', 'w');
plot_data = {X(idx_5s, 1), sig_delta(idx_5s, 1), sig_theta(idx_5s, 1), sig_alpha(idx_5s, 1), sig_beta(idx_5s, 1)};
titles_p1 = {'Raw EEG Signal', 'Delta Rhythm (1-4 Hz)', 'Theta Rhythm (4-8 Hz)', 'Alpha Rhythm (8-13 Hz)', 'Beta Rhythm (13-30 Hz)'};
colors_p1 = {'#000000', '#D95319', '#EDB120', '#7E2F8E', '#77AC30'}; % Modern colors

for p = 1:5
    subplot(5, 1, p);
    plot(t_axis, plot_data{p}, 'Color', colors_p1{p}, 'LineWidth', 1.2);
    title(titles_p1{p}, 'FontWeight', 'bold');
    xlabel('Time [s]'); ylabel('Voltage [\muV]');
    grid on; axis tight; box on;
end

%% Phase 2: 
numTrials = 200;
samplesPerTrial = 10 * sRate;

epoch_alpha = zeros(samplesPerTrial, numChans, numTrials);
epoch_beta  = zeros(samplesPerTrial, numChans, numTrials);
epoch_delta = zeros(samplesPerTrial, numChans, numTrials);
epoch_theta = zeros(samplesPerTrial, numChans, numTrials);

for tr = 1:numTrials
    idx_range = trial(tr) : (trial(tr) + samplesPerTrial - 1);
    epoch_alpha(:, :, tr) = sig_alpha(idx_range, :);
    epoch_beta(:, :, tr)  = sig_beta(idx_range, :);
    epoch_delta(:, :, tr) = sig_delta(idx_range, :);
    epoch_theta(:, :, tr) = sig_theta(idx_range, :);
end

%% Phase 3 & 4: 
numClasses = 5;

meanPwr_alpha = zeros(samplesPerTrial, numChans, numClasses);
meanPwr_beta  = zeros(samplesPerTrial, numChans, numClasses);
meanPwr_delta = zeros(samplesPerTrial, numChans, numClasses);
meanPwr_theta = zeros(samplesPerTrial, numChans, numClasses);

for cls = 1:numClasses
    
    class_mask = (y(:, 1) == cls); 
    
    meanPwr_alpha(:, :, cls) = mean(epoch_alpha(:, :, class_mask).^2, 3);
    meanPwr_beta(:, :, cls)  = mean(epoch_beta(:, :, class_mask).^2, 3);
    meanPwr_delta(:, :, cls) = mean(epoch_delta(:, :, class_mask).^2, 3);
    meanPwr_theta(:, :, cls) = mean(epoch_theta(:, :, class_mask).^2, 3);
end

%% Phase 5: 
winLen = 200;
smoothKernel = ones(winLen, 1) / sqrt(winLen); 

smooth_alpha = zeros(size(meanPwr_alpha));
smooth_beta  = zeros(size(meanPwr_beta));
smooth_delta = zeros(size(meanPwr_delta));
smooth_theta = zeros(size(meanPwr_theta));

for cls = 1:numClasses
    for ch = 1:numChans
        smooth_alpha(:, ch, cls) = conv(meanPwr_alpha(:, ch, cls), smoothKernel, 'same');
        smooth_beta(:, ch, cls)  = conv(meanPwr_beta(:, ch, cls), smoothKernel, 'same');
        smooth_delta(:, ch, cls) = conv(meanPwr_delta(:, ch, cls), smoothKernel, 'same');
        smooth_theta(:, ch, cls) = conv(meanPwr_theta(:, ch, cls), smoothKernel, 'same');
    end
end

%% Phase 6: 
targetChannel = 16; 
plotDur = 10;
numPoints = plotDur * sRate;
t_vec = linspace(0, plotDur, numPoints);

bandsToDisp = {'Alpha', 'Beta', 'Delta', 'Theta'};
dataToDisp = {smooth_alpha, smooth_beta, smooth_delta, smooth_theta};
classColors = lines(numClasses); 

for bnd = 1:4
    figure('Name', sprintf('Band Power: %s', bandsToDisp{bnd}), 'Color', 'w');
    sgtitle(sprintf('Mean %s Power Envelope | Channel: CPz', bandsToDisp{bnd}), 'FontSize', 12, 'FontWeight', 'bold');
    
    for cls = 1:numClasses
        subplot(5, 1, cls);
        plot(t_vec, dataToDisp{bnd}(1:numPoints, targetChannel, cls), ...
            'Color', classColors(cls,:), 'LineWidth', 1.5);
        
        title(sprintf('Task Condition: Class %d', cls), 'FontSize', 9);
        xlabel('Time [s]');
        ylabel('Power [\muV^2]');
        grid on; box on;
    end
end