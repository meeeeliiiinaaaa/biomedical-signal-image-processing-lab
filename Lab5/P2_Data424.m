%% Part2

%%

load("n_424.mat");


%%

Fs = 250;
time = (0:length(n_424(:,1)) - 1) / Fs;

%% Section3

window_length = 10 * Fs;
overlap = 5 * Fs;
step = window_length - overlap;
N = length(n_424(:, 1));

%%

% Read the file
filename = "atr_n424.txt";   
lines = readlines(filename, 'EmptyLineRule', 'skip');

% Preallocate cell arrays for results
timeStrs = {};
sampleVals = [];
eventStrs = {};

timePattern = '^\d{2}:\d{2}[.:]\d{3}';

for i = 1:length(lines)
    line = strtrim(lines(i));
    if isempty(line)
        continue;
    end
 
    if isempty(regexp(line, timePattern, 'once'))
        continue;
    end
    
    % Split line by whitespace (tabs/spaces)
    tokens = strsplit(line);
    % time, sample, '+', '0','0','0', event
    if length(tokens) < 7
        warning('Line %d has fewer tokens than expected: %s', i, line);
        continue;
    end
    
    timeStr = tokens{1};
    sampleStr = tokens{2};
    eventStr = strtrim(tokens{7});   % event is the 7th token
    
    sampleNum = str2double(sampleStr);
    if isnan(sampleNum)
        warning('Invalid sample number on line %d: %s', i, sampleStr);
        continue;
    end
    
    % Store
    timeStrs{end+1} = timeStr;
    sampleVals(end+1) = sampleNum;
    eventStrs{end+1} = eventStr;
end

% Display results
disp('Extracted data:');
for k = 1:length(timeStrs)
    fprintf('Time: %s, Sample: %d, Event: %s\n', timeStrs{k}, sampleVals(k), eventStrs{k});
end


%%

map_label = containers.Map( ...
    {'(N','(VFIB','(ASYS','(NOISE', '(NOD'}, ...
    [1,        2,     3,    4,  5] );


%%

windows     = {};
labels      = [];
end_samples = [];


sentinel = N + 1;
sampleVals_ext = [sampleVals, sentinel];

idx = 1;
w   = 1;

while idx + window_length - 1 <= N

    s1 = idx;
    s2 = idx + window_length - 1;
    end_samples(w) = s2;


    found_labels = [];

    for k = 1:length(eventStrs)        
        interval_start = sampleVals_ext(k);
        interval_end   = sampleVals_ext(k+1) - 1;

      
        if interval_start <= s2 && interval_end >= s1
            labelStr = eventStrs{k};
            if isKey(map_label, labelStr)
                found_labels(end+1) = map_label(labelStr);
            end
        end
    end

   
    unique_labels = unique(found_labels);

    if isempty(unique_labels)
        labels(w) = 0;                   
    elseif length(unique_labels) == 1
        labels(w) = unique_labels(1);    
    else
        labels(w) = 0;                   
    end

    windows{w} = n_422(s1:s2, 1);       
    idx = idx + step;
    w   = w + 1;
end

fprintf('\nTotal windows created = %d\n', length(windows));
fprintf('\n%-10s  %-12s  %-12s  %s\n', 'Window', 'Start Samp', 'End Samp', 'Label');
fprintf('%s\n', repmat('-', 1, 48));
for i = 1:length(windows)
    fprintf('%-10d  %-12d  %-12d  %d\n', ...
            i, end_samples(i) - window_length + 1, end_samples(i), labels(i));
end

save('labels_n424.mat', 'labels', 'end_samples', 'windows');

%% Section4

load('labels_n424.mat', 'labels', 'end_samples', 'windows');

numWindows = length(windows);

medFreqs   = zeros(numWindows, 1);
meanFreqs  = zeros(numWindows, 1);
bandPowers = zeros(numWindows, 1);

for w = 1:numWindows
    seg = windows{w};              
    
 
    if isrow(seg)
        seg = seg';
    end
    
    [seg_psd, f] = pwelch(seg, [], [], [], Fs); 
    % Median frequency (Hz)
    medFreqs(w) = medfreq(seg_psd, f);
    
    % Mean frequency (Hz)
    meanFreqs(w) = meanfreq(seg_psd, f);
    
    % Total band power (0 to Nyquist frequency = Fs/2)
    bandPowers(w) = bandpower(seg_psd, f, [0 Fs/2], 'psd');
end

T = table((1:numWindows)', end_samples', labels', ...
          medFreqs, meanFreqs, bandPowers, ...
          'VariableNames', {'Window', 'EndSample', 'Label', ...
                            'MedFreq_Hz', 'MeanFreq_Hz', 'BandPower'});
disp(T);

save('window_freq_features.mat', 'medFreqs', 'meanFreqs', 'bandPowers', 'T');
fprintf('Features saved to window_freq_features.mat\n');

%% Section5

load('window_freq_features.mat', 'medFreqs', 'meanFreqs', 'bandPowers', 'T');

labels = T.Label;

idx_class1 = (labels == 1);
idx_class2 = (labels == 2);

med_class1 = medFreqs(idx_class1);
med_class2 = medFreqs(idx_class2);

mean_class1 = meanFreqs(idx_class1);
mean_class2 = meanFreqs(idx_class2);

power_class1 = bandPowers(idx_class1);
power_class2 = bandPowers(idx_class2);

% plot
figure('Position', [100 100 1200 400]);

% Median Frequency
subplot(1,3,1);
histogram(med_class1, 'FaceColor', 'b', 'EdgeColor', 'k', 'Normalization', 'probability', 'BinWidth', 1);
hold on;
histogram(med_class2, 'FaceColor', 'r', 'EdgeColor', 'k', 'Normalization', 'probability', 'BinWidth', 1);
xlabel('Median Frequency (Hz)');
ylabel('Probability');
title('Median Frequency Distribution');
legend('Class 1 (N)', 'Class 2 (VFIB)');
grid on;

% Mean Frequency
subplot(1,3,2);
histogram(mean_class1, 'FaceColor', 'b', 'EdgeColor', 'k', 'Normalization', 'probability', 'BinWidth', 1);
hold on;
histogram(mean_class2, 'FaceColor', 'r', 'EdgeColor', 'k', 'Normalization', 'probability', 'BinWidth', 1);
xlabel('Mean Frequency (Hz)');
ylabel('Probability');
title('Mean Frequency Distribution');
legend('Class 1 (N)', 'Class 2 (VFIB)');
grid on;

% Band Power 
subplot(1,3,3);
histogram(log10(power_class1), 'FaceColor', 'b', 'EdgeColor', 'k', 'Normalization', 'probability', 'BinWidth', 0.2);
hold on;
histogram(log10(power_class2), 'FaceColor', 'r', 'EdgeColor', 'k', 'Normalization', 'probability', 'BinWidth', 0.2);
xlabel('Log10(Band Power)');
ylabel('Probability');
title('Total Band Power Distribution');
legend('Class 1 (N)', 'Class 2 (VFIB)');
grid on;

sgtitle('Comparison of Frequency Features for Class 1 vs Class 2');

%% Section6

[alarm_medfreq, t_medfreq] = va_detect_medfreq(n_424, Fs);
[alarm_meanfreq, t_meanfreq] = va_detect_meanfreq(n_424,Fs);

%% Section7

load('labels_n424.mat', 'labels'); 
evaluate(alarm_medfreq, labels, 'Median Frequency');
evaluate(alarm_meanfreq, labels, 'Mean Frequency');

%%

% Scatter plot
y_true = labels;
y_pred= alarm_medfreq + 1;


scatter(1:length(y_true), y_true,  40, 'b', 'o', 'filled', ...
    'DisplayName', 'y\_true', 'MarkerFaceAlpha', 0.5);
hold on;
scatter(1:length(y_pred), y_pred,  40, 'r', 'x', 'LineWidth', 1.5, ...
    'DisplayName', 'y\_pred');

xlabel('Window index');
yticks([1 2 3 4 5]);
ylabel('Label');
title('True vs Predicted Labels For Median Frequency');
legend('Location', 'best');
grid on;

%%

% Scatter plot 
y_true = labels;
y_pred= alarm_meanfreq + 1;


scatter(1:length(y_true), y_true,  40, 'b', 'o', 'filled', ...
    'DisplayName', 'y\_true', 'MarkerFaceAlpha', 0.5);
hold on;
scatter(1:length(y_pred), y_pred,  40, 'r', 'x', 'LineWidth', 1.5, ...
    'DisplayName', 'y\_pred');

xlabel('Window index');
yticks([1 2 3 4 5]);
ylabel('Label');
title('True vs Predicted Labels For Mean Frequency');
legend('Location', 'best');
grid on;

%% Section8

load('labels_n424.mat', 'labels', 'end_samples', 'windows');

numWindows = length(windows);

max_amplitude   = zeros(numWindows, 1);
min_amplitude  = zeros(numWindows, 1);
peaktopeaks = zeros(numWindows, 1);
mean_sig = zeros(numWindows, 1);
var_sig = zeros(numWindows, 1);

for w = 1:numWindows
    seg = windows{w};              
    
 
    if isrow(seg)
        seg = seg';
    end
    
    max_amplitude(w)      = max(seg);
    min_amplitude(w)      = min(seg);
    peaktopeaks(w) = max_amplitude(w)  - min_amplitude(w);
    mean_sig(w)     = mean(seg);
    var_sig(w)     = var(seg);
end

T = table((1:numWindows)', end_samples', labels', ...
          max_amplitude, min_amplitude, peaktopeaks, mean_sig, var_sig, ...
          'VariableNames', {'Window', 'EndSample', 'Label', ...
                            'max_amplitude', 'min_amplitude', 'peaktopeaks', 'mean_sig', 'var_sig'});

save('window_time_features.mat', 'max_amplitude', 'min_amplitude', 'peaktopeaks', 'mean_sig', 'var_sig', 'T');
fprintf('Features saved to window_time_features.mat\n');

%% Section9

load('window_time_features.mat', 'max_amplitude', 'min_amplitude', 'peaktopeaks', 'mean_sig', 'var_sig', 'T');

labels = T.Label;

idx_class1 = (labels == 1);
idx_class2 = (labels == 2);

max_class1 = max_amplitude(idx_class1);
max_class2 = max_amplitude(idx_class2);

min_class1 = min_amplitude(idx_class1);
min_class2 = min_amplitude(idx_class2);

peaktopeaks_class1 = peaktopeaks(idx_class1);
peaktopeaks_class2 = peaktopeaks(idx_class2);

mean_class1 = mean_sig(idx_class1);
mean_class2 = mean_sig(idx_class2);

var_class1 = var_sig(idx_class1);
var_class2 = var_sig(idx_class2); 

% plot
figure('Position', [100 100 1200 400]);

% Max
subplot(2,3,1);
histogram(max_class1, 'FaceColor', 'b', 'EdgeColor', 'k', 'Normalization', 'probability', 'BinWidth', 1);
hold on;
histogram(max_class2, 'FaceColor', 'r', 'EdgeColor', 'k', 'Normalization', 'probability', 'BinWidth', 1);
xlabel('Max');
ylabel('Probability');
title('Max Amplitude');
legend('Class 1 (N)', 'Class 2 (VFIB)');
grid on;

% Min
subplot(2,3,2);
histogram(min_class1, 'FaceColor', 'b', 'EdgeColor', 'k', 'Normalization', 'probability', 'BinWidth', 1);
hold on;
histogram(min_class2, 'FaceColor', 'r', 'EdgeColor', 'k', 'Normalization', 'probability', 'BinWidth', 1);
xlabel('Min');
ylabel('Probability');
title('Min Amplitude');
legend('Class 1 (N)', 'Class 2 (VFIB)');
grid on;

% Peaktopeaks
subplot(2,3,3);
histogram(peaktopeaks_class1, 'FaceColor', 'b', 'EdgeColor', 'k', 'Normalization', 'probability', 'BinWidth', 1);
hold on;
histogram(peaktopeaks_class2, 'FaceColor', 'r', 'EdgeColor', 'k', 'Normalization', 'probability', 'BinWidth', 1);
xlabel('Peak to Peak');
ylabel('Probability');
title('Amplitude Peak to Peak');
legend('Class 1 (N)', 'Class 2 (VFIB)');
grid on;

% Mean
subplot(2,3,4);
histogram(mean_class1, 'FaceColor', 'b', 'EdgeColor', 'k', 'Normalization', 'probability', 'BinWidth', 1);
hold on;
histogram(mean_class2, 'FaceColor', 'r', 'EdgeColor', 'k', 'Normalization', 'probability', 'BinWidth', 1);
xlabel('Mean');
ylabel('Probability');
title('Amplitude Mean');
legend('Class 1 (N)', 'Class 2 (VFIB)');
grid on;

% Var
subplot(2,3,5);
histogram(var_class1, 'FaceColor', 'b', 'EdgeColor', 'k', 'Normalization', 'probability', 'BinWidth', 1);
hold on;
histogram(var_class2, 'FaceColor', 'r', 'EdgeColor', 'k', 'Normalization', 'probability', 'BinWidth', 1);
xlabel('Variance');
ylabel('Probability');
title('Amplitude Variance');
legend('Class 1 (N)', 'Class 2 (VFIB)');
grid on;

subplot(2, 3, 6);
axis off;

sgtitle('Comparison of Time Features for Class 1 vs Class 2');

%% Section10

[alarm_mean, t_mean] = va_detect_mean(n_424, Fs);
[alarm_max, t_min] = va_detect_max(n_424,Fs);

%% Section11

load('labels_n424.mat', 'labels'); 
evaluate(alarm_mean, labels, 'Mean Amplitude');
evaluate(alarm_max, labels, 'Max Amplitude');

%%

% Scatter plot
y_true = labels;
y_pred= alarm_mean + 1;


scatter(1:length(y_true), y_true,  40, 'b', 'o', 'filled', ...
    'DisplayName', 'y\_true', 'MarkerFaceAlpha', 0.5);
hold on;
scatter(1:length(y_pred), y_pred,  40, 'r', 'x', 'LineWidth', 1.5, ...
    'DisplayName', 'y\_pred');

xlabel('Window index');
yticks([1 2 3 4]);
ylabel('Label');
title('True vs Predicted Labels For Mean Amplitude');
legend('Location', 'best');
grid on;

%%

% Scatter plot 
y_true = labels;
y_pred= alarm_max + 1;


scatter(1:length(y_true), y_true,  40, 'b', 'o', 'filled', ...
    'DisplayName', 'y\_true', 'MarkerFaceAlpha', 0.5);
hold on;
scatter(1:length(y_pred), y_pred,  40, 'r', 'x', 'LineWidth', 1.5, ...
    'DisplayName', 'y\_pred');

xlabel('Window index');
yticks([1 2 3 4]);
ylabel('Label');
title('True vs Predicted Labels For Max Amplitude');
legend('Location', 'best');
grid on;


%% 

function [alarm,t] = va_detect_medfreq(ecg_data,Fs)
%VA_DETECT  ventricular arrhythmia detection skeleton function
%  [ALARM,T] = VA_DETECT(ECG_DATA,FS) is a skeleton function for ventricular
%  arrhythmia detection, designed to help you get started in implementing your
%  arrhythmia detector.
%
%  This code automatically sets up fixed length data frames, stepping through 
%  the entire ECG waveform with 50% overlap of consecutive frames. You can customize 
%  the frame length  by adjusting the internal 'frame_sec' variable and the overlap by
%  adjusting the 'overlap' variable.
%
%  ECG_DATA is a vector containing the ecg signal, and FS is the sampling rate
%  of ECG_DATA in Hz. The output ALARM is a vector of ones and zeros
%  corresponding to the time frames for which the alarm is active (1) 
%  and inactive (0). T is a vector the same length as ALARM which contains the 
%  time markers which correspond to the end of each analyzed time segment. If Fs 
%  is not entered, the default value of 250 Hz is used. 

  %  Template Last Modified: 3/4/06 by Eric Weiss, 1/25/07 by Julie Greenberg


%  Processing frames: adjust frame length & overlap here
%------------------------------------------------------
frame_sec = 10;  % sec
overlap = 0.5;    % 50% overlap between consecutive frames


% Input argument checking
%------------------------
if nargin < 2
    Fs = 250;  % default sample rate
end;
if nargin < 1
    error('You must enter an ECG data vector.');
end;
ecg_data = ecg_data(:);  % Make sure that ecg_data is a column vector
 

% Initialize Variables
%---------------------
frame_length = round(frame_sec*Fs);  % length of each data frame (samples)
frame_step = round(frame_length*(1-overlap));  % amount to advance for next data frame
ecg_length = length(ecg_data);  % length of input vector
frame_N = floor((ecg_length-(frame_length-frame_step))/frame_step); % total number of frames
alarm = zeros(frame_N,1);	% initialize output signal to all zeros
t = ([0:frame_N-1]*frame_step+frame_length)/Fs;

% Analysis loop: each iteration processes one frame of data
%----------------------------------------------------------
for i = 1:frame_N
    %  Get the next data segment
    seg = ecg_data(((i-1)*frame_step+1):((i-1)*frame_step+frame_length));
    %  Perform computations on the segment . . .
    [seg_psd, f] = pwelch(seg, [], [], [], Fs);
    seg_medfreq = medfreq(seg_psd, f);

    if seg_medfreq  >= 3 && seg_medfreq  <= 5
     
        alarm(i) = 1;

    end
  
end


end

%%


function [alarm,t] = va_detect_meanfreq(ecg_data,Fs)
%VA_DETECT  ventricular arrhythmia detection skeleton function
%  [ALARM,T] = VA_DETECT(ECG_DATA,FS) is a skeleton function for ventricular
%  arrhythmia detection, designed to help you get started in implementing your
%  arrhythmia detector.
%
%  This code automatically sets up fixed length data frames, stepping through 
%  the entire ECG waveform with 50% overlap of consecutive frames. You can customize 
%  the frame length  by adjusting the internal 'frame_sec' variable and the overlap by
%  adjusting the 'overlap' variable.
%
%  ECG_DATA is a vector containing the ecg signal, and FS is the sampling rate
%  of ECG_DATA in Hz. The output ALARM is a vector of ones and zeros
%  corresponding to the time frames for which the alarm is active (1) 
%  and inactive (0). T is a vector the same length as ALARM which contains the 
%  time markers which correspond to the end of each analyzed time segment. If Fs 
%  is not entered, the default value of 250 Hz is used. 

  %  Template Last Modified: 3/4/06 by Eric Weiss, 1/25/07 by Julie Greenberg


%  Processing frames: adjust frame length & overlap here
%------------------------------------------------------
frame_sec = 10;  % sec
overlap = 0.5;    % 50% overlap between consecutive frames


% Input argument checking
%------------------------
if nargin < 2
    Fs = 250;  % default sample rate
end;
if nargin < 1
    error('You must enter an ECG data vector.');
end;
ecg_data = ecg_data(:);  % Make sure that ecg_data is a column vector
 

% Initialize Variables
%---------------------
frame_length = round(frame_sec*Fs);  % length of each data frame (samples)
frame_step = round(frame_length*(1-overlap));  % amount to advance for next data frame
ecg_length = length(ecg_data);  % length of input vector
frame_N = floor((ecg_length-(frame_length-frame_step))/frame_step); % total number of frames
alarm = zeros(frame_N,1);	% initialize output signal to all zeros
t = ([0:frame_N-1]*frame_step+frame_length)/Fs;

% Analysis loop: each iteration processes one frame of data
%----------------------------------------------------------
for i = 1:frame_N
    %  Get the next data segment
    seg = ecg_data(((i-1)*frame_step+1):((i-1)*frame_step+frame_length));
    %  Perform computations on the segment . . .
  
    [seg_psd, f] = pwelch(seg, [], [], [], Fs);
    seg_meanfreq = meanfreq(seg_psd, f);

    if  seg_meanfreq >= 3 && seg_meanfreq <= 4

        alarm(i) = 1;
       
    end

    
end


end

%%

function results = evaluate(alarm, labels, detector_name)
        
    mask = (labels == 1) | (labels == 2);
    y_true = labels(mask);
    y_pred = alarm(mask)';
    
    y_true_bin = (y_true == 2);   % 1 for VFIB, 0 for N
  
    % Confusion matrix
    TN = sum((y_true_bin == 0) & (y_pred == 0));
    FP = sum((y_true_bin == 0) & (y_pred == 1));
    FN = sum((y_true_bin == 1) & (y_pred == 0));
    TP = sum((y_true_bin == 1) & (y_pred == 1));
    
    confMat = [TN, FP; FN, TP];
    
    % Metrics
    accuracy = (TP + TN) / (TP + TN + FP + FN);
    sensitivity = TP / (TP + FN);   % recall for VFIB
    specificity = TN / (TN + FP);   % recall for N
    
    % Display
    fprintf(detector_name);
    fprintf('\n');
    fprintf('Accuracy     = %.3f\n', accuracy);
    fprintf('Sensitivity  = %.3f\n', sensitivity);
    fprintf('Specificity  = %.3f\n', specificity);
    
    % Plot confusion matrix
    figure;
    confusionchart(confMat, {'N', 'VFIB'}, ...
        'Title', sprintf('Confusion Matrix: %s Detector', detector_name));
    
    % Store results
    results = struct('accuracy', accuracy, 'sensitivity', sensitivity, ...
                 'specificity', specificity, 'confusionMatrix', confMat);
end


%% 

function [alarm,t] = va_detect_mean(ecg_data,Fs)
%VA_DETECT  ventricular arrhythmia detection skeleton function
%  [ALARM,T] = VA_DETECT(ECG_DATA,FS) is a skeleton function for ventricular
%  arrhythmia detection, designed to help you get started in implementing your
%  arrhythmia detector.
%
%  This code automatically sets up fixed length data frames, stepping through 
%  the entire ECG waveform with 50% overlap of consecutive frames. You can customize 
%  the frame length  by adjusting the internal 'frame_sec' variable and the overlap by
%  adjusting the 'overlap' variable.
%
%  ECG_DATA is a vector containing the ecg signal, and FS is the sampling rate
%  of ECG_DATA in Hz. The output ALARM is a vector of ones and zeros
%  corresponding to the time frames for which the alarm is active (1) 
%  and inactive (0). T is a vector the same length as ALARM which contains the 
%  time markers which correspond to the end of each analyzed time segment. If Fs 
%  is not entered, the default value of 250 Hz is used. 

  %  Template Last Modified: 3/4/06 by Eric Weiss, 1/25/07 by Julie Greenberg


%  Processing frames: adjust frame length & overlap here
%------------------------------------------------------
frame_sec = 10;  % sec
overlap = 0.5;    % 50% overlap between consecutive frames


% Input argument checking
%------------------------
if nargin < 2
    Fs = 250;  % default sample rate
end;
if nargin < 1
    error('You must enter an ECG data vector.');
end;
ecg_data = ecg_data(:);  % Make sure that ecg_data is a column vector
 

% Initialize Variables
%---------------------
frame_length = round(frame_sec*Fs);  % length of each data frame (samples)
frame_step = round(frame_length*(1-overlap));  % amount to advance for next data frame
ecg_length = length(ecg_data);  % length of input vector
frame_N = floor((ecg_length-(frame_length-frame_step))/frame_step); % total number of frames
alarm = zeros(frame_N,1);	% initialize output signal to all zeros
t = ([0:frame_N-1]*frame_step+frame_length)/Fs;

% Analysis loop: each iteration processes one frame of data
%----------------------------------------------------------
for i = 1:frame_N
    %  Get the next data segment
    seg = ecg_data(((i-1)*frame_step+1):((i-1)*frame_step+frame_length));
    %  Perform computations on the segment . . .
    seg_mean = mean(seg);

    if seg_mean  >= -37 && seg_mean  <= -26
     
        alarm(i) = 1;

    end
  
end


end

%%

function [alarm,t] = va_detect_max(ecg_data,Fs)
%VA_DETECT  ventricular arrhythmia detection skeleton function
%  [ALARM,T] = VA_DETECT(ECG_DATA,FS) is a skeleton function for ventricular
%  arrhythmia detection, designed to help you get started in implementing your
%  arrhythmia detector.
%
%  This code automatically sets up fixed length data frames, stepping through 
%  the entire ECG waveform with 50% overlap of consecutive frames. You can customize 
%  the frame length  by adjusting the internal 'frame_sec' variable and the overlap by
%  adjusting the 'overlap' variable.
%
%  ECG_DATA is a vector containing the ecg signal, and FS is the sampling rate
%  of ECG_DATA in Hz. The output ALARM is a vector of ones and zeros
%  corresponding to the time frames for which the alarm is active (1) 
%  and inactive (0). T is a vector the same length as ALARM which contains the 
%  time markers which correspond to the end of each analyzed time segment. If Fs 
%  is not entered, the default value of 250 Hz is used. 

  %  Template Last Modified: 3/4/06 by Eric Weiss, 1/25/07 by Julie Greenberg


%  Processing frames: adjust frame length & overlap here
%------------------------------------------------------
frame_sec = 10;  % sec
overlap = 0.5;    % 50% overlap between consecutive frames


% Input argument checking
%------------------------
if nargin < 2
    Fs = 250;  % default sample rate
end;
if nargin < 1
    error('You must enter an ECG data vector.');
end;
ecg_data = ecg_data(:);  % Make sure that ecg_data is a column vector
 

% Initialize Variables
%---------------------
frame_length = round(frame_sec*Fs);  % length of each data frame (samples)
frame_step = round(frame_length*(1-overlap));  % amount to advance for next data frame
ecg_length = length(ecg_data);  % length of input vector
frame_N = floor((ecg_length-(frame_length-frame_step))/frame_step); % total number of frames
alarm = zeros(frame_N,1);	% initialize output signal to all zeros
t = ([0:frame_N-1]*frame_step+frame_length)/Fs;

% Analysis loop: each iteration processes one frame of data
%----------------------------------------------------------
for i = 1:frame_N
    %  Get the next data segment
    seg = ecg_data(((i-1)*frame_step+1):((i-1)*frame_step+frame_length));
    %  Perform computations on the segment . . .
    seg_max = max(seg);

    if seg_max  >= 210 && seg_max  <= 300
     
        alarm(i) = 1;

    end
  
end


end









