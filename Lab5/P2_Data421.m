%% Part2

%%

load("n_421.mat");


%%

Fs = 250;
time = (0:length(n_421(:,1)) - 1) / Fs;

%% Section3

window_length = 10 * Fs;
overlap = 5 * Fs;
step = window_length - overlap;
N = length(n_421(:, 1));

%%

% Read the file
filename = "atr_n421.txt";  
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
    {'(N','(VFIB', '(NOISE'}, ...
    [1,        2,     3] );


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

save('labels_n421.mat', 'labels', 'end_samples', 'windows');


%% Section6

[alarm_medfreq, t_medfreq] = va_detect_medfreq_data424(n_421, Fs);


%% Section7

load('labels_n421.mat', 'labels'); 
evaluate(alarm_medfreq, labels, 'Median Frequency Data424');

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
yticks([1 2 3 4]);
ylabel('Label');
title('True vs Predicted Labels For Median Frequency');
legend('Location', 'best');
grid on;


%%

function [alarm,t] = va_detect_medfreq_data424(ecg_data,Fs)
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
