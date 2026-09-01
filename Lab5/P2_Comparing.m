%% Part2

%%

Fs = 250;

%%

[alarm_medfreq_422, t_medfreq_422] = va_detect_medfreq_data424(n_422, Fs);
[alarm_medfreq_424, t_medfreq_424] = va_detect_medfreq_data422(n_424, Fs);

load('labels_n422.mat', 'labels'); 
evaluate(alarm_medfreq_422, labels, 'Median Frequency Data424');

load('labels_n424.mat', 'labels'); 
evaluate(alarm_medfreq_424, labels, 'Median Frequency Data422');


%%

function [alarm,t] = va_detect_medfreq_data422(ecg_data,Fs)
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

    if seg_medfreq  >= 2
     
        alarm(i) = 1;

    end
  
end


end

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


