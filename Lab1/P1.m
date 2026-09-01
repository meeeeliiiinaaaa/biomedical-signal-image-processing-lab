%% PART 1
%% Q1

% load
load("EEG_sig.mat");   

% properties
fs = des.samplingfreq;   
channel_labels = des.channelnames;

% signal
channel_idx = 5;           % channel 5             
signal = Z(channel_idx, :);     

% time axis
N = length(signal);        
t = (0:N-1) / fs;              

% plot
figure;
set(gcf, 'Position', [80, 80, 1200, 400]);
plot(t, signal, 'b-', 'LineWidth', 0.5);
xlabel('Time (s)');
ylabel('Amplitude (\muV)');  
title(sprintf('EEG Channel %d (%s)', channel_idx, channel_labels{channel_idx}));
grid on;

%% Q2

% intervals
intervals = {
    '0–15 s',             [0, 15];
    '18–40 s',            [18, 40];
    '45–50 s',            [45, 50];
    '50 s to end',        [50, t(end)];
};


figure('Position', [100, 100, 1200, 600]);

for i = 1:length(intervals)
   
    t_range = intervals{i, 2};      % ex: i=1 --> [0 15]
    idx = t >= t_range(1) & t <= t_range(2);
    seg_signal = signal(idx);
    seg_time = t(idx);
    
    % Plot
    subplot(2, 2, i);
    plot(seg_time, seg_signal, 'b-', 'LineWidth', 0.5);
    xlabel('Time (s)');
    ylabel('Amplitude (\muV)');
    title(intervals{i, 1});
    grid on;
    xlim([t_range(1), t_range(2)]);
  
    
    % Compute properties
    mean_val = mean(seg_signal);
    std_val = std(seg_signal);
    peak2peak = max(seg_signal) - min(seg_signal);
       

    % Display on plot
    text_str = sprintf('Mean: %.2f µV\nStd: %.2f µV\nPk2Pk: %.2f', mean_val, std_val, peak2peak);
    text(0.02, 0.98, text_str, 'Units', 'normalized', 'VerticalAlignment', 'top','BackgroundColor', 'white', 'FontSize', 8);
end

sgtitle('EEG Channel 5 (C3) in 4 Intervals');

%% Q3

% signal
channel_idx2 = 17;           % channel 17    
signal2 = Z(channel_idx2, :);     

% time axis
N2 = length(signal2);        
t2 = (0:N2-1) / fs;              

% plot
figure;
set(gcf, 'Position', [80, 80, 1200, 400]);
plot(t2, signal2, 'b-', 'LineWidth', 0.5);
xlabel('Time (s)');
ylabel('Amplitude (\muV)');  
title(sprintf('EEG Channel %d (%s)', channel_idx2, channel_labels{channel_idx2}));
grid on;


%% Q4

offset = max(max(abs(Z)))/5 ; 
feq = 256 ; 
ElecName = des.channelnames ; 
disp_eeg(Z,offset,feq,ElecName) ; 

%% Q6

% intervals
intervals = {'2–7 s',    [2, 7];     
             '30–35 s',  [30, 35];   
             '42–47 s',  [42, 47];   
             '50–55 s',  [50, 55]};  

for i = 1:length(intervals)
   
    t_range = intervals{i, 2};                % ex: i=1 --> t_range=[2 7]
    idx = t >= t_range(1) & t <= t_range(2);
    seg_signal = signal(idx);
    seg_time = t(idx);
  
    
    % plot - time domain
    figure;
    subplot(2,1,1);
    plot(seg_time, seg_signal, 'b-', 'LineWidth', 0.7);
    xlabel('Time (s)');
    ylabel('Amplitude (\muV)');
    title(sprintf('EEG channel 5 (C3) : %s', intervals{i, 1}));
    grid on;
    xlim([t_range(1), t_range(2)]);
    
    % DFT
    N = length(seg_signal);

    % FFT (magnitude)
    Y = fft(seg_signal);
    P2 = abs(Y / N);             
    P1 = P2(1:floor(N/2)+1);       
    P1(2:end-1) = 2 * P1(2:end-1);
    
    % frequency axis
    f = fs * (0:(N/2)) / N;        
    
    % limit frequency (60 HZ)
    max_freq = 60;               
    idx_f = f <= max_freq;
    seg_f = f(idx_f);
    seg_P1 = P1(idx_f);
    
    % plot - DFT
    subplot(2,1,2);
    plot(seg_f, seg_P1, 'r-', 'LineWidth', 0.7);
    xlabel('Frequency (Hz)');
    ylabel('Magnitude');
    title(sprintf('Frequency spectrum (DFT)'));
    xlim([0, max_freq]);
    grid on;
    
end


%% Q7

% intervals
intervals = {'2–7 s',    [2, 7];     
             '30–35 s',  [30, 35];   
             '42–47 s',  [42, 47];  
             '50–55 s',  [50, 55]};  


% pwelch parameters
window = 4 * fs;                 
noverlap = round(0.5 * window);  
nfft = 4 * fs;                    

% frequency range 
freq_range = [0, 60];

for i = 1:length(intervals)

    t_range = intervals{i, 2};             % ex: i=1 --> t_range=[2 7]
    idx = t >= t_range(1) & t <= t_range(2);
    seg_signal = signal(idx);
    
    % compute PSD
    [pxx, f] = pwelch(seg_signal, hamming(window), noverlap, nfft, fs);
    
    % Limit frequency
    idx_f = f <= freq_range(2);
    seg_f = f(idx_f);
    seg_pxx = pxx(idx_f);
    
    % Plot - PSD
    figure;
    plot(seg_f, 10*log10(seg_pxx), 'b-', 'LineWidth', 1);
    xlabel('Frequencyseg_ (Hz)');
    ylabel('Power Spectral Density (dB/Hz)');
    title(sprintf('EEG channel 5 (C3) - PSD : %s', intervals{i, 1}));
    grid on;
    xlim(freq_range);
 
    hold off;
end

%% Q8

% intervals
intervals = {'2–7 s',    [2, 7];     
             '30–35 s',  [30, 35];  
             '42–47 s',  [42, 47];  
             '50–55 s',  [50, 55]};  

% spectrogram parameters 
L = 128;               
n_overlap = 64;          
nfft = L;                
window = hamming(L);     


for i = 1:length(intervals)

    t_range = intervals{i, 2};            % ex: i=1 --> t_range=[2 7]
    idx = t >= t_range(1) & t <= t_range(2);
    seg_signal = signal(idx);
    
    % compute spectrogram
    [S, F, T] = spectrogram(seg_signal, window, n_overlap, nfft, fs);
    
    % Plot - spectrogram 
    figure;
    imagesc(T, F, 10*log10(abs(S) + eps)); 
    axis xy;                                
    xlabel('Time (s)');
    ylabel('Frequency (Hz)');
    title(sprintf('EEG channel 5 (C3) - Spectrogram : %s', intervals{i, 1}));
    colorbar;
    colormap('jet');

end

%% Q9

fs_old = des.samplingfreq;             
channel_idx = 5;              % channel 5
signal = Z(channel_idx, :); 
N = length(signal);
t = (0:N-1) / fs_old;

% second interval: 30–35 seconds
t_start = 30;
t_end   = 35;
idx = t >= t_start & t <= t_end;
seg_sig_old = signal(idx);
seg_time_old = t(idx);

% --------Old signal-----------------

% plot - time
figure('Position', [50 50 1200 600]);

subplot(3,2,1);
plot(seg_time_old, seg_sig_old, 'b');
xlabel('Time (s)');
ylabel('Amplitude (\muV)');
title('Original EEG channel 5 (C3): 30 - 35 s');
xlim([t_start t_end]);
grid on;

% plot - DFT
N = length(seg_sig_old);
Y = fft(seg_sig_old);
P2 = abs(Y/N);
P1 = P2(1:floor(N/2)+1);
P1(2:end-1) = 2*P1(2:end-1);

f = fs_old * (0:(N/2)) / N;

max_freq = 60;  
idx_f = f <= max_freq;

subplot(3,2,3);
plot(f(idx_f), P1(idx_f), 'r');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title(sprintf('Original DFT (fs=%.1f Hz)', fs_old));
xlim([0 max_freq]);
grid on; 

% plot STFT
L = 128;                    
n_overlap = 64;
nfft = L;
window = hamming(L);

[S_old, F_old, T_old] = spectrogram(seg_sig_old, window, n_overlap, nfft, fs_old);

subplot(3,2,5);
imagesc(T_old + t_start, F_old, 10*log10(abs(S_old) + eps));
axis xy;
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title(sprintf('Original STFT(fs=%.1f Hz)', fs_old));
colorbar;
colormap('jet');

% --------Down sampled signal------------

down_factor = floor(fs_old / 100);   
if down_factor < 2
    down_factor = 2;                 
end

fs_new = fs_old / down_factor;

fprintf('Original fs = %.1f Hz, downsampling factor = %d, new fs = %.1f Hz\n', fs_old, down_factor, fs_new);

% design low pass filter
cutoff = fs_new / 2;              
trans_width = 5;                    
astop = 60;                        
order = firpmord([cutoff, cutoff+trans_width], [1 0], [0.01, 10^(-astop/20)], fs_old);

if order < 10, order = 30; end    
lp_filter = firpm(order, [0 cutoff cutoff+trans_width fs_old/2] / (fs_old/2), [1 1 0 0], [1 10^(-astop/20)]);


% Apply filter 
seg_filtered = filtfilt(lp_filter, 1, seg_sig_old);

% Downsample
seg_sig_down = resample(seg_filtered, 1, down_factor); 
seg_time_down = linspace(t_start, t_end, length(seg_sig_down));


% plot - time
subplot(3,2,2);
plot(seg_time_down, seg_sig_down, 'b');
xlabel('Time (s)');
ylabel('Amplitude (\muV)');
title(sprintf('Downsampled EEG channel 5 (C3): 30 - 35 s (M= %.1f)', down_factor));
xlim([t_start t_end]);
grid on; 

% plot - DFT
N2 = length(seg_sig_down);
Y2 = fft(seg_sig_down);
P2b = abs(Y2/N2);
P1b = P2b(1:floor(N2/2)+1);
P1b(2:end-1) = 2*P1b(2:end-1);

f2 = fs_new * (0:(N2/2)) / N2;
idx_f2 = f2 <= max_freq;

subplot(3,2,4);
plot(f2(idx_f2), P1b(idx_f2), 'r');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
title(sprintf('Downsampled DFT (fs=%.1f Hz)', fs_new));
xlim([0 max_freq]);
grid on; 

% plot - STFT

win_duration = L / fs_old;            
L_new = round(win_duration * fs_new);  
n_overlap_new = round(0.5 * L_new);     
nfft_new = L_new;
window_new = hamming(L_new);

[S_down, F_down, T_down] = spectrogram(seg_sig_down, window_new, n_overlap_new, nfft_new, fs_new);

subplot(3,2,6);
imagesc(T_down + t_start, F_down, 10*log10(abs(S_down) + eps));
axis xy;
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title(sprintf('Downsampled STFT (fs=%.1f Hz)', fs_new));
colorbar;
colormap('jet');

sgtitle('Original vs Downsampled EEG');


%% disp_eeg function

function t = disp_eeg(X,offset,feq,ElecName,titre)
% function t = disp_eeg(X,offset,feq,ElecName,titre)
%
% inputs
%     X: dynamics to display. (nbchannels x nbsamples) matrix
%     offset: offset between channels (default max(abs(X)))
%     feq: sapling frequency (default 1)
%     ElecName: cell array of electrode labels (default {S1,S2,...})
%     titre: title of the figure
%
% output
%     t: time vector
%
% G. Birot 2010-02


%% Check arguments
[N K] = size(X);

if nargin < 4
    for n = 1:N
        ElecName{n}  = ['S',num2str(n)];
    end
    titre = [];
end

if nargin < 5
    titre = [];
end

if isempty(feq)
    feq = 1;
end

if isempty(ElecName)
    for n = 1:N
        ElecName{n}  = ['S',num2str(n)];
    end
end

if isempty(offset)
    offset = max(abs(X(:)));
end


%% Build dynamic matrix with offset and time vector
X = X + repmat(offset*(0:-1:-(N-1))',1,K);
t = (1:K)/feq;
graduations = offset*(0:-1:-(N-1))';
shiftvec = N:-1:1;
Ysup = max(X(1,:)) + offset;
Yinf = min(X(end,:)) - offset;
% YLabels = cell(N+2) ElecName(shiftvec)

%% Display
figure1 = figure;
% a1 = axes('YAxisLocation','right');
a2 = axes('YTickLabel',ElecName(shiftvec),'YTick',graduations(shiftvec),'FontSize',7);
ylim([Yinf Ysup]);
box('on');
grid('on')
hold('all');
plot(t,X');
xlabel('Time (seconds)','FontSize',10);
ylabel('Channels','FontSize',10);
title(titre);
hold off

end






























