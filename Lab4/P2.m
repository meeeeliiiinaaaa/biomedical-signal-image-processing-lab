%% Part2

%%

load("SSVEP_EEG.mat");

%%

fs = 250;
sig = SSVEP_Signal;

%% Section1

time = (0:length(sig(1,:)) - 1) / fs;
channels = {'Pz','Oz','P7','P8','O2','O1'};
N_channels = length(channels);

low_cutoff = 1;  
high_cutoff = 40;   
order = 4; 
[b, a] = butter(order, [low_cutoff, high_cutoff] / (fs/2), 'bandpass');

figure;
for ch = 1:N_channels

    sig_filtered = filtfilt(b, a, sig(ch, :));

    subplot(N_channels,1,ch)
    plot(time, sig(ch,:), 'k', 'LineWidth', 1);
    hold on;
    plot(time, sig_filtered, 'r', 'LineWidth', 1.2);
    
    title(['Channel ', channels{ch}]);
    xlabel('Time(S)');
    ylabel('Amplitude(\muV)');
    legend('Raw','Denoised');
    grid on;

end

%% Section2

window_length = fs * 5;
N_events = 15;
N_channels = 6;
trial = zeros(N_channels, window_length, N_events);

for i = 1:N_events

    idx_start = Event_samples(i);
    idx_end = idx_start + window_length - 1;
    trial(:,:,i) = sig(:, (idx_start:idx_end));

end


%% Section3
clc; close all;
N_events = 15;
N_channels = 6;
channels = {'Oz','Pz','P7','P8','O1','O2'};

figure('Units', 'normalized', 'Position', [0.05, 0.05, 0.9, 0.9]); 

for n = 1:N_events
    subplot(3, 5, n); 
    hold on;
    for ch = 1:N_channels
        [Pxx, f] = pwelch(trial(ch,:,n), [], [], [], fs);
        plot(f, 10*log10(Pxx), 'LineWidth', 1);
    end
    
    xlim([0 40]); 
    
    xlabel('Frequency(Hz)');
    ylabel('Power(dB)');
    title(['Trial ', num2str(n)]);
    grid on;
    
    if n == 15
        legend(channels, 'Location', 'bestoutside');
    end
end




