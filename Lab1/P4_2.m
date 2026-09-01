
clearvars; close all force; clc;
load("EMG_sig.mat"); 
samp_rate = 4000; 

healthy_signal = emg_healthym;     
neuropathy_signal = emg_neuropathym;   
myopathy_signal = emg_myopathym;   


analysis_data = {
    healthy_signal,    'healthy';
    neuropathy_signal, 'neuropathy';
    myopathy_signal,   'myopathy'
};

num_classes = size(analysis_data, 1);

figure('Name', 'FFT', 'NumberTitle', 'off');
for k = 1:num_classes
    curr_sig = analysis_data{k, 1};
    sig_len = length(curr_sig);
    f_axis = (0:sig_len-1) * (fs / sig_len);
    fft_mag = abs(fft(curr_sig));  
    ax_fft = subplot(num_classes, 1, k);
    plot(ax_fft, f_axis, fft_mag);
    set(ax_fft, 'XLim', [0 1000]);
    title(ax_fft, ['frequency spectrum  ', analysis_data{k, 2}]);
    xlabel(ax_fft, 'Frequency (Hz)');
    ylabel(ax_fft, 'Magnitude');
end

figure('Name', 'PSD', 'NumberTitle', 'off');
for k = 1:num_classes
    subplot(num_classes, 1, k);
    pwelch(analysis_data{k, 1}, [], [], [], fs);
    title(['power spectral density  ', analysis_data{k, 2}]);
end

stft_cfg = struct('win', 128, 'overlap', 64, 'nfft', 128);
figure('Name', 'STFT', 'NumberTitle', 'off');
for k = 1:num_classes
    ax_stft = subplot(num_classes, 1, k); 
    spectrogram(analysis_data{k, 1}, stft_cfg.win, stft_cfg.overlap, stft_cfg.nfft, fs, 'yaxis'); 
    set(ax_stft, 'YLim', [0 2]);
    title(ax_stft, ['time-frequency spectrum  ', analysis_data{k, 2}]);
end
