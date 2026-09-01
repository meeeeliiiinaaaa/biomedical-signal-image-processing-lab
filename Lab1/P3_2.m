clc; close all;
data_path = 'EOG_sig.mat';
load(data_path);
left_eye = Sig(1, :); 
right_eye = Sig(2, :); 


sig_length = length(left_eye);
freq_axis = (0:sig_length-1) * (fs / sig_length);

channels_info = {
    left_eye,  'Left Eye';
    right_eye, 'Right Eye'
};


fig_fft = figure('Name', 'Frequency Domain');
for i = 1:2
    ax = subplot(2, 1, i);
    dft_mag = abs(fft(channels_info{i, 1}));
    plot(ax, freq_axis, dft_mag, 'LineWidth', 1.2);
    set(ax, 'XLim', [0 60]);
    title(ax, ['Frequency Spectrum - ', channels_info{i, 2}]);
    xlabel(ax, 'Frequency (Hz)');
    ylabel(ax, 'Magnitude');
    grid("on");
end

fig_psd = figure('Name', 'Power Spectral Density');

safe_win_psd = min(128, sig_length); 
safe_overlap_psd = floor(safe_win_psd / 2);

for i = 1:2
    subplot(2, 1, i);
    pwelch(channels_info{i, 1}, safe_win_psd, safe_overlap_psd, safe_win_psd, fs);
    title(['Power Spectral Density - ', channels_info{i, 2}]);
end

fig_stft = figure('Name', 'Spectrogram Analysis');
stft_win = min(128, sig_length);
stft_overlap = min(64, floor(stft_win/2));
stft_nfft = 128;

for i = 1:2
    ax_stft = subplot(2, 1, i);
    spectrogram(channels_info{i, 1}, stft_win, stft_overlap, stft_nfft, fs, 'yaxis');
    set(ax_stft, 'YLim', [0 60]);
    title(ax_stft, ['Time-Frequency Spectrum (STFT) - ', channels_info{i, 2}]);
end
