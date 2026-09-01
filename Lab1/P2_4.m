
clearvars; clc; close all;
data_path = 'ECG_sig.mat';
load(data_path);
idx_healthy = find(ANNOTD == 1);
idx_abnormal = find(ANNOTD ~= 1);

for k = 2:(numel(idx_healthy)-2)
    if all(diff(idx_healthy(k:k+2)) == 1)
        idx_start = floor(ATRTIMED(idx_healthy(k)-1) * sfreq);
        idx_end = floor(ATRTIMED(idx_healthy(k)+3) * sfreq);
        
        time_vec_norm = idx_start:idx_end;
        sig_healthy = Sig(idx_start:idx_end, :);
        
        figure('Name', 'healthy ECG');
        ax1 = subplot(2,1,1);
        plot(ax1, time_vec_norm / sfreq, sig_healthy(:, 1));
        title(ax1, 'Normal Signal Ch_1');
        xlabel(ax1, 'Time (s)'); ylabel(ax1, 'Amp');
        grid(ax1, 'minor'); axis(ax1, 'tight');
        
        ax2 = subplot(2,1,2);
        plot(ax2, time_vec_norm / sfreq, sig_healthy(:, 2));
        title(ax2, 'Normal Signal Ch_2');
        xlabel(ax2, 'Time (s)'); ylabel(ax2, 'Amp');
        grid(ax2, 'minor'); axis(ax2, 'tight');
        break;
    end
end

for m = 2:(length(idx_abnormal)-2)
    if all(diff(idx_abnormal(m:m+2)) == 1)
        idx_start_abn = floor(ATRTIMED(idx_abnormal(m)-1) * sfreq);
        idx_end_abn = floor(ATRTIMED(idx_abnormal(m)+3) * sfreq);
        
        time_vec_abn = idx_start_abn:idx_end_abn;
        sig_abnormal = Sig(idx_start_abn:idx_end_abn, :);
        
        figure('Name', 'Abnormal ECG');
        subplot(2,1,1);
        plot(time_vec_abn / sfreq, sig_abnormal(:, 1));
        title('Abnormal Signal Ch_1');
        xlabel('Time (s)'); ylabel('Amp');
        grid minor; axis tight;
        
        subplot(2,1,2);
        plot(time_vec_abn / sfreq, sig_abnormal(:, 2));
        title('Abnormal Signal Ch_2');
        xlabel('Time (s)'); ylabel('Amp');
        grid minor; axis tight;
        break;
    end
end

plot_data = {
    sig_healthy(:, 1),  'Normal (Ch 1)';
    sig_healthy(:, 2),  'Normal (Ch 2)';
    sig_abnormal(:, 1), 'Anomalous (Ch 1)';
    sig_abnormal(:, 2), 'Anomalous (Ch 2)'
};
figure('Name', 'Frequency Domains');
for j = 1:4
    subplot(2, 2, j);
    calc_and_plot_fft(plot_data{j, 1}, sfreq);
    title(['Freq Spectrum - ', plot_data{j, 2}]);
end

figure('Name', 'Time-Frequency Domains');
for j = 1:4
    subplot(2, 2, j);
    calc_and_plot_stft(plot_data{j, 1}, sfreq);
    title(['STFT - ', plot_data{j, 2}]);
end

function calc_and_plot_fft(data_seq, fs)
    N = length(data_seq);
    fft_res = fft(data_seq);
    magnitude = abs(fft_res / N);
    
    half_mag = magnitude(1:floor(N/2)+1);
    half_mag(2:end-1) = 2 * half_mag(2:end-1);
    
    freq_axis = fs * (0:floor(N/2)) / N;
    
    plot(freq_axis, half_mag);
    xlabel('Frequency (hz)');
    ylabel('Magnitude');
end

function calc_and_plot_stft(data_seq, fs)
    stft(data_seq, fs, 'Window', hamming(128), 'OverlapLength', 120, 'FFTLength', 256);
    colormap jet;
    colorbar;
end
