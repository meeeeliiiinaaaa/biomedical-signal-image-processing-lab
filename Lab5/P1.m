%% Part1

%% Section1
clearvars; clc; close all;

data_struct = load('normal.mat');
raw_data = data_struct.normal; 

t_vec = raw_data(:, 1)';
sig_vec = raw_data(:, 2)';

idx_healthy = 251:2500;
idx_faulty = 72851:75100;

t_healthy = t_vec(idx_healthy);
sig_healthy = sig_vec(idx_healthy);

t_faulty = t_vec(idx_faulty);
sig_faulty = sig_vec(idx_faulty);


fs = 250;
[psd_healthy, freq_h] = pwelch(sig_healthy, [], [], [], fs);
[psd_faulty, freq_f] = pwelch(sig_faulty, [], [], [], fs);


figure('Name', 'Initial Spectra Analysis');
plot(freq_h, 20*log10(psd_healthy), 'LineWidth', 1.5);
hold on;
plot(freq_f, 20*log10(psd_faulty), 'LineWidth', 1.5);
grid minor;
xlabel('Frequency [Hz]', 'FontWeight', 'bold');
ylabel('Magnitude [dB]', 'FontWeight', 'bold');
title('Spectral Comparison: Clean VS. Noisy ECG Signal');
legend('Clean', 'Noisy');
axis tight;
hold off;

%% Section2

E_total = sum(psd_healthy(freq_h > 0.5));

for upper_f = 1:0.5:250
    E_current = sum(psd_healthy((upper_f > freq_h) & (freq_h > 0.5)));
    if E_current > 0.99 * E_total
        break;
    end
end

fprintf('>> Calculated Upper Cutoff Frequency: %.2f Hz\n', upper_f);

N_ord = 4;
[num_coeff, den_coeff] = butter(N_ord, [0.004, upper_f/125], 'bandpass');

[H_resp, w_freq] = freqz(num_coeff, den_coeff, 1024, fs);
figure('Name', 'Bode Plot');
plot(w_freq, 20*log10(abs(H_resp)), 'LineWidth', 1.2, 'Color', '#0072BD');
title('Bode Diagram of the Band-pass Filter');
xlabel('Frequency [Hz]');
ylabel('Gain [dB]');
grid minor;
axis tight;

h_impulse = impz(num_coeff, den_coeff);
figure('Name', 'Impulse Behavior');
stem(h_impulse, 'filled', 'MarkerSize', 4, 'Color', '#D95319');
title('Filter Impulse Behavior');
xlabel('Sample Index');
ylabel('Amplitude [V]');
grid minor;

filtered_full_sig = bandpass(sig_vec, [0.004, upper_f/125]);

filt_sig_healthy = filtered_full_sig(idx_healthy);
filt_sig_faulty = filtered_full_sig(idx_faulty);


[psd_filt_h, f_filt_h] = pwelch(filt_sig_healthy, [], [], [], fs);
[psd_filt_f, f_filt_f] = pwelch(filt_sig_faulty, [], [], [], fs);

figure('Name', 'Filtered Spectra');
plot(f_filt_h, 20*log10(psd_filt_h), 'LineWidth', 1.5);
hold on;
plot(f_filt_f, 20*log10(psd_filt_f), 'LineWidth', 1.5, 'LineStyle', '-.');
xlabel('Frequency [Hz]', 'FontWeight', 'bold');
ylabel('Magnitude [dB]', 'FontWeight', 'bold');
title('Power Spectrum of Processed Signals');
legend('Clean (Filtered)', 'Noisy (Filtered)');
axis tight;
grid minor;
hold off;

%% Section3

figure('Name', 'Time Domain - Healthy');
subplot(2, 1, 1);
plot(t_healthy, sig_healthy, 'Color', '#0072BD');
title('Original Clean Data (Time-Domain)');
xlabel('Time [sec]');
ylabel('Voltage [V]');
axis tight;
grid minor;

subplot(2, 1, 2);
plot(t_healthy, filt_sig_healthy, 'Color', '#D95319');
title('Processed Clean Data (Time-Domain)');
xlabel('Time [sec]');
ylabel('Voltage [V]');
axis tight;
grid minor;

figure('Name', 'Time Domain - Faulty');
subplot(2, 1, 1);
plot(t_faulty, sig_faulty, 'Color', '#0072BD');
title('Original Noisy Data (Time-Domain)');
xlabel('Time [sec]');
ylabel('Voltage [V]');
axis tight;
grid minor;

subplot(2, 1, 2);
plot(t_faulty, filt_sig_faulty, 'Color', '#D95319');
title('Processed Noisy Data (Time-Domain)');
xlabel('Time [sec]');
ylabel('Voltage [V]');
axis tight;
grid minor;
