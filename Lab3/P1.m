%% Part1

%% Section1
clc; close all;

load('mecg1.dat');
load('fecg1.dat');
load('noise1.dat');

t_vec = linspace(0, 10, 2560);

total_signal = mecg1 + fecg1 + noise1;


figure('Name', 'Time Domain Signals');
subplot(4, 1, 1); plot(t_vec, mecg1, 'b'); 
title('Maternal ECG (mECG)'); xlabel('Time (s)'); ylabel('Amplitude');

subplot(4, 1, 2); plot(t_vec, fecg1, 'r'); 
title('Fetal ECG (fECG)'); xlabel('Time (s)'); ylabel('Amplitude');

subplot(4, 1, 3); plot(t_vec, noise1, 'k'); 
title('Background Noise'); xlabel('Time (s)'); ylabel('Amplitude');

subplot(4, 1, 4); plot(t_vec, total_signal, 'm'); 
title('Composite ECG Signal'); xlabel('Time (s)'); ylabel('Amplitude');

%% Section2 
clc; close all;

figure('Name', 'PSD Analysis', 'Position', [100, 100, 800, 800]);


t = tiledlayout(2, 2, 'TileSpacing', 'loose', 'Padding', 'normal');

nexttile; 
pwelch(mecg1); 
title('Power Spectral Density: Maternal ECG');

nexttile; 
pwelch(fecg1); 
title('Power Spectral Density: Fetal ECG');

nexttile; 
pwelch(noise1); 
title('Power Spectral Density: Noise Source');

nexttile; 
pwelch(total_signal); 
title('Power Spectral Density: Combined Signal');

%% Section3
clc; close all;

avg_m = mean(mecg1); var_m = var(mecg1);
avg_f = mean(fecg1); var_f = var(fecg1);
avg_n = mean(noise1); var_n = var(noise1);
fprintf('--- Statistical Properties (Mean & Variance) ---\n');
fprintf('Maternal Signal => Average: %.10f | Variance: %.10f\n', avg_m, var_m);
fprintf('Fetal Signal    => Average: %.10f | Variance: %.10f\n', avg_f, var_f);
fprintf('Noise Data      => Average: %.10f | Variance: %.10f\n', avg_n, var_n);
fprintf('----------------------------------------------\n');

%% Section4
clc; close all;

figure;
subplot(3, 1, 1); hist(mecg1, 50); title('Histogram of Maternal ECG Signal');
subplot(3, 1, 2); hist(fecg1, 50); title('Histogram of Fetal ECG Signal');
subplot(3, 1, 3); hist(noise1, 50); title('Histogram of Noise Signal');

kurtosis_mceg = kurtosis(mecg1);
kurtosis_fecg = kurtosis(fecg1);
kurtosis_noise = kurtosis(noise1);


fprintf('Kurtosis:\n');
fprintf('Maternal ECG: %.4f\n', kurtosis_mceg);
fprintf('Fetal ECG: %.4f\n', kurtosis_fecg);
fprintf('Noise: %.4f\n', kurtosis_noise);

