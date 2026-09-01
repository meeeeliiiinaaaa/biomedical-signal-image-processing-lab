
close all force; clc;
data_path = 'EOG_sig.mat';
load(data_path);

sample_rate = 256;
signal_matrix = Sig;
num_samples = size(signal_matrix, 2);

time_vector = linspace(0, (num_samples - 1) / sample_rate, num_samples);

eog_data = {
    signal_matrix(1, :), 'Left Eye';
    signal_matrix(2, :), 'Right Eye'
};

figure('Name', 'EOG In Time');
for ch = 1:2
    ax = subplot(2, 1, ch);
    plot(ax, time_vector, eog_data{ch, 1});
    title(ax, ['EOG Signal - ', eog_data{ch, 2}]);
    xlabel(ax, 'Time (seconds)');
    ylabel(ax, 'Amplitude');
    grid(ax, 'on');
end

figure('Name', 'Electrode Scheme');
hold on;

theta = linspace(0, 2*pi, 100);
radius = 0.5;
plot(-1 + radius*cos(theta), radius*sin(theta), 'b', 'LineWidth', 1.5);
plot(1 + radius*cos(theta), radius*sin(theta), 'r', 'LineWidth', 1.5); 

elec_pos_x = [-2, 2];
elec_pos_y = [0, 0];
scatter(elec_pos_x, elec_pos_y, 70, 'k', 'filled');

text(elec_pos_x(1) - 0.2, 0, 'left electrode', 'HorizontalAlignment', 'right');
text(elec_pos_x(2) + 0.2, 0, 'right electrode', 'HorizontalAlignment', 'left');

axis equal;
set(gca, 'XLim', [-3 3], 'YLim', [-1 1]);
title('electrode placement in EOG');
grid on;
hold off;
