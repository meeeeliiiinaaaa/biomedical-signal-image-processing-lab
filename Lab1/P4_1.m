
clearvars; close all force; clc;
data_in = load("C:\Users\melin\Downloads\codes\codes\EMG_sig.mat"); 
samp_rate = 4000; 

emg_collection = {
    data_in.emg_healthym,     'EMG signal healthy';
    data_in.emg_neuropathym,  'EMG signal neuropathy';
    data_in.emg_myopathym,    'EMG signal myopathy'
};

num_signals = size(emg_collection, 1);

figure('Name', 'full scale EMG', 'NumberTitle', 'off');
for i = 1:num_signals
    current_sig = emg_collection{i, 1};
    time_ax = (0:length(current_sig)-1) / samp_rate;    
    ax1 = subplot(num_signals, 1, i);
    plot(ax1, time_ax, current_sig);
    title(ax1, emg_collection{i, 2});
    xlabel(ax1, 'Time(sec)');
    ylabel(ax1, 'Amp');
    grid(ax1, 'on');
end

figure('Name', 'zoomed EMG', 'NumberTitle', 'off');
for i = 1:num_signals
    current_sig = emg_collection{i, 1};
    time_ax = (0:length(current_sig)-1) / samp_rate;  
    ax2 = subplot(num_signals, 1, i);
    plot(ax2, time_ax, current_sig);
    set(ax2, 'XLim', [0.5 1]); 
    title(ax2, ['Zoomed ', emg_collection{i, 2}]);
    xlabel(ax2, 'Time(sec)');
    ylabel(ax2, 'Amp');
    grid(ax2, 'on');
end
