clc; clearvars; close all;

rawData = load('X_org.mat');
signalData = rawData.X_org;

visualize_brain_signals(signalData);

function visualize_brain_signals(data_matrix)
    electrode_info = load('Electrodes');
    fs = 250; 
    spacing = max(abs(data_matrix(:))); % Offset between signals
    ch_labels = electrode_info.Electrodes.labels;
    render_eeg_traces(data_matrix, spacing, fs, ch_labels);
end

function time_vec = render_eeg_traces(data_mat, spacing, fs, ch_labels, fig_title)
    [num_channels, num_samples] = size(data_mat);
    if ~exist('ch_labels', 'var') || isempty(ch_labels)
        ch_labels = arrayfun(@(x) sprintf('Ch-%d', x), 1:num_channels, 'UniformOutput', false);
    end
    
    if ~exist('fig_title', 'var')
        fig_title = 'EEG Signal Overview';
    end
    
    if ~exist('fs', 'var') || isempty(fs)
        fs = 1;
    end
    
    if ~exist('spacing', 'var') || isempty(spacing)
        spacing = max(abs(data_mat(:)));
    end
    

    offset_multiplier = (0:-1:-(num_channels-1))';
    channel_offsets = spacing * offset_multiplier;

    data_offset = data_mat + (channel_offsets * ones(1, num_samples));

    time_vec = (1:num_samples) / fs;
    

    upper_bound = max(data_offset(1,:)) + (spacing * 1.1);
    lower_bound = min(data_offset(end,:)) - (spacing * 1.1);
    reverse_idx = num_channels:-1:1;

    eeg_fig = figure('Color', 'w', 'Name', 'EEG Viewer');
    ax = axes('Parent', eeg_fig);
    hold(ax, 'on');
    

    plot(ax, time_vec, data_offset', 'LineWidth', 0.8);

    set(ax, 'YTick', channel_offsets(reverse_idx), ...
            'YTickLabel', ch_labels(reverse_idx), ...
            'FontSize', 8, ...
            'XGrid', 'on', 'YGrid', 'on', ...
            'GridLineStyle', ':', ...
            'GridAlpha', 0.7);
            
    ylim(ax, [lower_bound upper_bound]);
    xlim(ax, [time_vec(1) time_vec(end)]);
    box(ax, 'on');

    xlabel(ax, 'Time [s]', 'FontWeight', 'bold', 'FontSize', 11);
    ylabel(ax, 'Electrode', 'FontWeight', 'bold', 'FontSize', 11);
    
    if ~isempty(fig_title)
        title(ax, fig_title, 'FontSize', 12, 'FontWeight', 'bold');
    end
    
    hold(ax, 'off');
end
