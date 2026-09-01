clc; clearvars; close all;
load('ElecPosXYZ.mat');

head_radii = [8, 8.5, 9.2];
conductivities = [3.3e-3, 8.25e-5, 3.3e-3];
model_struct.R = head_radii;
model_struct.Sigma = conductivities;
model_struct.Lambda = [.5979, .2037, .0237];
model_struct.Mu = [.6342, .9364, 1.0362];
grid_res = 1;

%% Part (Alef)
[source_locs, lead_field] = ForwardModel_3shell(grid_res, model_struct);

figure('Color', 'w');
scatter3(source_locs(1,:), source_locs(2,:), source_locs(3,:), 15, 'filled', 'MarkerFaceColor', [0.6 0.6 0.6]);
title('3D Source Space Grid Visualization');
xlabel('X-axis'); ylabel('Y-axis'); zlabel('Z-axis');
axis equal; grid on;
save('LeadFieldData.mat', 'lead_field');

%% Part (Be)
num_elec = 21;
elec_coords = zeros(3, num_elec);
elec_labels = cell(1, num_elec);

for k = 1:num_elec 
    elec_labels{k} = ElecPos{1,k}.Name;
    elec_coords(:,k) = ElecPos{1,k}.XYZ;
end

surface_radius = 9.2;
scaled_elec_coords = elec_coords * surface_radius;

figure('Color', 'w');
scatter3(source_locs(1,:), source_locs(2,:), source_locs(3,:), 10, 'b');
hold on;
scatter3(scaled_elec_coords(1,:), scaled_elec_coords(2,:), scaled_elec_coords(3,:), 45, 'r', 'filled');
text(scaled_elec_coords(1,:) + 0.3, scaled_elec_coords(2,:) + 0.3, scaled_elec_coords(3,:) + 0.3, elec_labels, 'FontSize', 9);
title('Sensor and Source Locations Overlay');
axis equal; grid on;

%% (Pe) to (Khe) 

test_indices = [100, 200, 500, 1000];
load('Interictal.mat');
signal_wave = Interictal(10,:);

for iter = 1:length(test_indices)
    current_id = test_indices(iter);
    figure('Color', 'w');
    scatter3(source_locs(1,:), source_locs(2,:), source_locs(3,:), 10, [0.75 0.75 0.75]);
    hold on;
    scatter3(scaled_elec_coords(1,:), scaled_elec_coords(2,:), scaled_elec_coords(3,:), 35, 'r', 'filled');
    text(scaled_elec_coords(1,:) + 0.3, scaled_elec_coords(2,:) + 0.3, scaled_elec_coords(3,:) + 0.3, elec_labels);
    
    target_pos = source_locs(:, current_id);
    scatter3(target_pos(1), target_pos(2), target_pos(3), 70, 'g', 'filled');
    
    mag_pos = norm(target_pos);
    radial_vec = target_pos / mag_pos;
    end_pt = target_pos + radial_vec;
    plot3([target_pos(1), end_pt(1)], [target_pos(2), end_pt(2)], [target_pos(3), end_pt(3)], 'g', 'LineWidth', 2.5);
    title(['Scenario ', num2str(iter), ' | Active Source Index: ', num2str(current_id)]);

    dipole_orientation = radial_vec; 
    cols = (current_id-1)*3 + 1 : current_id*3;
    eeg_matrix = lead_field(:, cols) * dipole_orientation * signal_wave;
    disp_eeg(eeg_matrix, max(abs(eeg_matrix(:))), [], elec_labels, sprintf('Simulated EEG (Node %d)', current_id));

    avg_spikes = zeros(21, 1);
    for ch = 1:21
        [~, peak_times] = findpeaks(eeg_matrix(ch,:), 'MinPeakHeight', 0.5 * max(eeg_matrix(ch,:)));
        offsets = -3:3;
        idx_matrix = peak_times' + offsets; 
        valid_indices = idx_matrix(:);
        avg_spikes(ch) = mean(eeg_matrix(ch, valid_indices));
    end
    
    figure;
    Display_Potential_3D(9.2, avg_spikes);

    reg_param = 0.1;
    I_mat = eye(size(lead_field, 1));
    w_mne = lead_field' * ((lead_field * lead_field' + reg_param * I_mat) \ avg_spikes);
    
    mne_reshaped = reshape(w_mne, 3, [])';
    source_magnitudes = sqrt(sum(mne_reshaped.^2, 2));
    [~, best_guess_idx] = max(source_magnitudes);

    guess_pos = source_locs(:, best_guess_idx);
    
    fprintf('----------------------------------------\n');
    fprintf('>> MNE RESULTS FOR EXPERIMENT %d \n', iter);
    fprintf('----------------------------------------\n');
    fprintf('True Source Node ID : %d\n', current_id);
    fprintf('Estimated Source ID : %d\n', best_guess_idx);
    
    pos_err_cm = norm(guess_pos - target_pos);
    
    est_dir = mne_reshaped(best_guess_idx, :)';
    est_dir = est_dir / norm(est_dir);
    
    dir_err_deg = acosd(dot(radial_vec, est_dir) / (norm(radial_vec) * norm(est_dir)));
    
    fprintf('Spatial Localization Error : %.3f cm\n', pos_err_cm);
    fprintf('Angular Orientation Error  : %.3f degrees\n\n', dir_err_deg);
end

%% Part (Dal)
fprintf('\n>>> Initializing Genetic Algorithm Search...\n');
true_moment = [1, 0, 0]; 
inter_signal = Interictal(1, :); 
total_nodes = size(source_locs, 2); 
random_target = randi(total_nodes); 
cols_true = (random_target-1)*3 + 1 : random_target*3;
true_potentials = lead_field(:, cols_true) * true_moment' * inter_signal;

ga_opts = optimoptions('ga', 'Display', 'iter', 'PopulationSize', 80, ...
    'MaxGenerations', 80, 'UseParallel', true);

bounds_lower = [1, -1, -1, -1];
bounds_upper = [total_nodes, 1, 1, 1];

fitness_fcn = @(x) evaluate_ga_cost(x, lead_field, true_potentials(:), total_nodes, inter_signal);

[optimal_params, min_cost] = ga(fitness_fcn, 4, [], [], [], [], bounds_lower, bounds_upper, [], ga_opts);

best_node_ga = round(optimal_params(1));
best_moment_ga = optimal_params(2:4);

fprintf('\n====== GA OPTIMIZATION SUMMARY ======\n');
fprintf('Identified Source Node: %d\n', best_node_ga);
fprintf('Coordinates Discovered: [%.2f, %.2f, %.2f]\n', source_locs(:, best_node_ga));
fprintf('Moment Vector Est: [%.2f, %.2f, %.2f]\n', best_moment_ga);
fprintf('Final Residual Cost: %.6f\n', min_cost);


function err_val = evaluate_ga_cost(vars, LFM, target_pot, max_nodes, time_signal)
    node_idx = round(vars(1));
    node_idx = max(1, min(max_nodes, node_idx)); 
    moment_vec = vars(2:4)';
    cols = (node_idx-1)*3 + 1 : node_idx*3;
    sim_pot = LFM(:, cols) * moment_vec * time_signal;
    err_val = norm(sim_pot(:) - target_pot(:));
end