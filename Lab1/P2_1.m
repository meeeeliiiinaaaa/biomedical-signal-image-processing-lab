clear; clc; close all;
load("ECG_sig.mat");

N = size(Sig, 1);

t = (0:N-1)/sfreq;
figure;
subplot(2, 1, 1);

plot(t, Sig(:, 1));

title('channel 1');
xlabel('time(s)');

ylabel('amp');
grid minor;

axis tight;

subplot(2, 1, 2);

plot(t, Sig(:, 2));
title('channel 2');

xlabel('time(s)');
ylabel('amp');
grid minor;

axis tight;

%%

indices = [14, 15]; 

for i = 1:2
    beat_time = ATRTIMED(indices(i));

    idx_center = round(beat_time * sfreq);
    
    beat_range = idx_center - round(0.3 * sfreq) : idx_center + round(0.4 * sfreq);
    
    figure('Name', ['Beat ' num2str(i)]);
    plot(t(beat_range), Sig(beat_range, 1), 'LineWidth', 1.5);
    hold on;
    
    text(t(idx_center), Sig(idx_center, 1), 'R', 'Color', 'r', 'FontSize', 12, 'FontWeight', 'bold');

     if i == 1
         
        text(8.72, -0.18, 'P', 'Color', 'r', 'FontSize', 12, 'FontWeight', 'bold');
        text(8.84, -0.27, 'Q', 'Color', 'r', 'FontSize', 12, 'FontWeight', 'bold');
        text(8.9, -0.235, 'S', 'Color', 'r', 'FontSize', 12, 'FontWeight', 'bold');
        text(9.12, -0.1, 'T', 'Color', 'r', 'FontSize', 12, 'FontWeight', 'bold');

    end

    if i==2

        text(9.5, -0.23, 'P', 'Color', 'r', 'FontSize', 12, 'FontWeight', 'bold');
        text(9.59, -0.27, 'Q', 'Color', 'r', 'FontSize', 12, 'FontWeight', 'bold');
        text(9.69, -0.22, 'S', 'Color', 'r', 'FontSize', 12, 'FontWeight', 'bold');
        text(9.88, -0.11, 'T', 'Color', 'r', 'FontSize', 12, 'FontWeight', 'bold');

    end
   

    title(['Sample Beat ' num2str(i) ' - Channel 1']);
    xlabel('Time(s)');
    ylabel('Amplitude');
    grid on;

end





