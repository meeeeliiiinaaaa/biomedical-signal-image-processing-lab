%% Part1

%% 

load("ERP_EEG.mat");

%%

sig = ERP_EEG;
fs = 240;

%% Section1

N_vals = 100:100:2500;
colors = jet(length(N_vals));
time = (0:length(sig(:,1)) - 1) /fs;

figure;
hold on;
for i = 1:length(N_vals)

    N = N_vals(i);
    ERP = mean(sig(:, 1:N), 2);
    plot(time, ERP, 'Color', colors(i,:));

end


xlabel('Time(s)');
ylabel('Amplitude(\muV)');
title('ERP for different number of trials');
legend(arrayfun(@(x) sprintf('N=%d', x), N_vals, 'UniformOutput', false));
grid on;

%% Section2

N_trials = length(sig(1, :));
MaxAbs = zeros(1, N_trials);

for N = 1:N_trials

    ERP = mean(sig(:, 1:N), 2);
    MaxAbs(N) = max(abs(ERP));
end

figure;
plot(1:N_trials, MaxAbs, 'LineWidth', 1, 'Color', 'b');
xlabel('Number of trials(N)');
ylabel('Amplitude(\muV)');
title('Maximum absolute ERP amplitude vs Number of averaged trials');
grid on;

%% Section3

N_trials = length(sig(1, :));
RMS = zeros(1, N_trials);

for N = 2:N_trials

    ERP1 = mean(sig(:, 1:N), 2);
    ERP2 = mean(sig(:, 1:N - 1), 2);
    RMS(N) = sqrt(mean((ERP1 - ERP2).^2));
end

figure;
plot(2:N_trials, RMS(2:end), 'LineWidth', 1, 'Color', 'b');
xlabel('Number of trials(N)');
ylabel('Amplitude(\muV)');
title('RMS of 2 consecutive trials vs Number of averaged trials');
grid on;

%% Section5

N0 = 400;
time = (0:length(sig(:,1)) - 1) /fs;
N_trials = length(sig(1, :));
idx_rand3 = randperm(N_trials,N0);    
idx_rand4 = randperm(N_trials,floor(N0/3));

ERP0 =  mean(sig(:, 1:N0), 2);
ERP1 = mean(sig(:, 1:N_trials), 2);
ERP2 = mean(sig(:, 1:floor(N0/3)), 2);
ERP3 = mean(sig(:, idx_rand3), 2);
ERP4 = mean(sig(:, idx_rand4), 2); 

figure;
hold on;
plot(time, ERP0, 'c--', 'LineWidth', 1);  
plot(time, ERP1, 'k', 'LineWidth', 1.5);          
plot(time, ERP2, 'r--', 'LineWidth', 1);   
plot(time, ERP3, 'b-.', 'LineWidth', 1);    
plot(time, ERP4, 'g-.', 'LineWidth', 1);    

xlabel('Time(s)');
ylabel('Amplitude(\muV)');
title('Comparison of ERP for different number of trials');
legend('First 400 trials','2550 trials','First 133 trials','Random 400 trials', 'Random 133 trials');
grid on;





