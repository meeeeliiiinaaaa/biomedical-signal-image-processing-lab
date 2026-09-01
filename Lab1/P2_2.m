clear; clc; close all;
load("ECG_sig.mat");

numSamples = size(Sig, 1);
timeAxis = (0:numSamples-1) ./ sfreq;

annoCodes = [0:14, 16, 18:41];
annoNames = {'NOTQRS', 'NORMAL', 'LBBB', 'RBBB', 'ABERR', 'PVC', 'FUSION', ...
             'NPC', 'APC', 'SVPB', 'VESC', 'NESC', 'PACE', 'UNKNOWN', ...
             'NOISE', 'ARFCT', 'STCH', 'TCH', 'SYSTOLE', 'DIASTOLE', ...
             'NOTE', 'MEASURE', 'PWAVE', 'BBB', 'PACESP', 'TWAVE', ...
             'RHYTHM', 'UWAVE', 'LEARN', 'FLWAV', 'VFON', 'VFOFF', ...
             'AESC', 'SVESC', 'LINK', 'NAPC', 'PFUS', 'WFON', 'WFOFF', 'RONT'};
labelDict = containers.Map(annoCodes, annoNames);

figure;
set(gcf, 'Position', [80, 80, 1200, 600]);
for ch = 1:2
    subplot(2, 1, ch);
    plot(timeAxis, Sig(:, ch));
    title(['channel ' num2str(ch)]);
    xlabel('time(s)');
    ylabel('amp');
    grid minor;
    hold on;
    
    for k = 1:length(ATRTIMED)
        t_val = ATRTIMED(k);
        code = ANNOTD(k);

        if isKey(labelDict, code)
            strLabel = labelDict(code);
        else
            strLabel = ['Anomaly ' num2str(code)];
        end

        if code == 1
            txtColor = 'green';
        else
            txtColor = 'red';
        end

        yPos = Sig(floor(t_val * sfreq), 1) - 0.2;

        text(t_val, yPos, strLabel, 'Color', txtColor, ...
             'FontSize', 6, 'HorizontalAlignment', 'center');
    end
    
    hold off;
    axis tight;
end
