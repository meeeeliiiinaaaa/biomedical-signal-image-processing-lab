clear; clc; close all;
load('ECG_sig.mat');

getSamp = @(idx) floor(ATRTIMED(idx) * sfreq);
getAmp  = @(sig) max(sig) - min(sig);

idxN = find(ANNOTD == 1);
idxA = find(ANNOTD ~= 1);

featNorm = struct('amp', [], 'qrs', [], 'slope', []);
featAnom = struct('amp', [], 'qrs', [], 'slope', []);

for k = 1:(length(idxN) - 1)
    segN = Sig(getSamp(idxN(k)-1) : getSamp(idxN(k)+1), 1);
    
    featNorm.amp(k, 1)   = getAmp(segN);
    featNorm.qrs(k, 1)   = extract_qrs_width(segN, sfreq);
    featNorm.slope(k, 1) = max(diff(segN));
end

for k = 2:length(idxA)
    segA = Sig(getSamp(idxA(k)-1) : getSamp(idxA(k)+1)+1, :); 
    idx_store = k - 1;
    featAnom.amp(idx_store, 1)   = getAmp(segA(:, 1));
    featAnom.qrs(idx_store, 1)   = extract_qrs_width(segA(:, 1), sfreq);
    featAnom.slope(idx_store, 1) = max(diff(segA(:, 1)));
end
print_metrics('Normal Beats Features:', featNorm);
print_metrics('Anomalous Beats Features:', featAnom);

figure;
t_norm = getSamp(idxN(1)-1) : getSamp(idxN(1)+1);
plot(t_norm / sfreq, Sig(t_norm, 1));
title('Normal'); xlabel('Time (seconds)'); ylabel('Amplitude');
grid minor; axis tight;

[~, unique_idx] = unique(ANNOTD(idxA(2:end)), 'stable');

for c = 1:length(unique_idx)
    targetIdx = idxA(unique_idx(c) + 1);
    disp(ANNOTD(targetIdx))
    t_anom = getSamp(targetIdx-1) : getSamp(targetIdx+1);
    
    figure;
    plot(t_anom / sfreq, Sig(t_anom, 1));
    title({'Anomally', num2str(c)});
    xlabel('time(s)'); ylabel('Amp');
    grid minor; axis tight;
end
function w = extract_qrs_width(signal_beat, fs)
    threshold = 0.5 * max(signal_beat);
    crossings = find(signal_beat > threshold);
    if isempty(crossings)
        w = 0; 
    else
        w = (crossings(end) - crossings(1)) / fs;
    end
end

function print_metrics(headerText, fStruct)
    disp(headerText);
    disp('Amplitude (Mean ± Std):');
    disp([mean(fStruct.amp), std(fStruct.amp)]);
    disp('QRS Duration (Mean ± Std):');
    disp([mean(fStruct.qrs), std(fStruct.qrs)]);
    disp('R Peak Slope (Mean ± Std):');
    disp([mean(fStruct.slope), std(fStruct.slope)]);
end
