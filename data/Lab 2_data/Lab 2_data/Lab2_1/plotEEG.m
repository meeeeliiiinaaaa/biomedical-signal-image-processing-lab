function plotEEG(X) 
load('Electrodes') ;

% Plot Data
% Use function disp_eeg(X,offset,feq,ElecName)
offset = max(abs(X(:))) ;
feq = 250 ;
ElecName = Electrodes.labels ;
disp_eeg(X,offset,feq,ElecName);