clearvars; clc; close all;

data_org = load('X_org.mat');
data_noise = load('X_noise.mat');
sigOriginal = data_org.X_org;
sigNoise = data_noise.X_noise;

pwrClean = sum(sigOriginal(:).^2);
pwrNoise = sum(sigNoise(:).^2);

alpha1 = sqrt((pwrClean / pwrNoise) * sqrt(10)); % 10^0.5
alpha2 = sqrt((pwrClean / pwrNoise) * (10 * sqrt(10))); % 10^1.5

mixedData_1 = sigOriginal + (alpha1 * sigNoise);
mixedData_2 = sigOriginal + (alpha2 * sigNoise);

[MixMat1, UnmixMat1, ~] = COM2R(mixedData_1, 32);
[MixMat2, UnmixMat2, ~] = COM2R(mixedData_2, 32);


src1 = UnmixMat1 * mixedData_1;
src2 = UnmixMat2 * mixedData_2;

render_eeg_components(src1, 'Extracted Sources - Scenario A');
render_eeg_components(src2, 'Extracted Sources - Scenario B');

keepIdx1 = [2, 5, 9, 10, 11, 12, 21, 23];
removeIdx1 = setdiff(1:32, keepIdx1);
MixMat1(:, removeIdx1) = 0;

keepIdx2 = [7, 8, 9, 15, 18, 19, 20, 22, 25, 30];
removeIdx2 = setdiff(1:32, keepIdx2);
MixMat2(:, removeIdx2) = 0;



function render_eeg_components(signalData, plotTitle)
   
    try
        load('Electrodes.mat', 'Electrodes');
        chLabels = Electrodes.labels;
    catch
        [numCh, ~] = size(signalData);
        chLabels = arrayfun(@(x) sprintf('Ch-%d', x), 1:numCh, 'UniformOutput', false);
    end
    
    fs = 250;
    gap = max(abs(signalData(:))) / 1.8;
    
    if nargin < 2
        plotTitle = 'EEG Source Signals';
    end
    
    [nChannels, nSamples] = size(signalData);
    timeSteps = (1:nSamples) / fs;
    
    offsets = gap * (0:-1:-(nChannels-1))';
    plotData = signalData + (offsets * ones(1, nSamples));

    f = figure('Color', [0.95 0.95 0.95], 'Name', 'Source Viewer', 'NumberTitle', 'off');
    ax = axes('Parent', f);
    
    plot(ax, timeSteps, plotData', 'LineWidth', 0.9, 'Color', [0.1 0.3 0.6]);
    
    yLimits = [min(plotData(end,:)) - gap, max(plotData(1,:)) + gap];
    ylim(ax, yLimits);
    xlim(ax, [timeSteps(1) timeSteps(end)]);

    revIdx = nChannels:-1:1;
    set(ax, 'YTick', offsets(revIdx), 'YTickLabel', chLabels(revIdx), ...
            'FontSize', 8, 'FontName', 'Arial', 'TickDir', 'out');
    
    grid(ax, 'on');
    set(ax, 'GridColor', [0.7 0.7 0.7], 'GridAlpha', 0.6, 'MinorGridLineStyle', '-');
    
    xlabel(ax, 'Time [sec]', 'FontSize', 10, 'FontWeight', 'bold');
    ylabel(ax, 'Components', 'FontSize', 10, 'FontWeight', 'bold');
    title(ax, plotTitle, 'FontSize', 12, 'FontWeight', 'bold');
end

function [F,W,K]=COM2R(Y,Pest)
    disp('COM2')
    % Comon, version 6 march 92
    % English comments added in 1994
    % [F,delta]=aci(Y)
    % Y is the observations matrix
    % This routine outputs a matrix F such that Y=F*Z, Z=pinv(F)*Y,
    % and components of Z are uncorrelated and approximately independent
    % F is Nxr, where r is the dimension Z;
    % Entries of delta are sorted in decreasing order;
    % Columns of F are of unit norm;
    % The entry of largest modulus in each column of F is positive real.
    % Initial and final values of contrast can be fprinted for checking.
    % REFERENCE: P.Comon, "Independent Component Analysis, a new concept?",
    % Signal Processing, Elsevier, vol.36, no 3, April 1994, 287-314.
    %
    [N,TT]=size(Y);T=max(N,TT);N=min(N,TT);
    if TT==N, Y=Y';[N,T]=size(Y);end; % Y est maintenant NxT avec N<T.
    %%%% STEPS 1 & 2: whitening and projection (PCA)
    [U,S,V]=svd(Y',0);tol=max(size(S))*norm(S)*eps;
    s=diag(S);I=find(s<tol);
    
    %--- modif de Laurent le 03/02/2009
    r = min(Pest,N);
    U=U(:,1:r);
    S=S(1:r,1:r);
    V=V(:,1:r);
    %---
    
    Z=U'*sqrt(T);L=V*S'/sqrt(T);F=L; %%%%%% on a Y=L*Z;
    %%%%%% INITIAL CONTRAST
    T=length(Z);contraste=0;
    for i=1:r,
     gii=Z(i,:)*Z(i,:)'/T;Z2i=Z(i,:).^2;;giiii=Z2i*Z2i'/T;
     qiiii=giiii/gii/gii-3;contraste=contraste+qiiii*qiiii;
    end;
    %%%% STEPS 3 & 4 & 5: Unitary transform
    S=Z;
    if N==2,K=1;else,K=1+round(sqrt(N));end;  % K= max number of sweeps
    Rot=eye(r);
    for k=1:K,                           %%%%%% strating sweeps
    Q=eye(r);
      for i=1:r-1,
      for j= i+1:r,
        S1ij=[S(i,:);S(j,:)];
        [Sij,qij]=tfuni4(S1ij);    %%%%%% processing a pair
        S(i,:)=Sij(1,:);S(j,:)=Sij(2,:);
        Qij=eye(r);Qij(i,i)=qij(1,1);Qij(i,j)=qij(1,2);
        Qij(j,i)=qij(2,1);Qij(j,j)=qij(2,2);
        Q=Qij*Q;
      end;
      end;
    Rot=Rot*Q';
    end;                                    %%%%%% end sweeps
    F=F*Rot;
    %%%%%% FINAL CONTRAST
    S=Rot'*Z;
    T=length(S);contraste=0;
    for i=1:r,
     gii=S(i,:)*S(i,:)'/T;S2i=S(i,:).^2;;giiii=S2i*S2i'/T;
     qiiii=giiii/gii/gii-3;contraste=contraste+qiiii*qiiii;
    end;
    %%%% STEP 6: Norming columns
    delta=diag(sqrt(sum(F.*conj(F))));
    %%%% STEP 7: Sorting
    [d,I]=sort(-diag(delta));E=eye(r);P=E(:,I)';delta=P*delta*P';F=F*P';
    %%%% STEP 8: Norming
    F=F*inv(delta);
    %%%% STEP 9: Phase of columns
    [y,I]=max(abs(F));
    for i=1:r,Lambda(i)=conj(F(I(i),i));end;Lambda=Lambda./abs(Lambda);
    F=F*diag(Lambda);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % CALCUL DE LA MATRICE DE FILTRAGE
    %---------------------------------
    W = pinv(F);
end
