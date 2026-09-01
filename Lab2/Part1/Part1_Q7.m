clearvars; clc; close all;

origStruct = load('X_org.mat');
noiseStruct = load('X_noise.mat');
cleanData = origStruct.X_org;
noiseData = noiseStruct.X_noise;

pwrClean = norm(cleanData, 'fro')^2;
pwrNoise = norm(noiseData, 'fro')^2;

ratio = pwrClean / pwrNoise;
scaleA = sqrt(ratio * sqrt(10));
scaleB = sqrt(ratio * (10 * sqrt(10)));

corruptedSig1 = cleanData + (scaleA * noiseData);
corruptedSig2 = cleanData + (scaleB * noiseData);

[MixA, ~, ~] = COM2R(corruptedSig1, 32);
[MixB, ~, ~] = COM2R(corruptedSig2, 32);


targetComps1 = [2, 5, 9, 10, 11, 12, 21, 23];
discardIdx1 = setdiff(1:32, targetComps1);
MixA(:, discardIdx1) = 0;

targetComps2 = [7, 8, 9, 15, 18, 19, 20, 22, 25, 30];
discardIdx2 = setdiff(1:32, targetComps2);
MixB(:, discardIdx2) = 0;

denoisedData1 = MixA * corruptedSig1;
denoisedData2 = MixB * corruptedSig2;


elecData = load('Electrodes.mat');
targetLabels = [elecData.Electrodes.labels(13); elecData.Electrodes.labels(24)];

fig1 = figure('Name', 'Analysis: Scenario A', 'Color', 'w', 'Position', [100 100 800 700]);
subplot(3, 1, 1);
renderEEGTraces([corruptedSig1(13,:); corruptedSig1(24,:)], targetLabels, 'Corrupted Data (Noise Level 1)');
subplot(3, 1, 2);
renderEEGTraces([cleanData(13,:); cleanData(24,:)], targetLabels, 'Reference Signal (Clean)');
subplot(3, 1, 3);
renderEEGTraces([denoisedData1(13,:); denoisedData1(24,:)], targetLabels, 'Reconstructed Data (ICA Filtered)');

fig2 = figure('Name', 'Analysis: Scenario B', 'Color', 'w', 'Position', [150 150 800 700]);
subplot(3, 1, 1);
renderEEGTraces([corruptedSig2(13,:); corruptedSig2(24,:)], targetLabels, 'Corrupted Data (Noise Level 2)');
subplot(3, 1, 2);
renderEEGTraces([cleanData(13,:); cleanData(24,:)], targetLabels, 'Reference Signal (Clean)');
subplot(3, 1, 3);
renderEEGTraces([denoisedData2(13,:); denoisedData2(24,:)], targetLabels, 'Reconstructed Data (ICA Filtered)');



function renderEEGTraces(dataMat, labels, plotTitle)
    [numCh, numSamples] = size(dataMat);
    fs = 250;
    tVec = (1:numSamples) / fs;
    gap = max(abs(dataMat(:)));
    offsets = gap * (0:-1:-(numCh-1))';
    shiftedData = dataMat + offsets;

    ax = gca;
    plot(ax, tVec, shiftedData', 'LineWidth', 1.5);

    set(ax, 'YTick', flipud(offsets), 'YTickLabel', flipud(labels), 'FontSize', 9);
    ylim([min(shiftedData(:)) - gap*0.3, max(shiftedData(:)) + gap*0.3]);
    xlim([tVec(1) tVec(end)]);
    
    xlabel('Time [s]', 'FontWeight', 'bold');
    ylabel('Sensors', 'FontWeight', 'bold');
    title(plotTitle, 'FontWeight', 'bold', 'Color', [0.1, 0.3, 0.5], 'FontSize', 11);
    
    grid on;
    ax.GridLineStyle = '--';
    ax.GridAlpha = 0.5;
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