
clearvars; clc; close all;

dataOrig = load('X_org.mat');
dataNoise = load('X_noise.mat');
cleanSig = dataOrig.X_org;
noiseSig = dataNoise.X_noise;

pwrClean = sum(cleanSig(:).^2);
pwrNoise = sum(noiseSig(:).^2);

scale1 = sqrt((pwrClean / pwrNoise) * sqrt(10));     
scale2 = sqrt((pwrClean / pwrNoise) * (10 * sqrt(10))); 

mixedSig1 = cleanSig + (scale1 * noiseSig);
mixedSig2 = cleanSig + (scale2 * noiseSig);

[MixMat1, UnmixMat1, ~] = COM2R(mixedSig1, 32);
[MixMat2, UnmixMat2, ~] = COM2R(mixedSig2, 32);

src_est_1 = UnmixMat1 * mixedSig1;
src_est_2 = UnmixMat2 * mixedSig2;

validComps1 = [2, 5, 9, 10, 11, 12, 21, 23];
removeIdx1 = setdiff(1:32, validComps1);
MixMat1(:, removeIdx1) = 0;

validComps2 = [7, 8, 9, 15, 18, 19, 20, 22, 25, 30];
removeIdx2 = setdiff(1:32, validComps2);
MixMat2(:, removeIdx2) = 0;


denoisedSig1 = MixMat1 * mixedSig1;
denoisedSig2 = MixMat2 * mixedSig2;


fig = figure('Color', 'w', 'Name', 'Denoised Signals Viewer', 'Position', [100, 100, 800, 600]);

ax1 = subplot(2, 1, 1);
plot(ax1, denoisedSig1(1:min(5, end), :)', 'LineWidth', 1.2);
title(ax1, 'Reconstructed Data - Configuration A', 'FontWeight', 'bold', 'Color', [0.15 0.15 0.55]);
xlabel(ax1, 'Time Samples', 'FontWeight', 'bold');
ylabel(ax1, 'Amplitude', 'FontWeight', 'bold');
grid(ax1, 'on'); set(ax1, 'GridLineStyle', ':', 'GridAlpha', 0.6);


ax2 = subplot(2, 1, 2);
plot(ax2, denoisedSig2(1:min(5, end), :)', 'LineWidth', 1.2);
title(ax2, 'Reconstructed Data - Configuration B', 'FontWeight', 'bold', 'Color', [0.15 0.15 0.55]);
xlabel(ax2, 'Time Samples', 'FontWeight', 'bold');
ylabel(ax2, 'Amplitude', 'FontWeight', 'bold');
grid(ax2, 'on'); set(ax2, 'GridLineStyle', ':', 'GridAlpha', 0.6);


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
