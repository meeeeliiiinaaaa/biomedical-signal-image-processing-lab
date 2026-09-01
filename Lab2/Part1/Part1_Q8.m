clearvars; clc; close all;

refData = load('X_org.mat');
errData = load('X_noise.mat');
sigClean = refData.X_org;
sigNoise = errData.X_noise;

powClean = norm(sigClean, 'fro')^2;
powNoise = norm(sigNoise, 'fro')^2;
powerRatio = powClean / powNoise;

coeff1 = sqrt(powerRatio * 10^0.5);
coeff2 = sqrt(powerRatio * 10^1.5);

mixedSignalA = sigClean + (coeff1 * sigNoise);
mixedSignalB = sigClean + (coeff2 * sigNoise);

[MixMatA, ~, ~] = COM2R(mixedSignalA, 32);
[MixMatB, ~, ~] = COM2R(mixedSignalB, 32);

keepIdxA = [16, 22, 23, 26];
removeIdxA = setdiff(1:32, keepIdxA);
MixMatA(:, removeIdxA) = 0;

keepIdxB = [3, 13, 16, 19, 24];
removeIdxB = setdiff(1:32, keepIdxB);
MixMatB(:, removeIdxB) = 0;

reconSignalA = MixMatA * mixedSignalA;
reconSignalB = MixMatB * mixedSignalB;

errorMetric1 = norm(sigClean - reconSignalA, 'fro') / norm(sigClean, 'fro');
errorMetric2 = norm(sigClean - reconSignalB, 'fro') / norm(sigClean, 'fro');

fprintf('\n==================================================\n');
fprintf('       PERFORMANCE EVALUATION (RRMSE RESULTS)       \n');
fprintf('==================================================\n');
fprintf('  [*] Scenario 1 (Lower Noise Level)  : %0.6f\n', errorMetric1);
fprintf('  [*] Scenario 2 (Higher Noise Level) : %0.6f\n', errorMetric2);
fprintf('==================================================\n\n');



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
