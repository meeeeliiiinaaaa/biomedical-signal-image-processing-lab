%% Part2 - Data2

%% 

load("NewData3.mat");
load("Electrodes.mat");
%%

sig = EEG_Sig;
fs = 250;
channel_labels = Electrodes.labels;
X_locs = Electrodes.X;
Y_locs = Electrodes.Y;

%% Q1

plotEEG(sig);

%% Q3

N = size(sig, 1);
[A, W, k] = COM2R(sig, N);
S = W * sig; 

%% Q4

plot_ica(S, A, fs, X_locs, Y_locs, channel_labels);

%% Q5

SelSources = [2, 3, 6, 8, 9, 11, 12, 14, 15, 17, 18, 19, 20, 21]; 

sig_denoised = A(:,SelSources) * S(SelSources,:);

%% Q6

plotEEG(sig_denoised);


%% plot EEG function

function plotEEG(X) 

load("Electrodes.mat") ;
% Plot Data
% Use function disp_eeg(X,offset,feq,ElecName)
offset = max(abs(X(:))) ;
feq = 250 ;
ElecName = Electrodes.labels ;
disp_eeg(X,offset,feq,ElecName);

end
%% disp_eeg function

function t = disp_eeg(X,offset,feq,ElecName,titre)
% function t = disp_eeg(X,offset,feq,ElecName,titre)
%
% inputs
%     X: dynamics to display. (nbchannels x nbsamples) matrix
%     offset: offset between channels (default max(abs(X)))
%     feq: sapling frequency (default 1)
%     ElecName: cell array of electrode labels (default {S1,S2,...})
%     titre: title of the figure
%
% output
%     t: time vector
%
% G. Birot 2010-02


%% Check arguments
[N K] = size(X);

if nargin < 4
    for n = 1:N
        ElecName{n}  = ['S',num2str(n)];
    end
    titre = [];
end

if nargin < 5
    titre = [];
end

if isempty(feq)
    feq = 1;
end

if isempty(ElecName)
    for n = 1:N
        ElecName{n}  = ['S',num2str(n)];
    end
end

if isempty(offset)
    offset = max(abs(X(:)));
end


%% Build dynamic matrix with offset and time vector
X = X + repmat(offset*(0:-1:-(N-1))',1,K);
t = (1:K)/feq;
graduations = offset*(0:-1:-(N-1))';
shiftvec = N:-1:1;
Ysup = max(X(1,:)) + offset;
Yinf = min(X(end,:)) - offset;
% YLabels = cell(N+2) ElecName(shiftvec)

%% Display
figure1 = figure;
% a1 = axes('YAxisLocation','right');
a2 = axes('YTickLabel',ElecName(shiftvec),'YTick',graduations(shiftvec),'FontSize',7);
ylim([Yinf Ysup]);
box('on');
grid('on')
hold('all');
plot(t,X');
xlabel('Time (seconds)','FontSize',10);
ylabel('Channels','FontSize',10);
title(titre);
hold off

end

%% plot_ica function

function plot_ica(S, A, Fs, elocsX, elocsY, elabels)
    
    num_comps = size(S, 1);
    N = length(S(1, :));
    disp(N);
    time = (0:N-1) / Fs; 

    figure('Name', 'All ICA Components', 'NumberTitle', 'off', 'Units', 'normalized', 'Position', [0 0 1 1]);
    
    for i = 1:num_comps
        row = ceil(i / 3);
        col_offset = mod(i-1, 3) * 3;

        subplot(7, 9, (row-1)*9 + col_offset + 1);
        plot(time, S(i, :), 'b-', 'LineWidth', 1);
        title(['Time IC ', num2str(i)]);
        grid("on");
        subplot(7, 9, (row-1)*9 + col_offset + 2);
        [pxx, f] = pwelch(S(i, :), [], [], [], Fs);
        plot(f, 10*log10(pxx), 'b-', 'LineWidth', 1);
        title(['PSD IC ', num2str(i)]);
        xlim([0 100]);
        grid("on");

        subplot(7, 9, (row-1)*9 + col_offset + 3);
        plottopomap(elocsX, elocsY, elabels, A(:, i)); 
        title(['Map IC ', num2str(i)]);
    end
end

%% tfuni function

function [S,A]=tfuni4(e)
% [S,A]=tfuni4(e)
% Comon, version 12 feb 1992
% English comments added in 1994
% Orthogonal real direct transform
% for separating 2 sources in presence of noise
% Sources are assumed zero mean
%
T=length(e);
%%%%% moments d'ordre 2
 g11=e(1,:)*e(1,:)'/T;%cv vers 1
 g22=e(2,:)*e(2,:)'/T;%cv vers 1
 g12=e(1,:)*e(2,:)'/T;%cv vers 0
%%%%% moments d'ordre 4
e2=e.^2;
 g1111=e2(1,:)*e2(1,:)'/T;
 g1112=e2(1,:).*e(1,:)*e(2,:)'/T;
 g1122=e2(1,:)*e2(2,:)'/T;
 g1222=e2(2,:).*e(2,:)*e(1,:)'/T;
 g2222=e2(2,:)*e2(2,:)'/T;
%%%%% cumulants croises d'ordre 4
 q1111=g1111-3*g11*g11;
 q1112=g1112-3*g11*g12;
 q1122=g1122-g11*g22-2*g12*g12;
 q1222=g1222-3*g22*g12;
	q2222=g2222-3*g22*g22;
%%%%% racine de Pw(x): si t est la tangente de l'angle, x=t-1/t.
u=q1111+q2222-6*q1122;v=q1222-q1112;z=q1111*q1111+q2222*q2222;
c4=q1111*q1112-q2222*q1222;
c3=z-4*(q1112*q1112+q1222*q1222)-3*q1122*(q1111+q2222);
c2=3*v*u;
c1=3*z-2*q1111*q2222-32*q1112*q1222-36*q1122*q1122;
c0=-4*(u*v+4*c4);
%c0=q1112*q2222-q1222*q1111-3*q1112*q1111+3*q1222*q2222-6*q1122*q1112+6*q1122*q1222;c0=4*c0
Pw=[c4 c3 c2 c1 c0];R=roots(Pw);I=find(abs(imag(R))<eps);xx=R(I);
%%%%% maximum du contraste en x
a0=q1111;a1=4*q1112;a2=6*q1122;a3=4*q1222;a4=q2222;
b4=a0*a0+a4*a4;
b3=2*(a3*a4-a0*a1);
b2=4*a0*a0+4*a4*a4+a1*a1+a3*a3+2*a0*a2+2*a2*a4;
b1=2*(-3*a0*a1+3*a3*a4+a1*a4+a2*a3-a0*a3-a1*a2);
b0=2*(a0*a0+a1*a1+a2*a2+a3*a3+a4*a4+2*a0*a2+2*a0*a4+2*a1*a3+2*a2*a4);
Pk=[b4 b3 b2 b1 b0];  % numerateur du contraste
Wk=polyval(Pk,xx);Vk=polyval([1 0 8 0 16],xx);Wk=Wk./Vk;
[Wmax,j]=max(Wk);Xsol=xx(j);
%%%%% maximum du contratse en theta
t=roots([1 -Xsol -1]);j=find(t<=1 & t>-1);t=t(j);
%%%%% test et conditionnement
if abs(t)<1/T,
  A=eye(2); %fprintf('pas de rotation plane pour cette paire\n');
else,
  A(1,1)=1/sqrt(1+t*t);A(2,2)=A(1,1);A(1,2)=t*A(1,1);A(2,1)=-A(1,2);
end;
%%%%% filtrage de la sortie
 S=A*e;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

%% Com2 function

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


%%  plottopomap function

function plottopomap(elocsX,elocsY,elabels,data)

% define XY points for interpolation
interp_detail = 100;
interpX = linspace(min(elocsX)-.2,max(elocsX)+.25,interp_detail);
interpY = linspace(min(elocsY),max(elocsY),interp_detail);

% meshgrid is a function that creates 2D grid locations based on 1D inputs
[gridX,gridY] = meshgrid(interpX,interpY);
% Interpolate the data on a 2D grid
interpFunction = TriScatteredInterp(elocsY,elocsX,data);
topodata = interpFunction(gridX,gridY);

% plot map
contourf(interpY,interpX,topodata);
hold on
scatter(elocsY,elocsX,10,'ro','filled');
for i=1:length(elocsX)
    text(elocsY(i),elocsX(i),elabels(i))
end
set(gca,'xtick',[])
set(gca,'ytick',[])

end

