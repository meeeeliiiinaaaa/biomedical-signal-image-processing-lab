%% Part3

%%

load("X.dat");

%%

sig = X;
sig_transposed = sig';
Fs = 256;

%% Section1

[W, Zhat]= ica(sig_transposed);
A = pinv(W);

%% Section2

plot_title = 'Raw ECG signal';
plot3ch(sig ,Fs , plot_title);

%%

colors={'r', 'g', 'b'};
s = 100;

for i=1:3
    plot3dv(A(:, i), s, colors{i});
end

xlabel('X');
ylabel('Y');
zlabel('Z');
title('3D Plot of Columns of Matrix A');
grid on;

%% Section3

for i=1:3
     subplot(3,1,i)
     plot(Zhat(i,:));
end

sgtitle('3 rows of matrix Zhat', ...
        'FontSize', 11, 'FontWeight', 'bold');

%%

 A_modified = A;
 A_modified(:, 1) = 0;  % Larger amplitude ---> Maternal ECG
 A_modified(:, 2) = 0;  % Seems noise

sig_separated_ica = A_modified * Zhat;

%% Section4

plot_title = "Separated ECG signal";
plot3ch(sig_separated_ica' ,Fs , plot_title);








%% ica function

function [W, Zhat] = ica(X)
%ICA  Perform independent component analysis
%  [W, ZHAT] = ICA(X) performs independent component analysis on data
%  observation matrix X. Matrix X is a transposed observation matrix, such that
%  each row of X represents an observed signal. This approach uses Cardoso's ICA
%  algorithm to estimate sources (ZHAT) and the de-mixing matrix W, an
%  approximation to A^{-1}, the original (unknown) mixing matrix. 

% Created by: G D Clifford 2004  gari AT mit DOT edu
% Last Modified: 5/7/06, documentation updated.

% Input argument checking
%------------------------
[a, b] = size(X);
if a > b
    fprintf('Warning - ICA works across the rows of the input data.\n');
    error('Please transpose input.');
end
Nsources = a;

if Nsources > min([a b])
    Nsources = min([a b]);
    fprintf('Warning - number of soures cannot exceed number of observation channels \n')
    fprintf(' ... reducing to %i \n',Nsources)
end

%tic
[Winv, Zhat] = jade(X,Nsources);
W = pinv(Winv);
%fprintf('algorithm timing ...  ')
%toc
end


%% plot3ch function

function plot3ch(X,Fs,plot_title)
%PLOT3CH  Plot 3 channel data in the time-domain and on a 3D scatter plot
%  PLOT3CH(X,FS,'TITLE') plots the three columns of data matrix X on a
%  time-domain plot with sample rate FS on and plots each column against the
%  other on a 3D scatter plot. The default value for FS is 256 Hz. The optional
%  'TITLE' input allows the user to specify the plot title string.
% Created by: G.D. Clifford 2004 gari AT mit DOT edu
% Modified 5/6/05, Eric Weiss. Documentation updated. Plot title input added.
% Input argument checking
%------------------------
if nargin < 2
    Fs = 256;
end;
if nargin < 3
    plot_title = '3 Channel Data';
end;
[M,N] = size(X);
if N ~= 3;
    error('Input matrix must have 3 columns!');
end;
% Generate time-domain plot
%--------------------------
t = [1/Fs:1/Fs:M/Fs];
figure;
for i = 1:N
    subplot(N,1,i)
    plot(t,X(:,i)); ylabel(['Ch',int2str(i)]);
    axis([0 max(t) min(X(:,i))-abs(0.1*max(X(:,i))) max(X(:,i))+abs(0.1*max(X(:,i)))]);
    %axis([0 max(t) min(X(:,i)) max(X(:,i))])
end;
xlabel('Time (sec)');
subplot(N,1,1); title(plot_title);
figure;
plot3(X(:,1), X(:,2), X(:,3),'.m');
xlabel('Ch1'); ylabel('Ch2'); zlabel('Ch3');
title(plot_title);
grid on;
end
%% plot3dv function

function plot3dv(v, s, col)
%PLOT3DV  Plots the specified vector onto a 3D scatter plot
%  PLOT3DV(V, S, 'COL') plots the eigenvector +/-V with singular value S and
%  color 'COL' onto a 3D plot of the currently displayed figure. The length of
%  the plotted eigenvector is equal to the square root of the singular value. If
%  the singular value S is not specified, the default scaling length is 10. If
%  the color 'COL' is not specified, the default color is 'k' (black).

% Created by: GD Clifford 2004 gari AT alum DOT mit DOT edu
% Last modified 5/7/06, Eric Weiss. Documentation updated.

% Input argument checking
%------------------------
if nargin < 2 | isempty(s)
    s = 100;
end;
if nargin < 3
    col = 'k';
end;
v = v(:); % ensure that eigenvector is in column format
[m, n] = size(v);
if (n ~= 1 | m ~= 3)
    error('vector must be 3x1')
end;
if s == 1  % legacy code: does not affect function
    ln = 1/sqrt((v(1)*v(1))+(v(2)*v(2))+(v(3)*v(3)));
end;

% Plot eigenvector on 3D plot
%----------------------------
sn = sqrt(s);
hold on;
plot3(sn*[-1*v(1) v(1)],sn*[-1*v(2) v(2)],sn*[-1*v(3) v(3)],col);
grid on;
view([1,1,1])

end

%% jade function
function [A,S]=jade(X,m)

% Source separation of complex signals with JADE.
% Jade performs `Source Separation' in the following sense:
%   X is an n x T data matrix assumed modelled as X = A S + N where
% 
% o A is an unknown n x m matrix with full rank.
% o S is a m x T data matrix (source signals) with the properties
%    	a) for each t, the components of S(:,t) are statistically
%    	   independent
% 	b) for each p, the S(p,:) is the realization of a zero-mean
% 	   `source signal'.
% 	c) At most one of these processes has a vanishing 4th-order
% 	   cumulant.
% o  N is a n x T matrix. It is a realization of a spatially white
%    Gaussian noise, i.e. Cov(X) = sigma*eye(n) with unknown variance
%    sigma.  This is probably better than no modeling at all...
%
% Jade performs source separation via a 
% Joint Approximate Diagonalization of Eigen-matrices.  
%
% THIS VERSION ASSUMES ZERO-MEAN SIGNALS
%
% Input :
%   * X: Each column of X is a sample from the n sensors
%   * m: m is an optional argument for the number of sources.
%     If ommited, JADE assumes as many sources as sensors.
%
% Output :
%    * A is an n x m estimate of the mixing matrix
%    * S is an m x T naive (ie pinv(A)*X)  estimate of the source signals
%
%
% Version 1.5.  Copyright: JF Cardoso.  
%
% See notes, references and revision history at the bottom of this file



[n,T]	= size(X);

%%  source detection not implemented yet !
if nargin==1, m=n ; end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A few parameters that could be adjusted
nem	= m;		% number of eigen-matrices to be diagonalized
seuil	= 1/sqrt(T)/100;% a statistical threshold for stopping joint diag


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% whitening
%
if m<n, %assumes white noise
 	[U,D] 	= eig((X*X')/T); 
	[puiss,k]=sort(diag(D));
 	ibl 	= sqrt(puiss(n-m+1:n)-mean(puiss(1:n-m)));
 	bl 	= ones(m,1) ./ ibl ;
 	W	= diag(bl)*U(1:n,k(n-m+1:n))';
 	IW 	= U(1:n,k(n-m+1:n))*diag(ibl);
else    %assumes no noise
 	IW 	= sqrtm((X*X')/T);
 	W	= inv(IW);
end;
Y	= W*X;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Cumulant estimation


R	= (Y*Y' )/T ;
C	= (Y*Y.')/T ;

Yl	= zeros(1,T);
Ykl	= zeros(1,T);
Yjkl	= zeros(1,T);

Q	= zeros(m*m*m*m,1) ;
index	= 1;

for lx = 1:m ; Yl 	= Y(lx,:);
for kx = 1:m ; Ykl 	= Yl.*conj(Y(kx,:));
for jx = 1:m ; Yjkl	= Ykl.*conj(Y(jx,:));
for ix = 1:m ; 
	Q(index) = ...
	(Yjkl * Y(ix,:).')/T -  R(ix,jx)*R(lx,kx) -  R(ix,kx)*R(lx,jx) -  C(ix,lx)*conj(C(jx,kx))  ;
	index	= index + 1 ;
end ;
end ;
end ;
end

%% If you prefer to use more memory and less CPU, you may prefer this
%% code (due to J. Galy of ENSICA) for the estimation the cumulants
%ones_m = ones(m,1) ; 
%T1 	= kron(ones_m,Y); 
%T2 	= kron(Y,ones_m);  
%TT 	= (T1.* conj(T2)) ;
%TS 	= (T1 * T2.')/T ;
%R 	= (Y*Y')/T  ;
%Q	= (TT*TT')/T - kron(R,ones(m)).*kron(ones(m),conj(R)) - R(:)*R(:)' - TS.*TS' ;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%computation and reshaping of the significant eigen matrices

[U,D]	= eig(reshape(Q,m*m,m*m)); 
[la,K]	= sort(abs(diag(D)));

%% reshaping the most (there are `nem' of them) significant eigenmatrice
M	= zeros(m,nem*m);	% array to hold the significant eigen-matrices
Z	= zeros(m)	; % buffer
h	= m*m;
for u=1:m:nem*m, 
	Z(:) 		= U(:,K(h));
	M(:,u:u+m-1)	= la(h)*Z;
	h		= h-1; 
end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% joint approximate diagonalization of the eigen-matrices


%% Better declare the variables used in the loop :
B 	= [ 1 0 0 ; 0 1 1 ; 0 -i i ] ;
Bt	= B' ;
Ip	= zeros(1,nem) ;
Iq	= zeros(1,nem) ;
g	= zeros(3,nem) ;
G	= zeros(2,2) ;
vcp	= zeros(3,3);
D	= zeros(3,3);
la	= zeros(3,1);
K	= zeros(3,3);
angles	= zeros(3,1);
pair	= zeros(1,2);
c	= 0 ;
s	= 0 ;


%init;
encore	= 1;
V	= eye(m); 

% Main loop
while encore, encore=0;
 for p=1:m-1,
  for q=p+1:m,

 	Ip = p:m:nem*m ;
	Iq = q:m:nem*m ;

	% Computing the Givens angles
 	g	= [ M(p,Ip)-M(q,Iq)  ; M(p,Iq) ; M(q,Ip) ] ; 
 	[vcp,D] = eig(real(B*(g*g')*Bt));
	[la, K]	= sort(diag(D));
 	angles	= vcp(:,K(3));
	if angles(1)<0 , angles= -angles ; end ;
 	c	= sqrt(0.5+angles(1)/2);
 	s	= 0.5*(angles(2)-j*angles(3))/c; 

 	if abs(s)>seuil, %%% updates matrices M and V by a Givens rotation
	 	encore 		= 1 ;
		pair 		= [p;q] ;
 		G 		= [ c -conj(s) ; s c ] ;
		V(:,pair) 	= V(:,pair)*G ;
	 	M(pair,:)	= G' * M(pair,:) ;
		M(:,[Ip Iq]) 	= [ c*M(:,Ip)+s*M(:,Iq) -conj(s)*M(:,Ip)+c*M(:,Iq) ] ;
 	end%% if
  end%% q loop
 end%% p loop
end%% while

%%%estimation of the mixing matrix and signal separation
A	= IW*V;
S	= V'*Y ;

return ;

end
