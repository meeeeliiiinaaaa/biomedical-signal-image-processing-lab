%% Part2

%% Section1

clc; close all;

load('X.dat');
plot3ch(X);
[U, S, V] = svd(X);

%% Section2

clc; close all;
colors={'r', 'g', 'b'};

for i = 1:3
    plot3dv(V(:,i), S(i,i), colors{i});
end


xlabel('X');
ylabel('Y');
zlabel('Z');
title('3D Plot of Columns of Matrix V');
grid on;

%% Section3


clc; close all;

for i = 1:3

    subplot(3,1,i)
    plot(U(:,i));

end

sgtitle('3 columns of matrix U', ...
        'FontSize', 11, 'FontWeight', 'bold');

figure
stem([S(1,1) S(2,2) S(3,3)], 'LineWidth', 1)
title('Eigenspectrum');
grid on;

%% Section4

clc; close all;

S(1,1) = 0;
S(3,3) = 0;

sig_separated_svd  = U * S * V';
plot3ch(sig_separated_svd );



