% This is the problem file for solving the 5 problems related to the
% example.

%% Problem 1: Basic Eigenvalue and Eigenvector calculation.

clc, clearvars

B = [5 4; 2 3]; % Matrix B
[V,D] = eig(B); % eigenvalues and eigenvectors.

% show the results
disp('PROBLEM 1 FINAL SOLUTIONS:')
disp('Eigenvalues:')
disp(diag(D)); %eigen values

disp('Eigenvectors:')
disp(V);

%% Problem 2: Eigenvalue Decomposition for Matrix exponentiation

clc, clearvars

B = [1 2 3; 0 1 4; 5 6 0] % matrix B
[V,D] = eig(B) % eigen calculations

t = 2 % time variable

exp_Bt = V * diag(exp(diag(D) * t)) * inv(V) % matrix exponential using eigenvalue decomposition

disp('PROBLEM 2 FINAL SOLUTIONS:');
disp('Matrix Exponential e^Bt:');
disp(exp_Bt);


%% Problem 3: Schur Decomposition
% Perform the Schir decomposition of matrix B using MATLAB

clc,clearvars

% matrix B
B = [4 1 2; -1 3 0; 1 1 2]

% Schur decomp
[U S] = schur(B) % (Orthogonal Matrix U) & (Upper Triangular Matrix S)

% Display the results
disp('PROBLEM 3 FINAL SOLUTIONS:')
disp('Orthogonal Matrix U:');
disp(U);

disp('Upper Triangular Matrix S:')
disp(S);

%% Problem 4: Stability Analysis using Eigenvalues
% given matrix B, determine if the system is stable using its eigenvalues

clc, clearvars

B = [-2 1; -1 3] % Matrix B

lambda = eig(B) % eigenvalues

% results
disp('PROBLEM 4 FINAL SOLUTIONS:')
disp('Eigenvalues:')
disp(lambda);

% Check stability based on eigenvalues
if all(real(lambda) < 0)
    disp('[+] The system is stable.')
else
    disp('[-] The system is unstable.')
end

%% Problem 5: Visualizing eigenvectors
% Visualize the eigenvectors of the Matrix B using MATLAB.

clc, clearvars

B = [3 1; 1 3] % B matrix

[V D] = eig(B) % eigenvalues

% Plotting original and eigenvectors.
figure;
hold on;

quiver(0, 0, 1, 0, 'r', 'LineWidth', 2); % x
quiver(0, 0, 0, 1, 'b', 'LineWidth', 2); % y
quiver(0, 0, V(1, 1), V(2, 1), 'g', 'LineWidth', 2); % First eigenvector
quiver(0, 0, V(1, 2), V(2, 2), 'm', 'LineWidth', 2); % Second eigenvector
axis equal;
legend('x-axis', 'y-axis', 'First eigenvector', 'Second eigenvector');
title('Visualizations of Eigenvectors from Matrix B');
hold off;
























