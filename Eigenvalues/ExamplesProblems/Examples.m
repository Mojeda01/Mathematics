%% EXAMPLE 1: Basic Eigenvalue and Eigenvector Calculation
% This example demonstrates how to calculate the eigenvalues and
% eigenvectors of a given matrix using the 'eig' function.

clc, clearvars

% define the matrix A
A = [4, -2; 1, 1];

% Calculate the eigenvalues and eigenvectors
[V,D]=eig(A);

% Display the results
disp('Eigenvalues:')
disp(diag(D)); % The eigenvalue are on the diagonal of D
disp('Eigenvectors:');
disp(V);

% Explanation:
% The 'eig' function is used to compute the eigenvalues and eigenvectors of
% A.
% V contains the eigenvectors, and D is a diagonal matrix with the
% eigenvalues on the diagonal.

%% EXAMPLE 2 - Eigenvalue Decomposition for Matrix Exponentiation
% This example shows how to use eigenvalue decomposition to compute the
% matrix exponential e^{At} for solving a differential equation.

clc, clearvars

% Define matrix A
A = [0 -6 -1; 6 2 -16; -5 20 -10];

% Calculate the eigenvalue decomposition
[V,D]=eig(A);

% Ddefine a time variable t
t = 1;

% Compute the matrix exponential using eigenvalue decomposition
exp_At = V * diag(exp(diag(D) * t)) * inv(V);

% Display the result
disp('Matrix Exponential e^At:');
disp(exp_At);

% Explanation:
% Eigenvalue decomposition (A = V*D*V^-') is used to compute the matrix
% exponent. We compute e^D (where D is diagonal) by exponentiating the
% diagonal elements, then reconstruct e^At.

%% Example 3: Schur Decomposition
% This example demonsrates the use of the 'schur' functiomn to perform
% schur decomposition, which is useful for matrices that do not have a full
% set of linearly independent eigenvectors.

clc, clearvars

% Define the matrix A
A = [6 12 19; -9 -20 -33; 4 9 15];

% Perform Schur decomposition
[U,S] = schur(A);

% Display the results
disp('Orthogonal Matrix U:');
disp(U);

disp('Upper Triangular Matrix S:');
disp(S);

% Explanation:
% The 'schur' function performs the schur decomposition, which factorizes A
% into U*S*U'.
% U is an orthogonal matrix, and S is an upper triangular matrix.

%% Example 4: Stability Analysis using Eigenvalues
% This example uses eigenvalues to analyze the stability of a system
% represented by a matrix.

clc, clearvars

% Define matrix A (system matrix)
A = [0 1; -2 -3];

% Calculate the eigenvalues
lambda = eig(A);

% Display the eigenvalues
disp('Eigenvalues:');
disp(lambda);

% Check stability based on eigenvalues
if all(real(lambda) < 0)
    disp('The system is stable.');
else
    disp('The system is unstable');
end

% Explanation:
% The 'eig' function is used to compute the eigenvalues of the system
% matrix A. Stability is determined by checking if the real parts of all
% eigenvalues are negative.

%% Example 5: Visualizing Eigenvectors
% This example visualizes the eigenvectors of a matrix to illustrate their
% direction.

clc, clearvars

% Define the matrix A
A = [2 1; 1 2];

% Calculate the eigenvalues and eigenvectors
[V,D] = eig(A);

% Plot the original vectors and eigenvectors
figure;
hold on;

quiver(0, 0, 1, 0, 'r', 'LineWidth', 2); % x-axis
quiver(0, 0, 0, 1, 'b', 'LineWidth', 2); % y-axis
quiver(0, 0, V(1, 1), V(2, 1), 'g', 'LineWidth', 2); % First eigenvector
quiver(0, 0, V(1, 2), V(2, 2), 'm', 'LineWidth', 2) % Second eigenvector
axis equal;
legend('x-axis', 'y-axis', 'First eigenvector', 'Second eigenvector');
title('Visualization of Eigenvectors');
hold off;

% Explanation:
% The 'eig' function is used to compute the eigenvalues and eigenvectors.
% The 'quiver' function is used to plot vectors. Eigenvectors are
% visualized to show their direction.