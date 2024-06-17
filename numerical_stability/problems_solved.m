
%% Problem 1: Load the historical returns for 5 different stocks
clc, clearvars

returns = [
    0.02, 0.04, -0.01, 0.05, 0.03;
    0.01, 0.02, 0.00, 0.04, 0.05;
    -0.02, 0.03, 0.01, 0.01, -0.02;
    0.04, -0.01, 0.02, 0.03, 0.04;
    0.03, 0.00, 0.04, -0.01, 0.03
];

disp('Returns:');
disp(returns);

%% Problem 2: Calculate the covariance matrix of the returns

clc;

covMatrix = cov(returns);

disp('Covariance Matrix:');
disp(covMatrix);

%% Problem 3: Calculate the condition number using differnet norms 
% (2-norm, infinity norm, and 1-norm)

clc;

cond_2 = cond(covMatrix); % 2-norm (default)
cond_inf = cond(covMatrix, 'inf'); % Infinity norm
cond_1 = cond(covMatrix,1); % 1-norm

disp('2-Norm:');
disp(cond_2)
disp('Infinity Norm:');
disp(cond_inf)
disp('1-Norm');
disp(cond_1)

%% Problem 4: Analyze numerical stability by perturbing the return data slightly and observing the changes in the covariance matrix

clc;

% Analyze numerical stability by pertubing the data by adding a small
% pertubation to the returns.
perturbed_returns = returns + 0.001 * randn(size(returns));
% perturbed covariance matrix
perturbed_matrix = cov(perturbed_returns);
% calculate the change in the covariance matrix
change_in_covMatriax = norm(perturbed_matrix - covMatrix, 'fro');

disp('Perturbed returns');
disp(perturbed_returns);
disp('Perturbed matrix');
disp(perturbed_matrix);
disp('Change in covMatrix');
disp(change_in_covMatriax);