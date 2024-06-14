%% Historical return data for 5 different stocks
clc

% Historical return data for 5 different stocks
returns = [
    0.02, 0.04, -0.01, 0.05, 0.03;
    0.01, 0.02, 0.00, 0.04, 0.05;
    -0.02, 0.03, 0.01, 0.01, -0.02;
    0.04, -0.01, 0.02, 0.03, 0.04;
    0.03, 0.00, 0.04, -0.01, 0.03
]

% Data loaded for the 5 different stocks, so problem 1 solved!

%% Problem 2: Calculate the covariance matrix of the returns

clc

%cov-matrix calculated
covMatrix = cov(returns)

% problem 2 solved - covMatrix is the covariance matrix calculated from
% returns matrix.

%% Problem 3: Use the Power Iteration method to find the dominant eigenvalue and its eigenvector

clc

% PIM to find the dominant eigenvalue and eigenvector
n = size(covMatrix, 1);
x = rand(n, 1); % Initial guess
tol = 1e-6;
maxIter = 1000;
lambda_old = 0;

for k = 1:maxIter 
    x = covMatrix * x;
    x = x / norm(x);
    lambda_new = x' * covMatrix * x;
    if abs(lambda_new - lambda_old) < tol 
        break;
    end
    lambda_old = lambda_new;
end

dominanteigenvalue = lambda_new
dominantEigenvector = x

%% Problem 4: Use the inverse Iteration Method to find the smallest eigenvalue approximation

clc

mu = 0; % Initial shift close to zero
x = rand(n, 1); % Initial guess
B = covMatrix - mu * eye(n);

for k = 1:maxIter
    y = B \ x;
    x = y / norm(y);
    mu = x' * covMatrix * x;
end

smallestEigenvalue = mu
smalelstEigenvector = x

%% Problem 5: Use the Rayleigh Quotient Iteration to refine the eigenvalue approximation
clc

x = rand(n, 1); % initial guess
lambda = x' * covMatrix * x;

for k = 1:maxIter
    B = covMatrix - lambda * eye(n);
    y = B\x;
    x = y / norm(y);
    lambda = x' * covMatrix * x;
    if norm(covMatrix * x - lambda * x) < tol
        break;
    end
end

rayleighEigenvalue = lambda
rayleighEigenvector = x

%% Interpretation

clc

disp('Largest variance direction in the portfolio:');
disp(dominanteigenvalue);
disp('Smallest variance direction in the portfolio:');
disp(smallestEigenvalue);
disp('Rayleighs refined eigenvalue estimation: ');
disp(rayleighEigenvalue);