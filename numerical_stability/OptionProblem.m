%% (1) generate synthetic historical return data for 5 different stocks over 1000 days
numDays = 1000;
numStocks = 5
returns = 0.001 * randn(numDays, numStocks);

%% (2) Calculate the mean and covaraince of the returns
meanReturns = mean(returns);
covMatrix = cov(returns);

%% (3) Perform monte carlo

numSimulations = 10000;
T = 1; 
r = 0.05;
S0 = 100;
K = 100;

% Preallocate array for simulated end prices
endPrices = zeros(numSimulations, numStocks);

% sim end prices using geometric brownian motion
for i = 1:numSimulations
    Z = randn(1, numStocks);
    endPrices(i,:) = S0 .* exp((r - 0.5 * diag(covMatrix)') * T + sqrt(T) * Z * chol(covMatrix));
end

% Calculate the basket prices
basketPrices = mean(endPrices, 2);

% Calculate the payoff for a call option
payoffs = max(basketPrices - K, 0)

% Discount the payoffs to present value
basketOptionPrice = exp(-r * T) * mean(payoffs)

%% (4) Calculate the condition number of the covairance matrix
