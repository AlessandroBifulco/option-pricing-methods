function [V0] = LSM(S0,K,T,r,sigma,numPaths,numSteps)
% Input:
% S0        = initial stock price
% K         = strike price
% T         = time to maturity
% r         = risk-free interest rate
% sigma     = volatility
% numPaths  = number of simulated paths
% numSteps  = number of time steps
%
% Output:
% V0        = estimated price of an American put option using the
%             Longstaff–Schwartz least-squares Monte Carlo method

% Time step and discount factor
dt = T/numSteps;
discountFactor = exp(-r*dt);

% Random shocks (antithetic variates)
Z = randn(numSteps,numPaths);
increments = exp((r - 0.5*sigma^2)*dt + sigma*sqrt(dt) * [Z, -Z]);

% Simulated price paths
paths = cumprod([S0*ones(1,2*numPaths); increments]);

% Terminal payoff (American put)
cashflows = max(K - paths(end,:), 0);

% Backward induction over exercise opportunities
for step = numSteps:-1:2
    ITM_idx = find(paths(step,:) < K); % In-the-money paths
    if isempty(ITM_idx)
        continue;
    end

    spots_ITM = paths(step, ITM_idx);
    discountedCF_ITM = cashflows(ITM_idx) * discountFactor;

    % Regression basis: 1, S, S^2
    regMatrix = [ones(size(spots_ITM)); spots_ITM; spots_ITM.^2]';
    regTargets = discountedCF_ITM';
    beta = regMatrix \ regTargets;

    % Estimated continuation value
    continuationValue = beta' * [ones(size(spots_ITM)); spots_ITM; spots_ITM.^2];

    % Immediate exercise value
    exerciseValue = K - spots_ITM;

    % Exercise decision
    exerciseIdx = ITM_idx(exerciseValue > continuationValue);

    cashflows(exerciseIdx) = exerciseValue(exerciseValue > continuationValue);
    nonExerciseIdx = setdiff(1:2*numPaths, exerciseIdx);
    cashflows(nonExerciseIdx) = cashflows(nonExerciseIdx) * discountFactor;
end

% Discounted expected value
V0 = mean(cashflows * discountFactor);
end

