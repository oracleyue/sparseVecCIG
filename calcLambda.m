function lambda = calcLambda(S, dL, lambdaRange, type, epsilon)
% CALCLAMBDA computes the BIC or AIC criterion to determine lambda,
% i.e. the regularization parameter for l0-penalized ML.
%
% INPUT:
%   S      :   (d x d) sample covariance matrix, normalised by T
%   dL     :   (p x 1) vector of positive integers, and Sum(dL) = d
%   lambda :   positive real number, regularisation parameter
%   epsilon:   positive value close to 0; tolerance to stop iteration
%   (The above arguments inherit from "bcdpMLcg.m")
%   type   :   string (default: 'BIC'); set 'AIC' or 'BIC'
%
% OUTPUT:
%   lambda :   positive real value, the first one that minimises BIC values.

% Copyright [2019] <oracleyue>
% Last modified on 24 Jun 2019


% search paths
addpath('./tools');

% data
d = sum(dL);
T = 10 * d;
X = mvnrnd(zeros(T,d), Sigma);
% sample covariance, normalized by T
S = cov(X, 1);

% Generate lamba's range (use log scale)
numL = 40;
lambdas = logspace(-2, 0, numL);
% Initialize criterion values
KL = zeros(numL, 1);
BIC = zeros(numL, 1);
AIC = zeros(numL, 1);

% Perform estimation and compute criterion values
fprintf('Computing KL and BIC values:\n')
algTimer = tic;
for k = 1:numL
    lambda = lambdas(k);
    [OmegaHat, SigmaHat] = bcdpMLcg(S, dL, lambda, 1e-12);

    % refer to /Goran Marjanovic & Victor Solo. ICASSP, 2018/
    KL(k) = trace(Sigma*OmegaHat) - log(det(Sigma*OmegaHat)) - d;
    BIC(k) = trace(S*OmegaHat) - log(det(OmegaHat)) + ...
             log(T)/T * l0norm(OmegaHat, dL);
    AIC(k) = trace(S*OmegaHat) - log(det(OmegaHat)) + ...
             1/T * l0norm(OmegaHat, dL);

    fprintf('    %2d-th iterate: lambda=%.4f, BIC=%.4f \n', ...
            k, lambda, BIC(k));
end
toc(algTimer)
fprintf('End.\n')

% visualization
figure
subplot(1,3,1)
semilogx(lambdas, KL, 'o-')
xlabel('lambdas'); ylabel('KL')
subplot(1,3,2)
semilogx(lambdas, BIC, 'o-')
xlabel('lambdas'); ylabel('BIC')
subplot(1,3,3)
semilogx(lambdas, AIC, 'o-')
xlabel('lambdas'); ylabel('AIC')