% The script computes the BIC criterion to determine lambda, the
% regularization parameter.

% Copyright [2019] <oracleyue>
% Last modified on 24 Jun 2019


clear all; close all;

% Search paths
addpath('./extern'); % external libraries
addpath('./Goran');  % using Goran's
addpath('./tools');

% Init seed
rng(2);

% Data
% load('./Goran/Omega_Goran.mat');
p = 6;
dL = randi(9, p, 1)*3;
Omega = sprandOm(dL, [.3 .8]);
Sigma = inv(Omega);
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