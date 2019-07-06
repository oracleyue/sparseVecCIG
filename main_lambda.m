% Code to generate multivariate white nosie and test the proposed
% algorithm.

% Copyright [2019] <oracleyue>
% Last modified on 18 Jun 2019


clear all; close all;

% search paths
addpath('./Goran');  % if using Goran's algorithm

% init seed
rng(2);

% Data
% load('./Goran/Omega_Goran.mat');
p = 10;
dL = randi(5, p, 1)*3;
Omega = sprandOm(dL, [.3 .8]);
Sigma = inv(Omega);
d = sum(dL);
N = 10 * d;
X = mvnrnd(zeros(N,d), Sigma);
S = cov(X, 1);  % sample cov, normalized by N

% Setup
lambdaList = logspace(-2, 0, 40);  % range of lambdas
algType = 'zyue';    % choose algorithm
icType = 'BIC';      % choose information criterion
tolOpt = [1e-4, 50]; % "options" from "bcdSpML.m"

% Estimation
algTimer = tic;
[lambda, OmegaHat, vecIC, ~] = calcLambda(S, dL, N, lambdaList, ...
                                          icType, algType, tolOpt);
toc(algTimer)

% Visualization
% information criterion curves
figure
semilogx(lambdaList, vecIC, 'o-')
xlabel('lambda'); ylabel(icType)

% matrix plot of Omega and its estimation
figure
set(gcf,'color','white');
subplot(1,2,1)
imshowOm(Omega, 'spy');
title('Original $\mathbf{\Omega}$', ...
      'FontSize', 15, 'Interpreter','latex');
subplot(1,2,2)
imshowOm(OmegaHat, 'spy');
title('Estimation $\mathbf{\widehat{\Omega}}$', ...
      'FontSize', 15, 'Interpreter','latex');
