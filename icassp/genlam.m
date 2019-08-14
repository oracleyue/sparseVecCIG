% Code to generate multivariate white nosie and test the proposed
% algorithm.

% Copyright [2019] <oracleyue>
% Last modified on 18 Jun 2019


clear all; close all;

addpath('../');       % project root
addpath('../Goran');  % Goran's algorithm

% init
rng(2);

% data
load('../Goran/Omega_Goran.mat');
% p = 100;
% dL = randi(5, p, 1)*3;
% Omega = sprandOm(dL, [.3 .8]);
Sigma = inv(Omega);
d = sum(dL);
N = 10 * d;
X = mvnrnd(zeros(N,d), Sigma);
S = cov(X, 1);  % sample cov, normalized by N

% Setup
lambdaList = logspace(-2, 0, 40);  % range of lambdas
algType = 'goran';    % choose algorithm
icType = 'BIC';      % choose information criterion
tolOpt = [1e-3, 10]; % "options" from "bcdSpML.m"

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
title('Ground truth $\mathbf{\Omega}$', ...
      'FontSize', 12, 'Interpreter','latex');
subplot(1,2,2)
imshowOm(OmegaHat, 'spy');
title('Estimation $\mathbf{\widehat{\Omega}}$', ...
      'FontSize', 12, 'Interpreter','latex');
