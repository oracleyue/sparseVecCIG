% Code to generate multivariate white nosie and test the proposed
% algorithm.

% Copyright [2019] <oracleyue>
% Last modified on 18 Jun 2019


clear all; close all;

% search paths
addpath('./goran');  % if using Goran's algorithm

% init seed (goran's demo: 2)
rng(2);

% Data
load('./data/Omega_Goran.mat'); p = 5;
% p = 6;
% dL = randi(5, p, 1)*3;
% Omega = sprandOm(dL, [.3 .8]);
Sigma = inv(Omega);
d = sum(dL);
N = round(.5 * d);
X = mvnrnd(zeros(N,d), Sigma);
S = cov(X, 1);  % sample cov, normalized by N

% Setup
lambdaList = logspace(log10(3), log10(4), 40);    % range of lambdas
algType = 'zyue';    % choose algorithm
icType = 'BIC';      % choose information criterion
perm = [];
% rng('shuffle'); perm = randperm(p);  % randomize iteration order
algOpt = setOptions('perm', perm, 'initType', 'fixed', ...
                            'precision', [1e-4, 100], ...
                            'errorType', {'rel', 'var'});
% if use all default, simply run "algOpt = setOptions()".

% Estimation
algTimer = tic;
[lambda, OmegaHat, vecIC, ~] = calcLambda(S, dL, N, lambdaList, ...
                                          icType, algType, algOpt);
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
