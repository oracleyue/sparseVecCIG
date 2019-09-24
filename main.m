% Code to generate multivariate white nosie and test the proposed
% algorithm.

% Copyright [2019] <oracleyue>
% Last modified on 18 Jun 2019


clear all; close all;

% search paths
addpath('./goran');  % if using Goran's algorithm

% init
rng(2);

% data
load('./data/Omega_Goran.mat');
% p = 6;
% dL = randi(5, p, 1)*3;
% Omega = sprandOm(dL, [.3 .8]);
Sigma = inv(Omega);
d = sum(dL);
N = round(.5 * d);
X = mvnrnd(zeros(N,d), Sigma);
S = cov(X, 1);  % sample cov, normalized by N
% setup
lambda = 4;   % 0.15 ~ 0.3
algName = 'zyue';
precision = [1e-4 50];
errorType = {'rel', 'var'};
perm = [];
% rng('shuffle'); perm = randperm(length(dL));

% estimation
algTimer = tic;
switch algName
  case 'zyue'
    % CG-embeded solver
    [OmegaHat, ~, optStatus] = spMLE(S, dL, lambda, perm, ...
                                     'precision', precision, ...
                                     'errorType', errorType);
    disp('Optimization status:'); disp(optStatus)
  case 'goran'
    % Goran's solver (use l1 penalty)
    [~, ~, OmegaHat] = Algorithm(speye(d), S, dL, lambda, ...
                                 precision(2), precision(1), 1);
end
toc(algTimer)

% visualization
figure
set(gcf,'color','white');
subplot(1,2,1)
imshowOm(Omega, 'spy');
title('Original $\mathbf{\Omega}$', ...
      'FontSize', 12, 'Interpreter','latex');
subplot(1,2,2)
imshowOm(OmegaHat, 'spy');
title('Estimation $\mathbf{\widehat{\Omega}}$', ...
      'FontSize', 12, 'Interpreter','latex');
