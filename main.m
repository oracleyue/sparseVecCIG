% Code to generate multivariate white nosie and test the proposed
% algorithm.

% Copyright [2019] <oracleyue>
% Last modified on 18 Jun 2019


clear all; close all;

% search paths
addpath('./Goran');  % if using Goran's algorithm

% init
rng(2);

% data
% load('./Goran/Omega_Goran.mat');
% load('~/Workspace/data/data_p500d4449.mat');
p = 200;
dL = randi(5, p, 1)*3;
Omega = sprandOm(dL, [.3 .8]);
Sigma = inv(Omega);
d = sum(dL);
N = 10 * d;
X = mvnrnd(zeros(N,d), Sigma);
S = cov(X, 1);  % sample cov, normalized by N

% setup
lambda = 0.5;  % 0.15 ~ 0.3
algName = 'zyue';

% estimation
algTimer = tic;
switch algName
  case 'zyue'
    % CG-embeded solver
    [OmegaHat, ~, optStatus] = bcdSpML(S, dL, lambda, [1e-3 20]);
    disp('Optimization status:'); disp(optStatus)
  case 'goran'
    % Goran's solver
    [~, ~, OmegaHat] = Algorithm(speye(d), S, dL, lambda, 20, 1e-3, 0);
end
toc(algTimer)

% visualization
figure
set(gcf,'color','white');
subplot(1,2,1)
imshowOm(Omega, 'spy', dL);
title('Original $\mathbf{\Omega}$', ...
      'FontSize', 15, 'Interpreter','latex');
subplot(1,2,2)
imshowOm(OmegaHat, 'spy', dL);
title('Estimation $\mathbf{\widehat{\Omega}}$', ...
      'FontSize', 15, 'Interpreter','latex');
