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
load('./Goran/Omega_Goran.mat');
% load('~/Workspace/data/data_p500d4449.mat');
% p = 100;
% dL = randi(5, p, 1)*3;
% Omega = sprandOm(dL, [.3 .8]);
Sigma = inv(Omega);
d = sum(dL);
N = 10 * d;
X = mvnrnd(zeros(N,d), Sigma);
S = cov(X, 1);  % sample cov, normalized by N

% setup
lambda = 0.5;   % 0.15 ~ 0.3
algName = 'zyue';
algOpt = {[1e-3 50], 'rel', 'var'};

% estimation
algTimer = tic;
switch algName
  case 'zyue'
    % CG-embeded solver
    [OmegaHat, ~, optStatus] = bcdSpML(S, dL, lambda, algOpt);
    disp('Optimization status:'); disp(optStatus)
  case 'goran'
    % Goran's solver
    [~, ~, OmegaHat] = Algorithm(speye(d), S, dL, lambda, 10, 1e-3, 0);
end
toc(algTimer)

%% visualization
figure
set(gcf,'color','white');
subplot(1,2,1)
imshowOm(Omega, 'raw');
title('Original $\mathbf{\Omega}$', ...
      'FontSize', 12, 'Interpreter','latex');
subplot(1,2,2)
imshowOm(OmegaHat, 'raw');
title('Estimation $\mathbf{\widehat{\Omega}}$', ...
      'FontSize', 12, 'Interpreter','latex');
