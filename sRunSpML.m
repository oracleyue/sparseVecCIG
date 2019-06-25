% Code to generate multivariate white nosie and test the proposed
% algorithm.

% Copyright [2019] <oracleyue>
% Last modified on 18 Jun 2019


clear all; close all;

% search paths
addpath('./extern');
addpath('./Goran');  % if using Goran's algorithm

% init seed
rng(2);

% data
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
% lambda
lambda = .4924;  % by calcLambda.m
% choose algorithm
algName = 'zyue';

% perform estimationr 3, K2, to be nonspar
algTimer = tic;
switch algName
  case 'zyue'
    % CG-embeded solver
    [OmegaHat, SigmaHat] = bcdSpML(S, dL, lambda, 1e-12);
  case 'goran'
    % Goran's solver
    [F, Time, OmegaHat] = Algorithm(speye(d), S, dL, lambda, 10, 1e-3, 0);
end
toc(algTimer)

%% visualization
figure
set(gcf,'color','white'); bt = 1;

subplot(1,2,1)
imOm = mat2gray(full(abs(Omega)));
imbOm = addborder(2*imOm, bt, 1, 'outer');
imshow(imbOm, 'InitialMagnification','fit');
colormap(1-colormap('gray'));
title('Original $\mathbf{\Omega}$', ...
      'FontSize', 15, 'Interpreter','latex');

subplot(1,2,2)
imOmHat = mat2gray(full(abs(OmegaHat)));
imbOmHat = addborder(2*imOmHat, bt, 1, 'outer');
imshow(imbOmHat, 'InitialMagnification','fit');
colormap(1-colormap('gray'));
title('Estimation $\mathbf{\widehat{\Omega}}$', ...
      'FontSize', 15, 'Interpreter','latex');
