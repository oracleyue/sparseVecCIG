% Code to generate multivariate white nosie and test the proposed
% algorithm.

% Copyright [2019] <oracleyue>
% Last modified on 18 Jun 2019


clear all; close all;

% search paths
addpath('./extern');
addpath('./Goran');  % if using Goran's algorithm

% flags
autoLambdaFlag = 1;

% init seed
rng(2);

% data
% load('./Goran/Omega_Goran.mat');
p = 8;
dL = randi(5, p, 1)*3;
Omega = sprandOm(dL, [.3 .8]);
Sigma = inv(Omega);
d = sum(dL);
N = 10 * d;
X = mvnrnd(zeros(N,d), Sigma);
% sample covariance, normalized by N
S = cov(X, 1);
% lambda
if autoLambdaFlag
    lambdaList = logspace(-2, 0, 40);
    bicTimer = tic;
    lambda = calcLambda(S, dL, N, lambdaList, 'AIC');
    toc(bicTimer)
else
    % fix a lambda
    lambda = .2;
end
% choose algorithm
algName = 'zyue';

% perform estimation
algTimer = tic;
switch algName
  case 'zyue'
    % CG-embeded solver
    [OmegaHat, ~] = bcdSpML(S, dL, lambda, 1e-12);
  case 'goran'
    % Goran's solver
    [~, ~, OmegaHat] = Algorithm(speye(d), S, dL, lambda, 10, 1e-3, 0);
end
toc(algTimer)

% save('result_large.mat');

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
