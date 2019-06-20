% Code to generate multivariate white nosie and test the proposed
% algorithm.

% Copyright [2019] <oracleyue>
% Last modified on 18 Jun 2019


clear all; close all;

% init seed
rng(2);

% dimensions
dL = [2, 3, 5, 1, 3]*3;
d = sum(dL);

% data
Omega = sparse(genOmega(dL));
Sigma = inv(Omega);
T = 10 * d;
X = mvnrnd(zeros(T,d), Sigma);
% sample covariance, normalized by T
S = cov(X, 1);
% lambda
lambda = 0.0522;

% perform estimation
[OmegaHat, SigmaHat] = bcdpML(S, dL, lambda, 1e-3);

% visualization
figure
set(gcf,'color','white'); bt = 1;

subplot(1,2,1)
imgOm = mat2gray(abs(Omega));
imshow(imgOm, 'InitialMagnification', 250);
colormap(1-colormap);
title('Original $\mathbf{\Omega}$', ...
      'FontSize',20, 'Interpreter','latex');

subplot(1,2,2)
imgOmHat = mat2gray(abs(OmegaHat));
imshow(imgOmHat, 'InitialMagnification',250);
colormap(1-colormap);
title('Estimator $\mathbf{\widehat{\Omega}}_{l_0}$', ...
      'FontSize', 20, 'Interpreter','latex');
