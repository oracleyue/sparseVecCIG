% Code to generate multivariate white nosie and test the proposed
% algorithm.

% Copyright [2019] <oracleyue>
% Last modified on 18 Jun 2019


clear all; close all;
addpath('./extern');

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
lambda = .2;

% perform estimationr 3, K2, to be nonspar
[OmegaHat, SigmaHat] = bcdpML(S, dL, lambda, 1e-6);

% visualization
figure
set(gcf,'color','white'); bt = 1;

subplot(1,2,1)
imOm = mat2gray(full(abs(Omega)));
imbOm = addborder(2*imOm, 1, 1, 'outer');
imshow(imbOm);
colormap(1-colormap('gray'));
title('Original $\mathbf{\Omega}$', ...
      'FontSize', 15, 'Interpreter','latex');

subplot(1,2,2)
imOmHat = mat2gray(abs(OmegaHat));
imbOmHat = addborder(2*imOmHat, 1, 1, 'outer');
imshow(imbOmHat);
colormap(1-colormap('gray'));
title('Estimation $\mathbf{\widehat{\Omega}}$', ...
      'FontSize', 15, 'Interpreter','latex');
