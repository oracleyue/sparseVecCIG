% Code to generate multivariate white nosie and test the proposed
% algorithm.

% Copyright [2019] <oracleyue>
% Last modified on 14 Aug 2019


clear all; close all;

% search paths
addpath('../');
addpath('../goran');
% Due to copyright issues, the github repository doesn't include Goran
% Marjanovic's codes (i.e. the folder "goran"). You may contact the
% authors in
% /Marjanovic, G., & Solo, V. (2018). Vector $l_0$ Sparse Conditional
%  Independence Graphs. 2018 IEEE International Conference on Acoustics,
%  Speech and Signal Processing (ICASSP), 2731–2735. IEEE./
% for their codes.

% init
rng(2);

% data
load('../data/Omega_Goran.mat');
Sigma = inv(Omega);
d = sum(dL);
N = 10 * d;
X = mvnrnd(zeros(N,d), Sigma);
S = cov(X, 1);  % sample cov, normalized by N

% setup
lambda0 = 0.5;
lambda1 = 0.14;
precision = [1e-3 10];

% estimation
algTimer = tic;
% quick solver
[OmHatL0, ~, optStatus] = spMLE(S, dL, lambda0, 'precision', precision);
disp('Optimization status:'); disp(optStatus)
% Goran's l0 solver (default)
[~, ~, OmHatL0g] = Algorithm(speye(d), S, dL, lambda0, ...
                             precision(2), precision(1), 0);
% Goran's solver using l1 modification
[~, ~, OmHatL1g] = Algorithm(speye(d), S, dL, lambda1, ...
                             precision(2), precision(1), 1);

% visualization
fig_hl = figure;
set(gcf,'color','white');
subplot(1,3,1)
imshowOm(OmHatL0, 'raw');
title('Estimation $\mathbf{\widehat{\Omega}}$', ...
      'FontSize', 12, 'Interpreter','latex');
subplot(1,3,2)
imshowOm(OmHatL0g, 'raw');
title('Estimation $\mathbf{\widehat{\Omega}}_{l_0}$', ...
      'FontSize', 12, 'Interpreter','latex');
subplot(1,3,3)
imshowOm(OmHatL1g, 'raw');
title('Estimation $\mathbf{\widehat{\Omega}}_{l_1}$', ...
      'FontSize', 12, 'Interpreter','latex');

filename = ['./goran-demo', '.pdf'];
pos = [2.0104 4.4688 8.1354 2.6771];
set(fig_hl,'Units','Inches', 'Position', pos);
set(fig_hl,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(fig_hl, filename, '-dpdf', '-r0')
