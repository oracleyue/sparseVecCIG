% Code to generate multivariate white nosie and test the proposed
% algorithm.

% Copyright [2019] <oracleyue>
% Last modified on 14 Aug 2019


clear all; close all;

% search paths
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

% setup
lambda0 = 0.24245;
lambda1 = 0.046416;
algOpt = {[1e-3 10], 'rel', 'var'};

% estimation
algTimer = tic;
% quick solver
[OmHatL0, ~, optStatus] = bcdSpML(S, dL, lambda0, algOpt);
disp('Optimization status:'); disp(optStatus)
% Goran's solver for l1
[~, ~, OmHatL1] = Algorithm(speye(d), S, dL, lambda1, 10, 1e-3, 1);

% visualization
fig_hl = figure;
set(gcf,'color','white');
subplot(1,3,1)
imshowOm(Omega, 'raw');
title('Ground truth $\mathbf{\Omega}$', ...
      'FontSize', 12, 'Interpreter','latex');
subplot(1,3,2)
imshowOm(OmHatL0, 'raw');
title('Estimation $\mathbf{\widehat{\Omega}}_{l_0}$', ...
      'FontSize', 12, 'Interpreter','latex');
subplot(1,3,3)
imshowOm(OmHatL1, 'raw');
title('Estimation $\mathbf{\widehat{\Omega}}_{l_1}$', ...
      'FontSize', 12, 'Interpreter','latex');

filename = ['./goran-demo', '.pdf'];
pos = [2.0104 4.4688 8.1354 2.6771];
set(fig_hl,'Units','Inches', 'Position', pos);
set(fig_hl,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(fig_hl, filename, '-dpdf', '-r0')
