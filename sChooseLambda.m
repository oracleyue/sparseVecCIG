% The script computes the BIC criterion to determine lambda, the
% regularization parameter.

% Copyright [2019] <oracleyue>
% Last modified on 24 Jun 2019


clear all; close all;
rng(2);

% use Goran's data
% addpath('./Goran');
% load('./Goran/Omega_Goran.mat');

% random data generation
p = 6;
dL = randi(9, p, 1)*3;
Omega = sprandOm(dL, [.3 .8]);
Sigma = inv(Omega);
d = sum(dL);
N = 10 * d;  % length of data
X = mvnrnd(zeros(N,d), Sigma);
S = cov(X, 1);  % sample cov, normalized by N

% generate lamba's range (use log scale)
numL = 40;
lambdas = logspace(-2, 0, numL);
BIC = zeros(numL, 1);

% perform estimation and compute criterion values
[lambda, BICs, eTime] = calcLambda(S, dL, N, lambdas, 'BIC');

% visualization
fig_hl = figure;
semilogx(lambdas, BICs, 'o-')
hold on
plot(lambda, BICs(find(lambdas==lambda)), '*', 'MarkerSize', 12);
xlabel('Lambda'); ylabel('BIC')