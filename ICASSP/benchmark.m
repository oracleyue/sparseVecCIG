% The script is to benchmark speeds of two implementations (CG-embeded
% vs. Goran's), by running the algorithms for a range of lambda.
%
% Notes: It may not be reasonable to benchmark time costs by running with an
% arbitrary fixed lambda, since the choice of lambda might effect the
% convergence of algorithms and hence time costs.

% Copyright [2019] <oracleyue>
% Last modified on 14 Aug 2019


clear all; close all;

% Search paths
addpath('../');       % project root
addpath('../goran');  % using Goran's

% Init seed
rng(2)

%% Benchmark for different dimensions
pList = 4:2:20;
dList = zeros(size(pList));
dLCell = cell(size(pList));
for k = 1:length(pList)
    if k == 1
        dLCell{k} = randi(5, pList(1), 1)*3;
    else
        dLCell{k} = [dLCell{k-1}; randi(5, pList(k)-pList(k-1), 1)*3];
    end
end

% set lambda range
numL = 40;
lambdaList = logspace(-2, 0, numL);

% init saving variables
lambdaBest = zeros(2, length(pList));     % 1: 'zyue'; 2: 'goran'
infoCrit = zeros(2*numL, length(pList));  % 1:numL: 'zyue';
                                          % numL+1:2*numL: 'goran'
eTimeMatrix = zeros(2*numL, length(pList));
eTime = zeros(2, length(pList));

fprintf('Speed benchmark (elapsed CPU time):\n')
for k = length(pList)
    % data
    p = pList(k);
    dL = dLCell{k};
    rng(2)
    Omega = sprandOm(dL, [.3 .8]);
    Sigma = inv(Omega);
    d = sum(dL); dList(k) = d;
    N = 10 * d;
    X = mvnrnd(zeros(N,d), Sigma);
    % sample covariance, normalized by T
    S = cov(X, 1);

    fprintf('  [%2d]: #blocks=%2d, dim=%3d\n', k, p, d);

    tolOpt = [1e-3, 10];
    % ML with different lambdas
    [lambdaZ, ~, ICsZ, eTimeZ] = calcLambda(S, dL, N, lambdaList, ...
                                            'BIC', 'zyue', tolOpt);
    fprintf('        CG   : %.6fs \n', sum(eTimeZ));
    [lambdaG, ~, ICsG, eTimeG] = calcLambda(S, dL, N, lambdaList, ...
                                            'BIC', 'goran', tolOpt);
    fprintf('        Goran: %.6fs \n', sum(eTimeG));

    lambdaBest(:,k) = [lambdaZ; lambdaG];
    infoCrit(1:numL,k) = ICsZ;
    infoCrit(numL+1:2*numL,k) = ICsG;
    eTimeMatrix(1:numL, k) = eTimeZ;
    eTimeMatrix(numL+1:2*numL, k) = eTimeG;
    eTime(1, k) = eTimeZ(find(lambdaList == lambdaZ));
    eTime(2, k) = eTimeG(find(lambdaList == lambdaG));
end
fprintf('End.\n')

% save benchmark results
save('bm_results.mat');
