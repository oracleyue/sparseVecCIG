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
% addpath('../goran');  % using Goran's codes

% Init seed
rng(2)

%% Benchmark for different dimensions
% pListZ = 4:2:20;
pListZ = [4:2:20 30:10:100 150 200 250];
dListZ = zeros(size(pListZ));
dLCellZ = cell(size(pListZ));
for k = 1:length(pListZ)
    if k == 1
        dLCellZ{k} = randi(5, pListZ(1), 1)*3;
    else
        dLCellZ{k} = [dLCellZ{k-1}; randi(5, pListZ(k)-pListZ(k-1), 1)*3];
    end
end

% set lambda range
numL = 40;
lambdaList = logspace(-2, 0, numL);

% init saving variables
lambdaBestZ = zeros(1, length(pListZ));     % 1: 'zyue'; 2: 'goran'
infoCritZ = zeros(numL, length(pListZ));  % 1:numL: 'zyue';
                                            % numL+1:2*numL: 'goran'
eTimeMatrixZ = zeros(numL, length(pListZ));
eTimeZ = zeros(1, length(pListZ));

fprintf('Speed benchmark (elapsed CPU time):\n')
for k = 1:length(pListZ)
    % data
    p = pListZ(k);
    dL = dLCellZ{k};
    rng(2)
    Omega = sprandOm(dL, [.3 .8]);
    Sigma = inv(Omega);
    d = sum(dL); dListZ(k) = d;
    N = 10 * d;
    X = mvnrnd(zeros(N,d), Sigma);
    % sample covariance, normalized by T
    S = cov(X, 1);

    fprintf('  [%2d]: #blocks=%2d, dim=%3d\n', k, p, d);

    algOpt = setOptions('precision', [1e-3, 10]);
    % ML with different lambdas
    [lambdaZ, ~, ICsZ, eTime] = calcLambda(S, dL, N, lambdaList, ...
                                            'BIC', 'zyue', algOpt);
    fprintf('        CG   : %.6fs \n', sum(eTimeZ));

    % [lambdaG, ~, ICsG, eTimeG] = calcLambda(S, dL, N, lambdaList, ...
    %                                         'BIC', 'goran', algOpt);
    % fprintf('        Goran: %.6fs \n', sum(eTimeG));
    % [Note]:
    % One run of Goran's for 40 lambda for case p = 20 has costed more
    % than 10 hours. Thus, we stop recording its time costs.

    lambdaBestZ(k) = lambdaZ;
    infoCritZ(:,k) = ICsZ;
    eTimeMatrixZ(:, k) = eTime;
    eTimeZ(k) = eTime(find(lambdaList == lambdaZ));
end
fprintf('End.\n')

% save benchmark results
% save('bm_results_zyue.mat', 'pListZ', 'dLCellZ', 'dListZ', ...
%      'lambdaBestZ', 'infoCritZ', 'eTimeMatrixZ', 'eTimeZ');
