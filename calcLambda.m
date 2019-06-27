function [lambda, Omega, icVec, eTime] = calcLambda(S, dL, N, lambdaList, ...
                                                    icType, algType, tolOptions)
% CALCLAMBDA computes the BIC or AIC criterion to determine lambda,
% i.e. the regularization parameter for l0-penalized ML.
%
% INPUT:
%   S          :   (d x d) sample covariance matrix, normalised by N
%   dL         :   (p x 1) vector of positive integers, and Sum(dL) = d
%   N          :   length of data
%   lambdaList :   positive real vector of lambdas
%   tolOptions :   [epsilon iterMax]
%    - epsilon :   0 < epsilon < 1; tolerance to stop iteration
%    - iterMax :   integer > 1; force to stop after iterMax iterations
%   (the above arguments inherit from "bcdSpML.m")
%   icType     :   string (default: 'BIC'); set 'AIC' or 'BIC'
%   algType    :   string; "zyue" or "goran"
%
% OUTPUT:
%   lambda     :   scalar, the best that minimises BIC values (if
% 				   multiple, choose the first)
%   Omega      :   PSD matrix, an estimation using the best lambda
%   icVec      :   values of AIC or BIC for lambdaList
%   eTime      :   real vector of time costs for each lambda

% Copyright [2019] <oracleyue>
% Last modified on 24 Jun 2019


% parse arguments
if nargin < 5
    icType = 'BIC';
end
assert(any(strcmp({'AIC', 'BIC'}, icType)), ...
       'Argument "icType" must to "AIC" or "BIC".');
if nargin < 6
    algType = 'zyue';
end
assert(any(strcmp({'zyue', 'goran'}, algType)), ...
       'Argument "algType" must to "zyue" or "goran".');
if strcmp(algType, 'goran')
    addpath('./Goran');
end

% debug flags
debugFlag = 0;

% data
p = length(dL);
d = sum(dL);
numL = length(lambdaList);

% initialize criterion values
icVec = zeros(numL, 1);
estOmega = cell(numL, 1);

% perform estimation and compute criterion values
if debugFlag
    fprintf('Computing KL and BIC values:\n')
end
eTime = zeros(numL, 1);
for k = 1:numL
    algTimer = tic;
    lambda = lambdaList(k);
    switch algType
      case 'zyue'
        [OmegaHat, ~] = bcdSpML(S, dL, lambda, tolOptions);
      case 'goran'
        [~, ~, OmegaHat] = Algorithm(speye(d), S, dL, lambda, ...
                                     tolOptions(2), tolOptions(1), 0);
    end
    estOmega{k} = OmegaHat;

    % refer to /Goran Marjanovic & Victor Solo. ICASSP, 2018/
    switch icType
      case 'BIC'
        icVec(k) = trace(S*OmegaHat) - log(det(OmegaHat)) + ...
            log(N)/N/2 * l0norm(OmegaHat, dL);
      case 'AIC'
        icVec(k) = trace(S*OmegaHat) - log(det(OmegaHat)) + ...
            1/N * l0norm(OmegaHat, dL);
    end

    if debugFlag
        fprintf('    %2d-th iterate: lambda=%.4f, %s=%.4f \n', ...
                k, lambda, icType, icVec(k));
    end
    eTime(k) = toc(algTimer);
end
if debugFlag
    fprintf('End.\n')
    fprintf('Total elapsed time: %f s\n', eTime);
end

% find first lambda that minimises the criterion
[~, idx] = min(icVec);
lambda = lambdaList(idx(1));
Omega = estOmega{idx(1)};

% visualization
if debugFlag
    figure
    semilogx(lambdaList, icVec, 'o-')
    xlabel('Lambda'); ylabel('Information Criterion')
end

end % END of calcLambda


% ================================================================
% Local Functions
% ================================================================
function val = l0norm(Omega, dL)
% PENALTYBIC computes the penalty value used in BIC cirterion for
% block-wise sparse inverse covariance matrix.
%     val = \sum_{i \neq j} I(Omega_{ij} \neq 0) d_i d_j

p = length(dL);
d = sum(dL);

val = 0;
for i = 1:p-1
    for j = i+1:p
        di = dL(i); dj = dL(j);
        iIdx = sum(dL(1:i-1))+1:sum(dL(1:i));
        jIdx = sum(dL(1:j-1))+1:sum(dL(1:j));
        if sum(sum(Omega(iIdx, jIdx)))
            val = val + di*dj;
        end
    end
end
val = val*2;

end % END of l0norm