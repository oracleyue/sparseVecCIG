function [Omega, Sigma] = bcdSpML(S, dL, lambda, tolOptions)
% BCDSPML the block-wise cyclic decent method with conjugate gradient
% embedded for group L0 penalised log-likelihood maximation problems. It
% optimises:
%
%    MAX(Omega) -log det(Omega) + tr(S Omega) + lambda \sum{I(Omega_ij \neq 0)}
%
% INPUT:
%   S          :   (d x d) sample covariance matrix, normalised by T
%   dL         :   (p x 1) vector of positive integers, and Sum(dL) = d
%   lambda     :   positive real number, regularisation parameter
%   tolOptions :   [epsilon iterMax] or epsilon or maxIter
%    - epsilon :   0 < epsilon < 1; tolerance to stop iteration
%    - iterMax :   integer > 1; force to stop after iterMax iterations
%
% OUTPUT:
%   Omega  :   inverse covariance matrix, i.e. inv(Sigma)
%   Sigma  :   covariance matrix of the Gaussian r.v. x(t)

% Copyright (c) 2019, Zuogong YUE
% Author: Zuogong YUE <oracleyue@gmail.com>
%         https://github.com/oracleyue
% Licensed under the GNU General Public License
%
% Last update on 20 Jun 2019


% Flags
debugFlag = 0;

% Preparations
d = size(S, 1);
p = length(dL);
assert(sum(dL) == d, 'Sample covariance S and partition dL fail to match!');
if nargin < 4
    epsilon = 1e-4;  % convergence precision
    iterMax = 50;
end
if nargin == 4
    if isscalar(tolOptions)
        if tolOptions < 1
            epsilon = tolOptions;
            iterMax = 50;
        else
            iterMax = tolOptions;
            epsilon = 10e-4;
        end
    else
        epsilon = tolOptions(1);
        iterMax = tolOptions(2);
    end
end

% Initialization
dLprev = dL;
Sigma = zeros(d,d);
Omega = zeros(d,d);
invSa = cell(p,1);
fvalPrev = d;  % init value of objective function
for k = 1:p
    iIdx = sum(dL(1:k-1)) + 1;
    jIdx = sum(dL(1:k));
    Sk = S(iIdx:jIdx, iIdx:jIdx);
    invSk = inv(Sk);
    invSa{k} = invSk;
    Sigma(iIdx:jIdx, iIdx:jIdx) = Sk;
    Omega(iIdx:jIdx, iIdx:jIdx) = invSk;
    fvalPrev = fvalPrev + log(det(Sk));
end
% pPos relates dLprev to dL, which allows us to rearrange the solution
% according to the original order of matrix blocks.
pPos = 1:p;

if debugFlag
    kDIter = -100:1:0;
    fvalList = zeros(size(kDIter));
    fvalList(end) = fvalPrev;
    fdvalList = zeros(size(kDIter));
    fdvalList(end) = log10(abs(fvalPrev));
    figure(1)
    subplot(1,2,1)
    pltHfval = plot(kDIter, fvalList, '*-');  % plH-: plot handler
    xlabel('Iterations');
    ylabel('Loss function')
    subplot(1,2,2)
    pltHfdval = plot(kDIter, fdvalList, '*-');  % plH-: plot handler
    xlabel('Iterations');
    % ylabel('Log abs. loss difference');
    ylabel('Log rel. loss difference');
end

% Cycle descent over diagonal blocks (i.e. OmegAa)
kIter = 0;  % iteration index
while 1  % cycle in sequence over diagonal block:
         % (p-1)-block to p position, and p-block to 1st position
    kIter = kIter + 1;
    % update dLnext and permutation matrix to update Mo
    % [P, dLnext, pPos] = circPerm(dLprev, pPos);
    dLnext = circshift(dLprev, 1);
    pPos = circshift(pPos, 1);


    % update Omega (being permuted)
    Omega = circshift(Omega, [dLprev(p) dLprev(p)]);
    Ba = Omega(1:end-dLnext(p), end-dLnext(p)+1:end);

    % permute S
    S = circshift(S, [dLprev(p) dLprev(p)]);  % one can use P to
                                              % permute, but this is faster
    So = S(1:end-dLnext(p), 1:end-dLnext(p));
    Soa = S(1:end-dLnext(p), end-dLnext(p)+1:end);
    Sa = S(d-dLnext(p)+1:d, d-dLnext(p)+1:d);
    % retrieve inv(Sa);
    Ta = invSa{pPos(p)};

    % permute Sigma
    Sigma = circshift(Sigma, [dLprev(p) dLprev(p)]);
    SigmAo = Sigma(1:end-dLnext(p), 1:end-dLnext(p));
    Ga = Sigma(1:end-dLnext(p), end-dLnext(p)+1:end);
    SigmAa = Sigma(end-dLnext(p)+1:end, end-dLnext(p)+1:end);

    % update Mo via Mo = SigmAo - Ga*inv(SigmAa)*Ga'
    % reduce inv(OmegAo) (d-dp x d-dp) to smaller inv(SigmAa) (dp x dp)
    Mo = SigmAo - Ga/SigmAa*Ga';

    % Nexted cycle descent loop to update off-diagonal elements
    % update OmegAa, Ba
    dl = dLnext(1:end-1);
    for k = 1:p-1
        % partition of Mo
        Mia = Mo(1:dl(1), 1:dl(1));
        Mmia = Mo(1:dl(1), dl(1)+1:end);
        Mmii = Mo(dl(1)+1:end, dl(1)+1:end);
        % partition of Ba and Soa
        Bia = Ba(1:dl(1), :);
        Bmia = Ba(dl(1)+1:end, :);
        Sioa = Soa(1:dl(1), :);
        Smioa = Soa(dl(1)+1:end, :);

        % use naive inverse
        % BiaStar = - Mia \ (Mmia*Bmia*Sa + Sioa) * Ta;
        % use vector CG method
        BiaStar = qpMatCG(-Mia, (Mmia*Bmia*Sa + Sioa) * Ta);

        % update Bia, i.e. $B_{ia}^+$
        lambdAia = .5 * trace(Sa*BiaStar'*Mia*BiaStar);
        if lambdAia > lambda
            BiaPlus = BiaStar;
        else
            BiaPlus = zeros(size(BiaStar));
        end

        % stacking based on BiaPlus
        Ba = [BiaPlus; Bmia];  % Ba+

        % permute
        Mo = circshift(Mo, [-dLnext(k), -dLnext(k)]);
        Ba = circshift(Ba, -dLnext(k));
        Soa = circshift(Soa, -dLnext(k));
        dl = circshift(dl, -1);
    end

    % update new Omega's and Simga's blocks
    iIdx = d-dLnext(p)+1;
    OmegAo = Omega(1:iIdx-1, 1:iIdx-1);
    OmegAa = Ba'*Mo*Ba + Ta;
    Ka = Mo * Ba;
    Ga = -Ka * Sa;
    SigmAo = Mo + Ka*Sa*Ka';

    % stopping criterion
    % fvalNext = evalObjFunc(Omega, S, lambda, dLnext);
    % dimension-reduced version; faster
    fvalNext = evalLossRe(OmegAo, OmegAa, Ba, Mo, ...
                             So, Sa, Soa, lambda, dLnext);
    % if abs(fvalPrev - fvalNext) < epsilon && kIter > p
    if abs(fvalPrev - fvalNext)/abs(fvalPrev) < epsilon && kIter > p
        if debugFlag
            fprintf('Block cyclic decent stops at the %d-th iteration,\n', ...
                    kIter);
            % fprintf('with difference of loss values: %d\n', ...
            %         abs(fvalNext - fvalPrev));
            fprintf('with relative difference of loss values: %d\n', ...
                    abs(fvalPrev - fvalNext)/abs(fvalPrev));
        end
        break
    end
    if kIter > iterMax * p
        warning(['Convergence is very slow, and ' ...
                 'you may set lambda too small.'])
        break
    end

    % debugging
    if debugFlag
        kDIter = kDIter + 1;
        fvalList = circshift(fvalList, -1);
        fvalList(end) = fvalPrev;
        fdvalList = circshift(fdvalList, -1);
        % fdvalList(end) = log10(abs(fvalNext - fvalPrev));
        fdvalList(end) = log10(abs(fvalPrev - fvalNext)/abs(fvalPrev));
        set(pltHfval, 'Xdata', kDIter, 'Ydata', fvalList);
        set(pltHfdval, 'Xdata', kDIter, 'Ydata', fdvalList);
        pause(.2);
    end

    % update Omega
    Omega(1:iIdx-1, iIdx:d) = Ba;
    Omega(iIdx:d, 1:iIdx-1) = Ba';
    Omega(iIdx:d, iIdx:d) = OmegAa;

    % update Sigma
    Sigma(1:iIdx-1, 1:iIdx-1) = SigmAo;
    Sigma(1:iIdx-1, iIdx:d) = Ga;
    Sigma(iIdx:d, 1:iIdx-1) = Ga';
    Sigma(iIdx:d, iIdx:d) = Sa;

    % backup current fval
    fvalPrev = fvalNext;
    % backup dLnext
    dLprev = dLnext;

end

Omega = retriOrder(Omega, dLnext, pPos);
Sigma = retriOrder(Sigma, dLnext, pPos);

end % END of bcdpML


% ================================================================
% Local Functions
% ================================================================

function [P, dLnext, pPos] = circPerm(dL, pPos)
% Permutation matrix for digonal blocks that moves the k-th block to the
% p-th position (i.e. OmegAa) and pushes the p-th block to the 1st.
% e.g. [1 2 p-1 p] --> [p 1 2 3 p-1]

    d = sum(dL);
    p = length(dL);
    Id = eye(d, d);

    dLnext = circshift(dL, 1);
    pPos = circshift(pPos, 1);

    pBi = d - dL(p) + 1;
    P = [Id(pBi:d,:); Id(1:pBi-1,:)];

end % END of bPerm

function [P, dLnext, pPos] = kPerm(dL, k, pPos)
% Permutation matrix for digonal blocks that moves the k-th block to the
% p-th position (i.e. OmegAa) and pushes the p-th block to the 1st.
% e.g. [1 2 k 3 p] --> [p 1 2 3 k]

    d = sum(dL);
    p = length(dL);
    Id = eye(d, d);

    dLnext = [dL(p) dL(1:k-1) dL(k+1:p-1) dL(k)];
    pPos = [pPos(p) pPos(1:k-1) pPos(k+1:p-1) pPos(k)];

    kBi = sum(dL(1:k-1)) + 1;  % start index of k-th block
    kBj = sum(dL(1:k));        % end index of k-th block
    pBj = d;
    pBi = d - dL(p) + 1;

    P = [Id(pBi:pBj,:); Id(1:kBi-1,:); Id(kBj+1:pBi-1,:); Id(kBi:kBj,:)];

end % END of bPerm

function matNew = retriOrder(mat, dL, pPos)
% Retrieve the shape indicated by the original dL.

    d = sum(dL);
    p = length(dL);
    assert(p == length(pPos), ...
           'Error: the tracking pPos has a different length from dL!');
    idx1 = find(pPos == 1);
    upleft = sum(dL(1:idx1-1));
    matNew = circshift(mat, [-upleft -upleft]);

end % END of retriOrder

function fval = evalObjFunc(Omega, S, lambda, dL)
% Evaluate the cost function.
% Faster version: see "evalLossRe()"

    fval = -logdet(Omega) + trace(S * Omega);

    p = length(dL);
    d = sum(dL);
    zNorm = 0;
    for i = 1:p
        for j = 1:p
            iIdxU = sum(dL(1:i-1)) + 1;
            iIdxD = sum(dL(1:i));
            jIdxL = sum(dL(1:j-1)) + 1;
            jIdxR = sum(dL(1:j));
            if (i ~= j) && sum(sum(Omega(iIdxU:iIdxD, jIdxL:jIdxR)))
                zNorm = zNorm + 1;
            end
        end
    end

    fval = fval + lambda * zNorm;

end % END of evalObjFunc

function fval =  evalLossRe(OmegAo, OmegAa, Ba, invOmegAo, ...
                               So, Sa, Soa, lambda, dL)
% Evaluate the cost function via its dimensionally reducted form.

    fval = -logdet(OmegAo) - logdet(OmegAa - Ba'*invOmegAo*Ba) + ...
           trace(So*OmegAo) + 2*trace(Soa'*Ba) + trace(Sa*OmegAa);

    p = length(dL);
    d = sum(dL);
    zNormOmega = 0;
    zNormBa = 0;
    for i = 1:p-1
        for j = 1:p-1
            iIdxU = sum(dL(1:i-1)) + 1;
            iIdxD = sum(dL(1:i));
            jIdxL = sum(dL(1:j-1)) + 1;
            jIdxR = sum(dL(1:j));
            if (i ~= j) && sum(sum(OmegAo(iIdxU:iIdxD, jIdxL:jIdxR)))
                zNormOmega = zNormOmega + 1;
            end
        end
        idxU = sum(dL(1:i-1)) + 1;
        idxD = sum(dL(1:i));
        if sum(sum(Ba(idxU:idxD, :)))
            zNormBa = zNormBa + 1;
        end
    end

    fval = fval + lambda * (zNormOmega + 2*zNormBa);

end % END of evalObjFunc

function val = logdet(X)
% Reliably compute log(det(X)), where X must be positive definte matrix.

    L = chol(X);
    val = 2*sum(log(diag(L)));

end % END of logdet

function X = qpMatCG(M, W, epsilon)
% Conjugate gradient method (matrix version) for the matrix quadratic
% programming:
%    MAX{X} tr(X'MX) + 2 tr(X'W)
% equivalently, solving MX = -W, where M are sysmetric positivie definite.
% Moreover, in our case which maximizes tr(SX'MX) + 2 tr(X'W), it is
% equivalent to solve MXS = -W, where S is symmetric and positive definite.

    [m, m1] = size(M);
    [m2, n] = size(W);
    assert((m == m1) & (m == m2), ...
           'The dimensions of M and W fail to match!')
    if nargin < 3
        epsilon = 1e-20;
    end

    % Initialization
    R = W;
    D = R;
    X = zeros(m, n);
    rho = trace(R'*R);

    % CG update
    for k = 0:1:m
        % update X
        alpha = rho / trace(D'*M*D);
        % stop if possible
        % [Warning]: should guarantee the decrease of loss function
        if k == 0
            alpha0rho0 = alpha*rho;
        elseif abs(alpha*rho/alpha0rho0) < epsilon
            break
        end
        X = X + alpha*D;  % X_k+1

        % update variables for (k+1) iterate
        R = R - alpha*M*D;
        rho_p = trace(R'*R);
        gamma = rho_p / rho;
        D = R + gamma*D;
        rho = rho_p;
    end

end % END of qpMatCG