function [Omega, Sigma, optStatus] = bcdSpML(S, dL, lambda, tolOptions)
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
%   optStatus  : structure; information about optimization processes

% Copyright (c) 2019, Zuogong YUE
% Author: Zuogong YUE <oracleyue@gmail.com>
%         https://github.com/oracleyue
% Licensed under the GNU General Public License

% Log of updates:
%
% - Adding several stopping conditions due to observations of different
% converence cases.
%
% - Update the way to compute log(det(X)), by noticing that MATLAB
% default det() cannot handle large matrices, which gives +/-Inf and
% fails our algorithm. By taking advantages of positive definitiness of
% X, we use chol() to compute log(det(X)).
%
% - Replace "circshift" when updating Bia, since "circshift" is pretty
% costly when dealing with large matrices. Now we directly retrieve and
% update Bia by directly element positioning.
%
% - Matrix CG method may not be necessary when submatrix is not large
% enough, e.g. size < 100x100. Here we add a threshold "maxDiagMatrix"
% to choose between Matlab default "mldivide" and our CG method.
%
% - By using MATLAB "profile", we observe that the major cost of time is
% not the algorithm body, but "evalLoss()" (computing the value of
% objective function) and children funtion "l0norm()". Thus, we update it
% to recursively compute ||Omega||_0, i.e. "zNorm", instead of calling
% "l0norm()" separately. Another time consuming operation in
% "evalLoss()" is "trace(S*Omega)", which is pending to be updated.

% Last update on 04 Jul 2019


% Flags
debugFlag = 0;

% Preparations
d = size(S, 1);
p = length(dL);
maxDiagMatrix = max(dL);  % maximal size of submatrix on diagonal
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
OmTemp = zeros(d,d);   % temporary Omega
SigTemp = zeros(d,d);  % temporary Sigma
invSa = cell(p,1);
fval = d;  % init value of objective function
zNorm = 0; % block-wise zero norm of Omega (no diagonal)
zNormPm = zNorm; % block-wise zNorm of Omega when rem(kIter, p) == 0
zNormFull = p^2 - p; % block-wise zero norm of full Omega matrix
for k = 1:p
    iIdx = sum(dL(1:k-1)) + 1;
    jIdx = sum(dL(1:k));
    Sk = S(iIdx:jIdx, iIdx:jIdx);
    invSk = inv(Sk);
    invSa{k} = invSk;
    Sigma(iIdx:jIdx, iIdx:jIdx) = Sk;
    Omega(iIdx:jIdx, iIdx:jIdx) = invSk;
    fval = fval + log(det(Sk));
end
% pPos relates dLprev to dL, which allows us to rearrange the solution
% according to the original order of matrix blocks.
pPos = 1:p;

if debugFlag
    kDIter = -100:1:0;
    fvalList = zeros(size(kDIter));
    fvalList(end) = fval;
    fdvalList = zeros(size(kDIter));
    fdvalList(end) = log10(abs(fval));
    figure(1)
    subplot(1,2,1)
    pltHfval = plot(kDIter, fvalList, '*-');  % plH-: plot handler
    xlabel('Iterations');
    ylabel('Loss function')
    subplot(1,2,2)
    pltHfdval = plot(kDIter, fdvalList, '*-');  % plH-: plot handler
    xlabel('Iterations');
    ylabel('Incremental loss value (log10)');
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
    % Mo = SigmAo - Ga/SigmAa*Ga';
    % use Cholesky to accelerate
    R = chol(SigmAa);
    GRinv = Ga / R;
    Mo = SigmAo - GRinv*GRinv';

    % Nexted cycle descent loop to update off-diagonal elements
    % update OmegAa, Ba
    dl = dLnext(1:end-1);
    for k = 1:p-1
        % calculate block indices (to avoid permutation via circshift)
        ibIdx = sum(dl(1:k-1))+1 : sum(dl(1:k));
        jbIdxPre  = 1:sum(dl(1:k-1));
        jbIdxPost = sum(dl(1:k))+1 : sum(dl);

        % partition of Mo
        Mia = Mo(ibIdx, ibIdx);
        Mmia = [Mo(ibIdx, jbIdxPost) Mo(ibIdx, jbIdxPre)];
        % Mia = Mo(1:dl(1), 1:dl(1));
        % Mmia = Mo(1:dl(1), dl(1)+1:end);

        % partition of Ba and Soa
        Bia = Ba(ibIdx, :);
        Bmia = [Ba(jbIdxPost,:); Ba(jbIdxPre,:)];
        % Bmia = Ba(dl(1)+1:end, :);
        Sioa = Soa(ibIdx,:);
        % Sioa = Soa(1:dl(1), :);

        % compute Bia*
        if maxDiagMatrix < 100  % use naive inverse
            BiaStar = - Mia \ (Mmia*Bmia*Sa + Sioa) * Ta;
        else  % use matrix CG method
            BiaStar = qpMatCG(-Mia, (Mmia*Bmia*Sa + Sioa) * Ta);
        end

        % update Bia+, i.e. $B_{ia}^+$
        lambdAia = .5 * trace(Sa*BiaStar'*Mia*BiaStar);
        if lambdAia > lambda
            BiaPlus = BiaStar;
        else
            BiaPlus = zeros(size(BiaStar));
        end

        % update (block) zNorm recursively, since
        % l0norm() in evalLoss() is costly for large matrices
        isNonZeroBlk = any(Bia, 'all');
        if isNonZeroBlk && lambdAia <= lambda
            zNorm = zNorm - 2;
        elseif ~isNonZeroBlk && lambdAia > lambda
            zNorm = zNorm + 2;
        end

        % update Ba with Bia+
        Ba(ibIdx,:) = BiaPlus;
        % Ba = [BiaPlus; Bmia];

        % % permute
        % Mo = circshift(Mo, [-dLnext(k), -dLnext(k)]);
        % Ba = circshift(Ba, -dLnext(k));
        % Soa = circshift(Soa, -dLnext(k));
        % dl = circshift(dl, -1);
    end

    % update new Omega's and Simga's blocks
    iIdx = d-dLnext(p)+1;
    OmegAo = Omega(1:iIdx-1, 1:iIdx-1);
    OmegAa = Ba'*Mo*Ba + Ta;
    Ka = Mo * Ba;
    Ga = -Ka * Sa;
    SigmAo = Mo + Ka*Sa*Ka';

    % setup Omega
    OmTemp(1:iIdx-1, 1:iIdx-1) = Omega(1:iIdx-1, 1:iIdx-1);
    OmTemp(1:iIdx-1, iIdx:d) = Ba;
    OmTemp(iIdx:d, 1:iIdx-1) = Ba';
    OmTemp(iIdx:d, iIdx:d) = OmegAa;

    % setup Sigma
    SigTemp(1:iIdx-1, 1:iIdx-1) = SigmAo;
    SigTemp(1:iIdx-1, iIdx:d) = Ga;
    SigTemp(iIdx:d, 1:iIdx-1) = Ga';
    SigTemp(iIdx:d, iIdx:d) = Sa;

    % Stopping Criteria
    stopCase = 0;
    tolType = 'rel';  % "abs" "rel"
    [fvalNext, errsig] = evalLoss(OmTemp, S, lambda*zNorm);
    % dimensionally reduced version
    % [fvalNext, zNormNext] = evalLossRe(OmegAo, OmegAa, Ba, Mo, ...
    %                                    So, Sa, Soa, lambda, dLnext);
    % use abs./rel. error in convergence checking
    switch tolType
      case {'absolute', 'abs', 'abs.'}
        tolType = 'absolute';
        fvalInc = abs(fval - fvalNext);
      case {'relative', 'rel', 'rel.'}
        tolType = 'relative';
        fvalInc = abs(fval - fvalNext)/abs(fval);
    end
    % case I: successful convergence
    if kIter > p ... % at least iterate over all diagonal lements
        && fval >= fvalNext ...  % fval is decreasing
        && fvalInc < epsilon
        if debugFlag
            fprintf('Block cyclic decent stops at the %d-th iteration,\n', ...
                    kIter);
            fprintf('with %s decrease of loss values: %d\n', ...
                    tolType, fvalInc);
        end
        stopCase = 1;
        stopMsg = 'Succeed';
    end
    % case II: run out of maximal iterations
    if kIter > iterMax * p
        warning(['Convergence is very slow, and ' ...
                 'you may set lambda too small.'])
        stopCase = 2;
        stopMsg = 'Run out of maximal iterations';
    end
    % case III: non-sparsity acquired due to too small lambdas
    if zNormPm == zNormFull && zNorm == zNormFull ...
            && ~rem(kIter, p)
        stopCase = 3;
        stopMsg = 'Stop due to non-sparsity; return inv(S) as MLE(Omega)';
    end
    % case IV: divergence
    if fvalNext > fval
        % stopCase = 4;
        stopMsg = 'Divergence observed';
    end
    % case V: fail in evaluating loss functions (like chol() or det())
    if errsig
        stopCase = 5;
        stopMsg = 'Failure in evaluating likelihood';
    end
    % debugging
    if debugFlag
        kDIter = kDIter + 1;
        fvalList = circshift(fvalList, -1);
        fvalList(end) = fvalNext;
        fdvalList = circshift(fdvalList, -1);
        fdvalList(end) = log10(fvalInc);
        set(pltHfval, 'Xdata', kDIter, 'Ydata', fvalList);
        set(pltHfdval, 'Xdata', kDIter, 'Ydata', fdvalList);
        pause(.1);
    end
    % report optimization status
    if stopCase
        optStatus.Iterations = kIter;
        optStatus.OptimalValue = fval;
        optStatus.LastIncrement = -fvalInc;
        optStatus.ToleranceType = [tolType ' error'];
        optStatus.l0Norm = l0norm(Omega, dLprev, 'element');
        optStatus.Groupl0Norm = zNorm;
        optStatus.Status = stopMsg;
    end
    % handle stopping cases
    switch stopCase
      case {1, 2}
        Omega = OmTemp;
        Sigma = SigTemp;
        break
      case 3
        Omega = inv(S);
        Sigma = S;
        return
      case 4
        optStatus.LastIncrement = fvalInc;
        break
      case 5
        switch errsig
          case 1
            warning('"chol()" failed in evaluating likelihood!')
          case 2
            warning('"det()" failed in evaluating likelihood!')
        end
        break
    end

    % proceed to the next iteration
    Omega = OmTemp;
    Sigma = SigTemp;
    % update fval
    fval = fvalNext;
    % save zNorm when iterating at multiplicity of p
    if ~rem(kIter, p)
        zNormPm = zNorm;
    end
    % update dLnext
    dLprev = dLnext;

end

% restore original shape of Omega and Sigma
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

% function [fval, zNorm, errsig] = evalLoss(Omega, S, lambda, dL)
% % Evaluate the cost function.
% % Faster version: see "evalLossRe()"

%     [lval, errsig] = logdet(Omega, 'chol');
%     fval = -lval + trace(S * Omega);

%     p = length(dL);
%     d = sum(dL);
%     zNorm = l0norm(Omega, dL, 'block');
%     fval = fval + lambda * zNorm;

% end % END of evalObjFunc

function [fval, errsig] = evalLoss(Omega, S, lambNorm)
% Evaluate the cost function.

    [lval, errsig] = logdet(Omega, 'chol');
    fval = -lval + trace(S * Omega) + lambNorm;

end % END of evalObjFunc

function [fval, zNorm] =  evalLossRe(OmegAo, OmegAa, Ba, invOmegAo, ...
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

    zNorm = zNormOmega + 2*zNormBa;
    fval = fval + lambda * zNorm;

end % END of evalObjFunc

function [val, errsig] = logdet(X, method)
% Reliably compute log(det(X)), where X must be positive definte matrix.
% Notes:
%   errsig = 0 : no errors
%   errsig = 1 : chol() fails
%   errsig = 2 : det() fails

    if nargin < 2
        method = 'det';
    end
    assert(any(strcmpi({'eig', 'chol', 'det'}, method)), ...
           'The argument "method" must be "eig", "chol" or "det".');
    errsig = 0;  % no errors

    switch method
      case 'chol'
        % MATLAB 2019a has a severe bug on "chol"!
        % It fails on matrices with minimal eigenvalue larger than 5e-4!
        [L, flag] = chol(X);
        if flag
            errsig = 1;  % chol on a non-positive definite matrix
            msg = sprintf(['An approximation of log(det(Omega)) is computed'...
                           ', due to %d number of eigenvalues '...
                           'that are particularly close to 0.'], ...
                          size(X,1)-flag+1);
            warning(msg);
        end
        val = 2*sum(log(diag(L)));
      case 'eig'
        val = log(prod(eig(X)));
      case 'det'
        % This method is not reliable when Omega is large, e.g. dim > 400.
        val = log(det(X));
        if isinf(val)
            errsig = 2;  % det gives Inf value
        end
    end

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

function val = l0norm(Omega, dL, type)
% PENALTYBIC computes the penalty value used in BIC cirterion for
% element/block-wise sparse inverse covariance matrix.
% element-wise:
%     val = \sum_{i \neq j} I(Omega_{ij} \neq 0) d_i d_j
% block-wise:
%     val = \sum_{i \neq j} I(Omega_{ij} \neq 0)
% INPUTS:
% 	  type  :  string; 'element' or 'block'

if nargin < 3
    type = 'element';
end
assert(any(strcmp({'element', 'block'}, type)), ...
       'The argument "type" has to be either "element" or "block".');

p = length(dL);
d = sum(dL);

val = 0;
for i = 1:p-1
    for j = i+1:p
        di = dL(i); dj = dL(j);
        iIdx = sum(dL(1:i-1))+1:sum(dL(1:i));
        jIdx = sum(dL(1:j-1))+1:sum(dL(1:j));
        if any(Omega(iIdx, jIdx), 'all')
            switch type
              case 'element'
                val = val + di*dj;
              case 'block'
                val = val + 1;
            end
        end
    end
end
val = val*2;

end % END of l0norm
