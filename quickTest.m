%% Test algorithms for ML
% [OmegaHat, SigmaHat] = bcdpML(S, dL, lambda, 1e-3);


%% Test CG local function
% m = 6; n = 1;
% M = symdec(m, 1.8);
% X = rand(m, n);
% W = M*X;
% Xhat = qpMatCG(M, W);
% norm(W - M*(M\W), 2)
% norm(W - M*Xhat, 2)


%% Test random Omega generation
% dL = [2, 3, 5, 1, 3, 4, 2, 1]*3;
% % Omega = sparse(genOmega(dL));
% Omega = sprandOm(dL, [.3, .8]);
% imshowOm(Omega);
