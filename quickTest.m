% [OmegaHat, SigmaHat] = bcdpML(S, dL, lambda, 1e-3);

m = 6; n = 1;
M = symdec(m, 1.8);
X = rand(m, n);
W = M*X;

Xhat = qpMatCG(M, W);

norm(W - M*(M\W), 2)
norm(W - M*Xhat, 2)