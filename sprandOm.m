function Omega = sprandOm(dL, density, plotEnable)
% SPRANDOM generates a block-wise sparse inverse covariance matrix
% (i.e. precision matrix). The partition is specified by a vector `dL`.

% INPUT:
%   dL      :   vector of positive integers;
%               sum(dL) is the dimension of Omega.
%   density :   scalar or 2-dim vector with elements in interval (0,1).
%               The density of sparsity of non-zero blocks and each
%               blocks. If scalar, the same value for both; if vector,
%               the 1st for block density and the 2nd for submatrix sparsity.
%   plotEnable : Boolean (default: 0); set 1 to plot Omega as image.
%
% OUTPUT:
%   Omega   :   precision matrix with dim (sum(dL) x sum(dL))
%
% Examples:
%   Omega = sprandOm([3 2 5 4 2], .4)
%   Omega = sprandOm([3 2 5 4 2], [.4 .8])

% Copyright (c) 2015-2019, Zuogong YUE
% Author: Zuogong YUE <oracleyue@gmail.com>
%         https://github.com/oracleyue
% Licensed under the GNU General Public License
%
% Last update on 22 Jun 2019


d = sum(dL);     % dim of Omega
n = length(dL);  % num of blocks on diagonal, i.e. in total n^2
Omega = zeros(d,d);

% Parsing arguments
if nargin < 1
    error('At least one argument is required!')
elseif nargin == 1
    blkDensity = (n+2)/n^2;
    density = [blkDensity .8];
end
if nargin < 3
    plotEnable = 0;
end
if isscalar(density)
    blkDensity = density;
    subBlkDensity = density;
else
    blkDensity = density(1);
    subBlkDensity = density(2);
end
checkDensity = round(n^2*blkDensity) > n;
msg = ['The density is to low to make Omega be positive definite!' ...
       '\n' 'Notes: the density should specify #(nonzero blocks) being ' ...
       'larger than #(diagonal blocks).'];
assert(checkDensity, msg);

% Generate a block diagonal matrix
for k=1:length(dL)
    U = sprand(dL(k),dL(k), subBlkDensity)*eye(dL(k));
    U(U~=0 & U>0.5) = 1;
    U(U~=0 & U<=0.5) = -1;

    Om_kk = U'*U;
    while min(eig(Om_kk))<=0
        Om_kk = Om_kk + eye(dL(k));
    end

    if k==1
        Omega(1:dL(k),1:dL(k)) = Om_kk;
    else
        Omega(sum(dL(1:k-1))+1:sum(dL(1:k)),...
           sum(dL(1:k-1))+1:sum(dL(1:k))) = Om_kk;
    end
end

% Randomly choose nonzero off-diagonal blocks
numNZblk = round((round(n^2*blkDensity) - n)/2);
iter = 0;
idxList = zeros(numNZblk, 2);
% % This method is too slow when n is particularly large
% while iter < numNZblk
%     while 1
%         iIdx = randi(n);
%         jIdx = randi([iIdx n]);
%         if iIdx < jIdx && isempty(intersect(idxList, [iIdx jIdx], 'rows'))
%             idxList(iter+1, :) = [iIdx jIdx];
%             iter = iter + 1;
%             break
%         end
%     end
% end
% alternative way: fast!
[iVec, jVec] = flatUpperMatrix(n);
selected = randperm(length(iVec), numNZblk);
idxList = [iVec(selected)' jVec(selected)'];

% Set nonzero off-diagonal blocks
for k = 1:numNZblk
    i = idxList(k, 1);
    j = idxList(k, 2);

    Om_ij = sprand(dL(i),dL(j), subBlkDensity)*eye(dL(j));
    Om_ij(Om_ij~=0 & Om_ij>0.5) = 1;
    Om_ij(Om_ij~=0 & Om_ij<=0.5) = -1;

    iBegin = sum(dL(1:i-1))+1; iEnd = sum(dL(1:i));
    jBegin = sum(dL(1:j-1))+1; jEnd = sum(dL(1:j));
    Omega(iBegin:iEnd, jBegin:jEnd) = Om_ij;
    Omega(jBegin:jEnd, iBegin:iEnd) = Om_ij';
end

% Force positiveness
minEigOm = min(eig(Omega));
while minEigOm <= 0
    Omega = Omega + (abs(minEigOm) + 1e-4)*eye(d);
    minEigOm = min(eig(Omega));
end

% Return in sparse data structure
Omega = sparse(Omega);

% Plot if required
if plotEnable
    imshowOm(Omega);
end

end  % END of sprandOm

% ================================================================
% Local Functions
% ================================================================

function [xVec, yVec] = flatUpperMatrix(N)
% It builds up two vectors, which [x,y] is an index specify an element
% in the upper triangle part of a square matrix. For example,
%    N = 4, then xVec = [1, 1, 1, 2, 2, 3],
%                yVec = [2, 3, 4, 3, 4, 4].

    xVec = zeros(1, N*(N-1)/2);
    yVec = zeros(1, N*(N-1)/2);
    pos = 1;
    for k = 1:N
        xVec(pos:pos-1+(N-k)) = k;
        yVec(pos:pos-1+(N-k)) = k+1:N;
        pos = pos+(N-k);
    end

end  % END of flatUpperMatrix