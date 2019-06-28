function val = l0norm(Omega, dL, type)
% PENALTYBIC computes the penalty value used in BIC cirterion for
% block-wise sparse inverse covariance matrix.
%     val = \sum_{i<j} I(Omega_{ij} \neq 0) d_i d_j

% Copyright [2019] <oracleyue>
% Last modified on 24 Jun 2019


if nargin < 3
    type = 'element';
end
if ~any(strcmp({'element', 'block', 'group'}, type))
    error('The argument "type" is wrong, which must be: element, block, group.')
end

p = length(dL);
d = sum(dL);

val = 0;
for i = 1:p-1
    for j = i+1:p
        di = dL(i); dj = dL(j);
        iIdx = sum(dL(1:i-1))+1:sum(dL(1:i));
        jIdx = sum(dL(1:j-1))+1:sum(dL(1:j));
        if sum(sum(Omega(iIdx, jIdx)))
            switch type
              case 'element'
                val = val + di*dj;
              case {'block', 'group'}
                val = val + 1;
            end
        end
    end
end
val = val*2;

end % END of l0norm
