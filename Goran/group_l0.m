% group l0 of matrix K

% inputs:
%--------
% K     = symmetric matrix (inverse covariance matrix)
%         in sparse format.
% d_1_p = vector of dimensions d_1,...,d_p.
%         Note, there are p^2 blocks all up.

% outputs:
%---------
% sum_{i~=j}|K_ij|_0.

function penalty = group_l0(K,d_1_p)

p = length(d_1_p);

penalty = 0;
for i=1:p
    for j=i:p
        if j~=i
            if i==1
                lim_i_st = 1;
            else
                lim_i_st = sum(d_1_p(1:i-1))+1;
            end
            lim_i_fin = sum(d_1_p(1:i));
            
            lim_j_st = sum(d_1_p(1:j-1))+1;
            lim_j_fin = sum(d_1_p(1:j));
            
            K_ij = K(lim_i_st:lim_i_fin,lim_j_st:lim_j_fin);
            if sum(K_ij(:) ~= 0)~=0
                penalty = penalty + 2;
            end
        end
    end
end

end