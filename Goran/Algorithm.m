
% SPARSE BLOCK INVERSE COVARIANCE MATRIX ALGORITHM

% INPUTS:
%--------
% K0 = initial DIAGONAL Inverse Covariance matrix 
%      (in sparse format). 
% S = sample Covariance matrix.
% d_1_p = vector of dimensions d_1,...,d_p.
%     Note that sum(d) is the dimension of K.
% lam = penalty parameter.
% N = number of iterations of the Algorithm 
%     (complete sweeps (through all blocks) of K).
% tol = tolerance for early termination.
% pen = the subalgorithm used:
%       pen == 0 is group l0, pen == 1 is group l1.

% OUTPUTS:
%---------
% F = objective function.
% The function looks like:
% F = -log(det(K)) + 0.5*tr(S*K) + lam*sum_{i~=j}p(K_ij), 
% where
% p(K_ij) = |K_ij|_0 if pen == 0
% p(K_ij) = |K_ij|_F if pen == 1.
% Time is the running time of the algorithm
% as a function of iteration.
% K = final Inverse Covariance matrix.

function [F,Time,K] = Algorithm(K0,S,d_1_p,lam,N,tol,pen)

time = cputime;
diagK = diag(K0);
log_detK = sum(log(diagK));
diagS = diag(S);
trSK = diagS'*diagK;

K = K0;
d = sum(d_1_p);
time = cputime-time;

F = zeros(1,N);
Time = F;
if pen == 0
    F(1) = -log_detK + 0.5*trSK + lam*group_l0(K,d_1_p);
elseif pen == 1
    F(1) = -log_detK + 0.5*trSK + lam*group_l1(K,d_1_p);
end
    
Time(1) = time;

clear diagS K0 diagK

for k=2:N
    
    time = cputime;
    
    % inner loop: cycling through d_1,...,d_p
    for i=1:length(d_1_p)
        
        if i==1
            Ka = K(1:d_1_p(1),1:d_1_p(1));
            Ko = K(d_1_p(1)+1:d,d_1_p(1)+1:d);
            Ba = K(d_1_p(1)+1:d,1:d_1_p(1));
            
            Sa = S(1:d_1_p(1),1:d_1_p(1));
            Soa = S(d_1_p(1)+1:d,1:d_1_p(1));
            
        elseif i>1 && i < d
            sum_d_i_1 = sum(d_1_p(1:i-1));
            Ka = K(sum_d_i_1+1:sum_d_i_1+d_1_p(i),...
                   sum_d_i_1+1:sum_d_i_1+d_1_p(i));
            % constructing Ko:
            % put together 3 matrices
            % Ko11, Ko12, Ko12', Ko22
            Ko11 = K(1:sum_d_i_1,1:sum_d_i_1);
            Ko12 = K(1:sum_d_i_1,sum_d_i_1+d_1_p(i)+1:d);
            Ko22 = K(sum_d_i_1+d_1_p(i)+1:d,...
                     sum_d_i_1+d_1_p(i)+1:d);
            Ko = [Ko11,Ko12;Ko12',Ko22];
            clear Ko11 Ko12 Ko22
            
            % constructing Ba:
            % put together 2 matrices
            % B_1a, B_2a
            B_1a = K(1:sum_d_i_1,sum_d_i_1+1:sum_d_i_1+d_1_p(i));
            B_2a = K(sum_d_i_1+d_1_p(i)+1:d,...
                     sum_d_i_1+1:sum_d_i_1+d_1_p(i));
            Ba = [B_1a;B_2a];
            clear B_1a B_2a
            
            Sa = S(sum_d_i_1+1:sum_d_i_1+d_1_p(i),...
                   sum_d_i_1+1:sum_d_i_1+d_1_p(i));
            
            % constructing Soa:
            % put together 3 matrices
            % So_1a, So_2a
            So_1a = S(1:sum_d_i_1,sum_d_i_1+1:sum_d_i_1+d_1_p(i));
            So_2a = S(sum_d_i_1+d_1_p(i)+1:d,...
                      sum_d_i_1+1:sum_d_i_1+d_1_p(i));
            Soa = [So_1a;So_2a];
            clear So_1a So_2a
            
        else % when d_1_p(i) = d_p (the last d_i)
            Ka = K(d-d_1_p(i)+1:d,d-d_1_p(i)+1:d);
            Ko = K(1:d-d_1_p(i),1:d-d_1_p(i));
            Ba = K(1:d-d_1_p(i),d-d_1_p(i)+1:d);
            
            Sa = S(d-d_1_p(i)+1:d,d-d_1_p(i)+1:d);
            Soa = S(1:d-d_1_p(i),d-d_1_p(i)+1:d);
        end
        
        KoI = Ko\speye(d-d_1_p(i));
        % fixing numerical error:
        KoI = 0.5*(KoI + KoI);
        
        % now we need to update Ba:
        % for this we use a CD procedure.
        % letting beta_a = vec(Ba) = Ba(:), and
        % MoI = kron(Sa,KoI);
        % soa = Soa(:);
        % the objective function is
        % J = beta_a'*MoI*beta_a + 2*beta_a'*soa + 2*lam*sum_i p(beta_ia),
        % where
        % p(beta_ia) = |beta_ia|_0 if pen == 0
        % p(beta_ia) = |beta_ia|_F if pen == 1
        
        L = chol(kron(Sa,KoI));
        index_matrix = sparse(get_index_matrix(i,d_1_p));
        ba_new = group_PLS(L,Ba(:),Soa(:),index_matrix,lam,1e2,1e-2,pen);
        
        % update Ba:
        Ba_new = reshape(ba_new,[d-d_1_p(i),d_1_p(i)]);
        clear ba_new;
        % update Ka:
        SaI = Sa\speye(d_1_p(i));
        % fixing numerical error:
        SaI = 0.5*(SaI + SaI');
        Ka_new = (Ba_new'*KoI*Ba_new) + SaI;
        clear SaI;
        
        % fixing numerical error:
        Ka_new = 0.5*(Ka_new' + Ka_new);
            
        log_detK = log_detK - log(det(Sa)) ... 
                            - log(det(Ka-(Ba'*KoI*Ba)));
        clear KoI;               
        trSK = trSK + 2*trace(Soa'*(Ba_new-Ba)) ...
                    + trace(Sa*(Ka_new-Ka));
        clear Soa Sa Ba;
                
        % replace Ba with Ba_new and 
        % Ka with Ka_new in matrix K:
        if i==1
            K(1:d_1_p(1),1:d_1_p(1)) = Ka_new;
            
            K(d_1_p(1)+1:d,1:d_1_p(1)) = Ba_new;
            
            K(1:d_1_p(1),d_1_p(1)+1:d) = Ba_new';
            
        elseif i>1 && i<d
            sum_d_i_1 = sum(d_1_p(1:i-1));
            K(sum_d_i_1+1:sum_d_i_1+d_1_p(i),...
              sum_d_i_1+1:sum_d_i_1+d_1_p(i)) = Ka_new;
          
            K(1:sum_d_i_1,sum_d_i_1+1:sum_d_i_1+d_1_p(i)) = ...
            Ba_new(1:sum_d_i_1,:);
            K(sum_d_i_1+d_1_p(i)+1:d,sum_d_i_1+1:sum_d_i_1+d_1_p(i)) = ...
            Ba_new(sum_d_i_1+1:d-d_1_p(i),:);
        
            K(sum_d_i_1+1:sum_d_i_1+d_1_p(i),1:sum_d_i_1) = ...
            Ba_new(1:sum_d_i_1,:)';
            K(sum_d_i_1+1:sum_d_i_1+d_1_p(i),sum_d_i_1+d_1_p(i)+1:d) = ...
            Ba_new(sum_d_i_1+1:d-d_1_p(i),:)';
        
        else
            K(d-d_1_p(i)+1:d,d-d_1_p(i)+1:d) = Ka_new;
            
            K(1:d-d_1_p(i),d-d_1_p(i)+1:d) = Ba_new;
            
            K(d-d_1_p(i)+1:d,1:d-d_1_p(i)) = Ba_new';
        end
        
        clear K_new Ba_new;
        
    end
    
    if pen == 0
        F(k) = -log_detK + trSK + lam*group_l0(K,d_1_p);
    elseif pen == 1
        F(k) = -log_detK + trSK + lam*group_l1(K,d_1_p);
    end
    
    time = cputime - time;
    
    Time(k) = Time(k-1) + time;
    
    % exit condition:
    if abs(F(k-1) - F(k))/abs(F(k-1))<=tol
        F = F(1:k);
        Time = Time(1:k);
        break;
    end
    
end

end