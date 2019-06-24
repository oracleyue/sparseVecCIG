
function beta_a = group_PLS(L,beta_a_0,sig_oa,index_matrix,lam,Iter,tol,pen)

beta_a = beta_a_0;

for k=2:Iter
    
    [~,num_blocks] = size(index_matrix);
    
    v1 = beta_a;
    
    for i=1:num_blocks
       
        index_ia = index_matrix(:,i);
        
        beta_mia = beta_a(index_ia==0);
        
        L_ia = L(:,index_ia==1);
        L_mia = L(:,index_ia==0);
        
        z_ia = (L_ia'*L_mia)*beta_mia + sig_oa(index_ia==1);
        
        beta_ias = -(L_ia'*L_ia)\z_ia;
        
        if pen == 0
        
            if (z_ia'*beta_ias) + 2*lam < 0
                beta_ia_p = beta_ias;
            else
                beta_ia_p = zeros(sum(index_ia),1);
            end
        elseif pen == 1
            
            LtL = L_ia'*L_ia;
            
            beta_ias = -LtL\z_ia;
            for r=1:3
                beta_ias = -LtL\(z_ia+lam*beta_ias/norm(beta_ias,2));
            end
            
            val = (z_ia'*beta_ias)/norm(beta_ias,2)+lam;
            
            if val<0
                beta_ia_p = beta_ias;
            else
                beta_ia_p = zeros(sum(index_ia),1);
            end
            
            beta_a(index_ia==1) = beta_ia_p;
            
        end
        
        beta_a(index_ia==1) = beta_ia_p;
        
    end
    
    v2 = beta_a;
    
    if norm(v1-v2,2)/norm(v1)<=tol
        break;
    end
    
end
end



