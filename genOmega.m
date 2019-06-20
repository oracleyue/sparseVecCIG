function Omega = genOmega(dL)

d = sum(dL);
Omega = zeros(d,d);

for k=1:length(dL)

    U = sprand(dL(k),dL(k),0.3)*eye(dL(k));
    U(U~=0&U>0.5) = 1;
    U(U~=0&U<=0.5) = -1;

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

i = 2; j = 4;
Om_ij = sprand(dL(i),dL(j),0.8)*eye(dL(j));
Om_ij(Om_ij~=0&Om_ij>0.5) = 1;
Om_ij(Om_ij~=0&Om_ij<=0.5) = -1;

Omega(sum(dL(1:i-1))+1:sum(dL(1:i)),...
   sum(dL(1:j-1))+1:sum(dL(1:j))) = Om_ij;

Omega(sum(dL(1:j-1))+1:sum(dL(1:j)),...
   sum(dL(1:i-1))+1:sum(dL(1:i))) = Om_ij';

i = 3; j = 1;
Om_ij = sprand(dL(i),dL(j),0.2)*eye(dL(j));
Om_ij(Om_ij~=0&Om_ij>0.5) = 1;
Om_ij(Om_ij~=0&Om_ij<=0.5) = -1;

Omega(sum(dL(1:i-1))+1:sum(dL(1:i)),1:dL(j)) = Om_ij;
Omega(1:dL(j),sum(dL(1:i-1))+1:sum(dL(1:i))) = Om_ij';

% extra:
i = 5; j = 1;
Om_ij = sprand(dL(i),dL(j),0.2)*eye(dL(j));
Om_ij(Om_ij~=0&Om_ij>0.5) = 1;
Om_ij(Om_ij~=0&Om_ij<=0.5) = -1;

Omega(sum(dL(1:i-1))+1:sum(dL(1:i)),1:dL(j)) = Om_ij;
Omega(1:dL(j),sum(dL(1:i-1))+1:sum(dL(1:i))) = Om_ij';

while min(eig(Omega))<=0
    Omega = Omega + eye(d);
end

end