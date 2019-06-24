
function index_matrix = get_index_matrix(a,d_1_p)

p = length(d_1_p);
d = sum(d_1_p);
d_a = d_1_p(a);
d_ma = d-d_a;

index_matrix = zeros(d_a*d_ma,p-1);

if a == 1
    for i = 1:p-1
        Bia_indicator = zeros(d_ma,d_a);
        Bia_indicator(sum(d_1_p(1:i))-d_a+1:sum(d_1_p(1:i+1))-d_a,1:d_a) = ...
                      ones(d_1_p(i+1),d_a);
        index_matrix(:,i) = Bia_indicator(:);
    end
    
elseif a>1 && a<p
    for i = 1:p-1
        Bia_indicator = zeros(d_ma,d_a);
        if i < a
            if i == 1
                Bia_indicator(1:d_1_p(i),1:d_a) = ...
                      ones(d_1_p(i),d_a);
            else
                Bia_indicator(sum(d_1_p(1:i-1))+1:sum(d_1_p(1:i)),1:d_a) = ...
                      ones(d_1_p(i),d_a);
            end
        elseif i == a
            Bia_indicator(sum(d_1_p(1:i-1))+1:sum(d_1_p(1:i+1))-d_a,1:d_a) = ...
                      ones(d_1_p(i+1),d_a);
        else
            Bia_indicator(sum(d_1_p(1:i))-d_a+1:sum(d_1_p(1:i+1))-d_a,1:d_a) = ...
                      ones(d_1_p(i+1),d_a);
        end              
        index_matrix(:,i) = Bia_indicator(:);
    end
    
else
    for i = 1:p-1
        Bia_indicator = zeros(d_ma,d_a);
        if i == 1
            Bia_indicator(1:d_1_p(i),1:d_a) = ...
                      ones(d_1_p(i),d_a);
        else
            Bia_indicator(sum(d_1_p(1:i-1))+1:sum(d_1_p(1:i)),1:d_a) = ...
                      ones(d_1_p(i),d_a);
        end
        index_matrix(:,i) = Bia_indicator(:);
    end
end

end

