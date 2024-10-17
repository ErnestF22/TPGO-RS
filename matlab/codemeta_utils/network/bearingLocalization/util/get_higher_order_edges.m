function F = get_higher_order_edges(E,Adj)


nE = size(E,1); % number of edges

nV= size(Adj,1);% number of views

for ie=1:nE
    
    i = E(ie,1);
    j = E(ie,2);

    F{ie} = [];
    for k=1:nV
        if (i~=k) && (j~=k) && (Adj(i,k))  && (Adj(j,k)) 
            F{ie} = [F{ie} k];
        end
    end
end    


end

