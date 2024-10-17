function [Tij1j2, Tij1j2_tilde] = make_Tij1j2s(node_id, ~,T,Tijs,edges,params)

num_rows_T = size(T,1);
% nrs = params.nrs;
d = params.d;
node_deg = params.node_degrees(node_id); % usually low_deg

Tij1j2 = zeros(d, node_deg);
Tij1j2_tilde = zeros(num_rows_T, node_deg);

found = 1;
for e = 1:size(edges,1)
    e_i = edges(e,1);
    e_j = edges(e,2);
    if e_i == node_id
        Tij1j2(:,found) = Tijs(:,e);
        Tij1j2_tilde(:,found) = T(:,e_j) - T(:,e_i); 
        found = found + 1;
    end

end


end %file function