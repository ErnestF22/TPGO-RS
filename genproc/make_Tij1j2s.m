function [Tij1j2, Tij1j2_tilde] = make_Tij1j2s(node_id, R,T,Tijs,edges,params)

nrs = params.nrs;
d = params.d;
node_deg = params.node_degrees(node_id); % usually low_deg

Tij1j2 = zeros(d, node_deg);
Tij1j2_tilde = zeros(nrs, node_deg);

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


%additional checks

disp("Tij1j2 inside make_Tij1j2s()");
disp(Tij1j2);

disp("Tij1j2_tilde inside make_Tij1j2s()");
disp(Tij1j2_tilde);

R_i = R(:,:,node_id);
disp("[R_i * Tij1j2, Tij1j2_tilde] inside make_Tij1j2s()");
disp([R_i * Tij1j2, Tij1j2_tilde]);

disp("max(abs(R_i * Tij1j2 - Tij1j2_tilde), [], ""all"")");
disp(max(abs(R_i * Tij1j2 - Tij1j2_tilde), [], "all"));

end %file function