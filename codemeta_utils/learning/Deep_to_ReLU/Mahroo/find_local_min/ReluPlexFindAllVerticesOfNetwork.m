function [V_total,Z_total] = ReluPlexFindAllVerticesOfNetwork(V,Ab_set,L,flag)
z=[];
for i=1:L
    z = [z (Ab_set(i).z)'];
end

Q = priority_prepare();
for i=1:size(V,1)
    Q = priority_insert(Q,V(i,:),z);
end

[A_c,b_c,z] = compute_A_cum_b_cum_all_combinations(Ab_set,L,size(z,2));
V_total = [];
Z_total = [];
while ~isempty(Q)
    [Q,v_t,z_t]=priority_minExtract(Q);
    v_t = round(v_t,4);
    V_total = [V_total;v_t];
    Z_total = [Z_total;z_t];
    
    z_neighbor = find_neighbors(A_c,b_c,z,v_t,flag);
    
    for j = 1:size(z_neighbor,1)
        if ~ismember(z_neighbor(j,:),Z_total,'rows')
%             Z_total = [Z_total;z_neighbor(j,:)];
%             V_total = [V_total;v_t];
%             Q = priority_insert(Q, v_t,z_neighbor(j,:));
            V_new = ReluPlexFindAllVerticesOfRegion((z_neighbor(j,:))',Ab_set,L);
            V_new = round(V_new,4);
            for i=1:size(V_new,1)
                if ~ismember(V_new(i,:),V_total,'rows')
                    Q = priority_insert(Q, V_new(i,:),z_neighbor(j,:));
                    Z_total = [Z_total;z_neighbor(j,:)];
                    V_total = [V_total;V_new(i,:)];
                end
            end
        end
    end
end

end