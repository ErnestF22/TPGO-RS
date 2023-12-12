

function z_neighbor = find_neighbors(A_c,b_c,z,v,flag)
idx_neighbor=[];
V = v(1,1:2);
% if flag
%     v(3) = v(3)-100;
% end
for i = 1:size(A_c,1)
%     z(i,:)
%     A_c(i,:)*V'+b_c(i)
    if abs(A_c(i,:)*V'+b_c(i)-v(3))<=0.001
       idx_neighbor = [idx_neighbor;i];
    end
end
z_neighbor = z(idx_neighbor,:);
end