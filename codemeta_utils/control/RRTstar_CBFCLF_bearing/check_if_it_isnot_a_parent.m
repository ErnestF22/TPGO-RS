% This function check if node_i in a tree T is a leaf (a node with degree
% 1) or not. If it is a leaf, then it returns True. 
function flag = check_if_it_isnot_a_parent(T,i)
flag = 1;
for j =1: size(T,2)
    if T(j).parent == i
        flag =0;
    end
end
end
