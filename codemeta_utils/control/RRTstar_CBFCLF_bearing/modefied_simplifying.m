function T = modefied_simplifying(T,i,j)
p1 = round(T(i).position,1);
q1 = round(T(T(i).parent).position,1);
p2 = round(T(j).position,1);
q2 = round(T(T(j).parent).position,1);
array = [p1 q1 p2 q2];
x = find_intersecition(p1,q1,p2,q2);
x = round(x,1);
flag = ismember(array,x);
flag = any(flag(1,:).*flag(2,:));
if ~flag
    T(end+1).position = x;
    c_i = norm(x-q1)+T(T(i).parent).cost;
    c_j = norm(x-q2)+T(T(j).parent).cost;
    if c_i <= c_j
        T(end).parent = T(i).parent;
        T(end).cost = c_i;
        T(end).children = [i,j];
        parent_i_children = T(T(i).parent).children;
        idx_i = find(parent_i_children==i);
        parent_i_children(idx_i)=size(T,2);
        T(T(i).parent).children = parent_i_children;
        
        parent_j_children = T(T(j).parent).children;
        idx_j = find(parent_j_children==j);
        parent_j_children(idx_j)=[];
        T(T(j).parent).children = parent_j_children;
    else
        T(end).parent = T(j).parent;
        T(T(j).parent).children = [T(T(j).parent).children,size(T,2)];
        T(end).cost = c_j;
        T(end).children = [i,j];
        parent_j_children = T(T(j).parent).children;
        idx_j = find(parent_j_children==j);
        parent_j_children(idx_j)=size(T,2);
        T(T(j).parent).children = parent_j_children;
        parent_i_children = T(T(i).parent).children;
        idx_i = find(parent_i_children==i);
        parent_i_children(idx_i)=[];
        T(T(i).parent).children = parent_i_children;
    end
    
    T(i).parent = size(T,2);
    T(j).parent = size(T,2);
    
    
end
end