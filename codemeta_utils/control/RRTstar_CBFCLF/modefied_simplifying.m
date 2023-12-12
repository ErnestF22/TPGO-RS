function T = modefied_simplifying(T,i,j)
p1 = round(T(i).position,3);
q1 = round(T(T(i).parent).position,3);
p2 = round(T(j).position,3);
q2 = round(T(T(j).parent).position,3);
array = [p1 q1 p2 q2];
x = find_intersecition(p1,q1,p2,q2);
x = round(x,3);
flag = ismember(array,x);
flag = any(flag(1,:).*flag(2,:));
if ~flag
    T(end+1).position = x;
    c_i = norm(x-q1)+T(T(i).parent).cost;
    c_j = norm(x-q2)+T(T(j).parent).cost;
    if c_i <= c_j
        T(end).parent = T(i).parent;
        T(end).cost = c_i;
    else
        T(end).parent = T(j).parent;
        T(end).cost = c_j;
    end
    T(i).parent = size(T,2);
    T(j).parent = size(T,2);
end
end