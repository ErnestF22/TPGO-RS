function k_path = navigate(start,tree,obs,n,st,l1,l2,l3,l4)
x_n=start;
u=[100;100];
e = 0.5;
c=0;

k_idx  = find_nearest_k(tree,x_n,obs);
k_path = k_idx;
while norm(u)>10^-3
    c = c+1;
    x_p=x_n;
    y = tree(k_idx).L;
    if y == 1
        y = l1;
    elseif y==2
        y = l2;
    elseif y==3
        y = l3;
    else
        y = l4;
    end
    k = tree(k_idx).K;
    Y=(y-x_p*ones(1,size(y,2)));
    M = zeros(2*size(y,2),1);
    for s=1:size(Y,2)
        M(2*s-1,1)=Y(1,s);
        M(2*s,1)=Y(2,s);
    end
    dis = norm(x_n-tree(tree(k_idx).parent).position);
    Y=M;
    u=k*Y;
    if dis>1
        u = u/norm(u);
    end
    x_n=x_p+e*u;
    if (dis<10^0.2) && ((tree(k_idx).parent)~=1)
        k_idx  = tree(k_idx).parent;
        k_path = [k_path k_idx];
        norm(u)
    end
    figure(n)
    plot(x_n(1),x_n(2),st,'MarkerSize',4)
    hold on
%     figure(3)
%     plot(c,dis,st)
%     hold on
end
norm(u)
end
