function k_path = navigate_exp(start,tree,n,st,l1,l2,l3,l4)
x_n=start;
u=[100;100];
e = 0.1;
c=0;
if start(1)==215
    k_idx = 13;
else
k_idx  = find_nearest_k(tree,x_n,[]);
end
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
    Y=M;
    u=k*Y;
%     %%%%%%%%%%%%
%     K1 = k(1:2,1:2);
%     K2 = k(1:2,3:4);
%     K3 = k(1:2,5:6);
%     P1 = y(:,1);
%     P2 = y(:,2);
%     P3 = y(:,3);
%     K_bias = K2*(P2-P1)+K3*(P3-P1);
%     u = (K1+K2+K3)*(P1-x_p)+K_bias;
    %%%%%%%%%%%%%
    dis = norm(x_n-tree(tree(k_idx).parent).position);
    dis = norm(x_n-tree(k_idx).CLF.xe);
    x_n=x_p+e*u;
    %x_n-xp
    if (norm(u)<1) && ((tree(k_idx).parent)~=1)
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
