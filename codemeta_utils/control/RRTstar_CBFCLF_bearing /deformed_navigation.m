% navigate after deformation:
function deformed_navigation(start,k_path,tree,l1,l2,l3,l4,n,st)
x_n=start;
u=[100;100];
e = 0.1;
k_idx = k_path(1);
dx = 1;
U1 = 10;
DIS = 10;
while norm(U1)>0.001 && DIS>2.5
    x_p=x_n;
    y = tree(k_idx).L;
    DIS = norm(x_p-[0;0]);
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
    U1=u;
%     if dis>2
%         u = u/norm(u);
%     end
    x_n=x_p+e*u;
    if (norm(U1)<0.1) && (dx<size(k_path,2))
        dx = dx +1;
        k_idx  = k_path(dx);
    end
    figure(n)
    plot(x_n(1),x_n(2),st,'MarkerSize',4)
    hold on
end
norm(u)
end
