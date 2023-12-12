function k_path = navigate(start,K,L,XE,y,st)
x_n=start;
u=[100;100];
e = 0.1;
c=0;

idx  = 1;
k_path = idx;
while norm(u)>10^-3 && c<1000
    c = c+1;
    x_p=x_n;
    k = K(idx).K;
    Y=(L-x_p*ones(1,size(L,2)));
    Y = reshape(Y,[],1);
    u=k*Y+K(idx).added;
    u = u/norm(u);
    a = y(idx).y-x_p;
    b = y(idx).y-XE(:,idx);
    crossP = cross([a;1],[b;1]);
    sin_a  = abs(crossP(3)/(norm(a)*norm(b)));
    dis = norm(a)*sin_a;
    x_n=x_p+e*u;
    if (dis<0.1) && idx<=2
        plot([x_p(1) y(idx).y(1)],[x_p(2) y(idx).y(2)],'r-')
        hold on
        plot([y(idx).y(1) XE(1,idx)],[y(idx).y(2) XE(2,idx)],'g-')
        hold on
        idx  = idx+1;
        k_path = [k_path idx];
    end
    if idx>2
        break
    end
    plot(x_n(1),x_n(2),st,'MarkerSize',2)
    hold on
end
end

