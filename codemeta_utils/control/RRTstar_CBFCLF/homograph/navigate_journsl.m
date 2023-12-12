function k_path = navigate_journsl(start,K,M,st,flag_loop,fig)
figure(fig)
x_n=start;
u=[100;100];
e = 0.1;
idx  = 1;
k_path = idx;
loop = 0;
if flag_loop
    counter = 5;
else
    counter =1;
end
while norm(u)>10^-3 && loop<counter
    x_p=x_n;
    k = K(idx).K;
    y = M(idx).y;
    L = M(idx).L;
    Y=(L-x_p*ones(1,size(L,2)));
    Y = reshape(Y,[],1);
    u = k*Y+K(idx).added;
    if norm(x_p-[28;25])<0.07 && ~flag_loop
        break
    end
    u = u/norm(u);
    a = y(:,1)-x_p;
    b = y(:,1)-M(idx).xe;
    crossP = cross([a;1],[b;1]);
    sin_a  = abs(crossP(3)/(norm(a)*norm(b)));
    dis = norm(a)*sin_a;
    x_n=x_p+e*u;
    if (dis<0.1) && idx<=6
%         plot([x_p(1) y(1,4)],[x_p(2) y(2,4)],'b-')
%         hold on
%         plot([y(1,4) M(idx).xe(1)],[y(2,4) M(idx).xe(2)],'g-')
%         hold on
        idx  = idx+1;
        k_path = [k_path idx];
    end
    if idx>6
        %break
        idx = 1;
        loop = loop+1;
    end
    plot(x_n(1),x_n(2),st,'MarkerSize',5)
    hold on
end
end

