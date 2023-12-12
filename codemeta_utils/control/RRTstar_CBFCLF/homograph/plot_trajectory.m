function plot_trajectory()
close all
clear 
start = [10,30]';
load('res.mat')
x_n = start;
c = 0;
e = 0.5;
sty = ["b*","b*","g*","m*","c*","y*"];
for i=1:6
    L = res(i).landmarks;
    y = res(i).y;
    k = res(i).K;
    xe = res(i).xe;
    plot_env(y,xe,y)
    hold on
%     plot_controllers(res(i).Ax,res(i).bx,res(i).K,res(i).k_added,res(i).landmarks,2,1,0,1)
end
% set(gca,'XTick',[], 'YTick', [])
i = 1;
flag =1;
while i<7 && flag
    L = res(i).landmarks;
    k = res(i).K;
    xe = res(i).xe;
    smax = res(i).smax;
    x_p=x_n;
    k_added = res(i).k_added;
    Y=(L-x_p*ones(1,size(L,2)));
    Y = reshape(Y,[],1);
    S = (2-1/2)*rand(1,4)+1/2;
    r = randperm(4);
    I = ones(1,8);
    for j=1:4
        I(r(j)) = S(j);
    end
    S = I;
    S = diag(S);
    Y = S*Y;
    u=k*Y+k_added;
    u = u/norm(u);
    x_n=x_p+e*u;
    plot(x_n(1),x_n(2),"b*",'MarkerSize',0.8)
    hold on
    d = point_to_line([x_n;0], [L(:,1);0], [L(:,4);0]);
    if d<0.5
        i=i+1;
        if i>6
            i=1;
            c = c+1;
            if c == 5
                flag = 0;
            end
        end
    end
end
% set(gca,'XTick',[], 'YTick', [])
end