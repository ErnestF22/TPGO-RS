clear all
close all

nodes = [10 30;20 60;45 80;90 60;70 20;30 0;30 40;40 55;50 58;65 50;60 30;40 25];
sets = [ 8 7 1 2;
    9 8 2 3;
    10 9 3 4;
    11 10 4 5;
    12 11 5 6;
    7 12 6 1];
node_set = sets(1,:);
y1 = nodes(node_set,:)';
node_set = sets(2,:);
y2 = nodes(node_set,:)';
L = y1;

Y(1).y = y1;
Y(2).y = y2;
M = struct([]);
flag_plot = 1;

if flag_plot
    
    for t = 1:2
        y = Y(t).y;
        for i=1:size(y,2)
            plot(y(1,i),y(2,i),'b*')
            hold on
        end
        
        for i=1:size(L,2)
            plot(L(1,i),L(2,i),'g*')
            hold on
        end
        
        plot( y(1,:), y(2,:), '.-r');
        hold on
        plot( [y(1,end) y(1,1)],[y(2,end) y(2,1)], '.-r');
    end
end
XE = [];
for i=1:2
    y = Y(i).y;
    p = (y(:,2)+y(:,4))/2;
    [A,b] = convexSet(y);
    [Ax,bx] = LPSet(A,b,p);
    Ah = -Ax;
    bh = bx;
    Ah(4,:)=[];
    bh(4)=[];
    exitDir =  -Ax(4,:)';
    % exitDir = exitDir/norm(exitDir);
    
    % xe = intersectionTwoLines(A(2,1),b(2),A(4,1),b(4));
    xe = (y(:,1)+y(:,4))/2;
    XE = [XE xe];
    M(i).Ax = Ax;
    M(i).bx = bx;
    M(i).Ah = Ah;
    M(i).bh = bh;
    M(i).xe = xe;
    M(i).z = exitDir;
end

epsilon = 0.7;
E = [0;0];

Cb = 1;
Cl = 0.5;
% 
xc1 = y1(:,1);
xc2 = y1(:,4);
plot(xc1(1),xc1(2),'m*')
hold on
plot(xc2(1),xc2(2),'m*')
hold on 

K_new = optFirstorderU_newCost(Cb,Cl,L,M,xc1,xc2);
Wb = 1*ones(size(Ah,1),1);
Wl = 1;
for i=1:2
    [K,~,~,~,~,k_added] = optFirstorderWithU(Wb,Wl,Cb,Cl,M(i).z,M(i).Ah,M(i).Ax,M(i).bx,M(i).bh,L,XE(:,i));
    K_old(i).K = K;
    K_old(i).added = k_added;
end
figure(1)
start = [20;40];
Ye(1).y = y1(:,end);
Ye(2).y = y2(:,end);
navigate(start,K_new,L,XE,Ye,'b*');
navigate(start,K_old,L,XE,Ye,'r*');
title('red: old, blue: new')
axis equal
return 

for i=1:2
    k = K_new(i).K;
    k_added = K_new(i).added;
    Ax = M(i).Ax;
    bx = M(i).bx;
    xe = M(i).xe;
    plot_controllers(Ax,bx,k,k_added,L,epsilon,0,E)
    plot(xe(1),xe(2),'c*')
    hold on;
    grid on
    
end

