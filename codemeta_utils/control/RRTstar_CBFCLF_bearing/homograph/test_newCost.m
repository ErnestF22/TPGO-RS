clear all
close all

y1 = [2 3;3.5 8;11.5 6.1;11.4 4.5]';
y2 = [11.4 4.5;11.5 6.1;17 3;15 2.5]';
L = [12 8;12 1;16 10;18 1]'+[5;0];

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

for i=1:2
    y = Y(i).y;
    p = (y(:,2)+y(:,4))/2;
    [A,b] = convexSet(y);
    [Ax,bx] = LPSet(A,b,p);
    Ah = -Ax;
    bh = bx;
    Ah(3,:)=[];
    bh(3)=[];
    exitDir =  -Ax(3,:)';
    % exitDir = exitDir/norm(exitDir);
    
    % xe = intersectionTwoLines(A(2,1),b(2),A(4,1),b(4));
    xe = (y(:,3)+y(:,4))/2;
    
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

xc1 = [11.4;4.5];
xc2 = [11.5;6.1];

K = optFirstorderU_newCost(Cb,Cl,L,M,xc1,xc2);
% display(K)
% display(k_added)
figure(1)
for i=1:2
    k = K(i).K;
    k_added = K(i).added;
    Ax = M(i).Ax;
    bx = M(i).bx;
    xe = M(i).xe;
    plot_controllers(Ax,bx,k,k_added,L,epsilon,0,E)
    plot(xe(1),xe(2),'c*')
    hold on;
    grid on
    
end

