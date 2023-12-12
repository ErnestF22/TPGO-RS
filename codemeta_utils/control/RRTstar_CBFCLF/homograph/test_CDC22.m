clear all
close all
sty = ["k*","b*","g*","m*","c*","y*"];

%% Settings:
epsilon =5; % step size to plot the controller
% refder to the optimization problem in the paper
Wb = 1*ones(3,1);
Wl = 1;

% 0: translation
% 1: bearing
flag_controller = 1;
flag_navigation = 1;
%% Environment
nodes = [10 30;20 60;45 80;90 60;70 20;30 0;30 40;40 55;50 58;65 50;60 30;40 25];
sets = [ 8 7 1 2;
    9 8 2 3;
    10 9 3 4;
    11 10 4 5;
    12 11 5 6;
    7 12 6 1];
%% Find controller:
res = struct([]);
Cl = [0.001,0.001,0.0001,0.001,0.001,0.001];
Cb = [100,100,100,100,100,100];%[20,500,1000,1,10,200];
s = [1.43, 1.43,1.29,1.057,1.2,1.5];%1.5,10,0.1
for i=1:length(s)
    node_set = sets(i,:);
    y = nodes(node_set,:)';
    
    p = (y(:,2)+y(:,4))/2; % point in a convex set
    [A,b] = convexSet(y); % Ax = b
    [Ax,bx] = LPSet(A,b,p); % Ax<b
    %plot_convex(Ax,bx)
    
    % barrier function A_hx+b_h>0
    Ah = -Ax;
    bh = bx;
    Ah(4,:)=[];
    bh(4)=[];
%     
    exitDir =  (-Ax(4,:)');
%     xe = intersectionTwoLines(A(1,1),b(1),A(3,1),b(3));
    xe = (y(:,1)+y(:,4))/2;
%     xe = p;
    plot([xe(1) xe(1)+exitDir(1)],[xe(2) xe(2)+exitDir(2)],'-b')
    hold on
    plot(xe(1)+exitDir(1),xe(2)+exitDir(2),'*b')
    a = [0.818380832892274,0.394046627022142,0.575430144636146,0.684437372895086;0.753466398507143,0.334776309440223,0.312565283713393,0.632149144742378];
    L = y+a;
    plot_env(y,xe,y)
%     hold on
    
   
    [K,k_added] = find_controller(exitDir,Wb,Wl,Ah,Ax,bx,bh,xe,y,L,Cl(i),Cb(i),0,flag_controller,s(i));
    
    plot_controllers(Ax,bx,K,k_added,L,epsilon,flag_navigation,0,1)
    
    res(end+1).K = K;
    res(end).k_added = k_added; 
    res(end).xe = xe;
    res(end).Ax = Ax;
    res(end).bx = bx;
    res(end).landmarks = L;
    res(end).y = y;
    res(end).smax = s(i);
    
end

% save('res')
% plot_trajectory([20,40]')
% plot_trajectory()
% plot_trajectory(1/1.2,2)
% plot_trajectory(1.5,3)