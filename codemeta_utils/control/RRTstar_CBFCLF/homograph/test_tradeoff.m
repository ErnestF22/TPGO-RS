clear
close all


%% Settings:
epsilon = 5; % step size to plot the controller
% refder to the optimization problem in the paper
Wb = 1*ones(3,1);
Wl = 1;

%% Environment
y = [134 25;17 10;15 80;115 60]';
L = y;
%% Find controller:
res = struct([]);
Cl = 0.0001;
Cb = 1;
s = 1.3;

%Cb = 1,k = 2
%cl_s = [4 1.01;3 1.1;2 1.18;1 1.25;0.5 1.3;0.1 1.36;0.05 1.4;0.01 1.55]

%cl = 0.01, cb = 1 
%c_s = [1 1.36;2 1.55;10 1.59]

p = (y(:,2)+y(:,4))/2; % point in a convex set
[A,b] = convexSet(y); % Ax = b
[Ax,bx] = LPSet(A,b,p); % Ax<b
%plot_convex(Ax,bx)

% barrier function A_hx+b_h>0
Ah = -Ax;
bh = bx;
Ah(4,:)=[];
bh(4)=[];

exitDir =  (-Ax(4,:)');
xe = (y(:,1)+y(:,4))/2;

plot([xe(1) xe(1)+exitDir(1)],[xe(2) xe(2)+exitDir(2)],'-b')

hold on
plot(xe(1)+exitDir(1),xe(2)+exitDir(2),'*b')
plot_env(y,xe,L)
hold on


[K,k_added] = find_controller(exitDir,Wb,Wl,Ah,Ax,bx,bh,xe,y,L,Cl,Cb,0,1,s);

plot_controllers(Ax,bx,K,k_added,L,epsilon,1,0,1)
plot([xe(1) xe(1)+k_added(1)],[xe(2) xe(2)+k_added(2)],'-r')
% plot_trajectory([20,40]')
% plot_trajectory(1,1)
% plot_trajectory(1/1.2,2)
% plot_trajectory(1.5,3)