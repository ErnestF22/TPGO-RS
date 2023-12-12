clear 
close all
w = load('solved_experiment_tree.mat');
T = w.T;
[~,obs,obs_offset,l1,l2,l3,l4] = experiment_tree;
T = fitLandmarks_exp(T);
T(9).L = 3;
t = -10;
c =[300;200];
figure(1)
plot_exp(T,obs,obs_offset,l1,l2,l3,l4)
start = [245 155;240 265;263 360;165 265 ]';%
Ss = ["r","m","g","b","k"];
Sp = ["r.-","m.-","g.-","b.-","k.-"];
for i=1:4
    %     '`MarkerEdgeColor',Ss(i),...
    plot(start(1,i),start(2,i),'p','MarkerSize',13,...
    'MarkerFaceColor',Ss(i),'LineWidth',1.5)
    hold on
    k_path(i).idx = navigate_exp(start(:,i),T,1,Sp(i),l1,l2,l3,l4);
                                 
end

grid on 
xlabel('x')
ylabel('y')

obs2 = rotating_obstacles(obs,c,t);
obs_offset2 = rotating_obstacles(obs_offset,c,t);

l1 = rotating_obstacles(l1,c,t)+[0 -2 -5;0 -25 -25];%[0 -9 -9;0 -30 -55];
l2 = rotating_obstacles(l2,c,t);
l3 = rotating_obstacles(l3,c,t)+[0 -10 0;20 40 0];
l4 = rotating_obstacles(l4,c,t)+[0 0 0;0 0 0];
l1(:,1) = [0;0];
l2(:,1) = [250;0];
axis equal
figure(2)
plot_exp(T,obs2,obs_offset2,l1,l2,l3,l4)

start = rotating_obstacles(start,c,t);
for i=1:4
    plot(start(1,i),start(2,i),'p','MarkerSize',13,...
    'MarkerEdgeColor',Ss(i),...
    'MarkerFaceColor',Ss(i),'LineWidth',1.5)
    hold on
    deformed_navigation(start(:,i),k_path(i).idx,T,l1,l2,l3,l4,2,Sp(i))
end
grid on
xlabel('x')
ylabel('y')
axis equal