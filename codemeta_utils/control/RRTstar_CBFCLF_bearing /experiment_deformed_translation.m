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
    plot(start(1,i),start(2,i),'p','MarkerSize',13,...
    'MarkerEdgeColor',Ss(i),...
    'MarkerFaceColor',Ss(i),'LineWidth',1.5)
    hold on
    k_path(i).idx = navigate_exp(start(:,i),T,1,Sp(i),l1,l2,l3,l4);
                                 
end
grid on 
xlabel('x')
ylabel('y')
axis equal


obs2 = [183 90;
    214.5 90;
    214.5 191.5;
    279.5 191.5;
    279.5 238;
    214.5 238;
    214.5 368;
    183 368;
    183 318-88;
    91 318-88;
    91 286.5-88;
    183 286.5-88]';

obs_offset2 = [167.25 74.25;
    230.25 74.25;
    230.25 175.75;
    295.25 175.75;
    295.25 253.75;
    230.25 253.75;
    230.25 383.75;
    167.25 383.75;
    167.25 333.75-88;
    75.25 333.75-88;
    75.25 270.75-88;
    167.25 270.75-88]';

l4 = l4 + [0 0 0;0 -88 -88];
l1 = l1 + [0 0 0;0 -30 -70];

figure(2)
plot_exp(T,obs2,obs_offset2,l1,l2,l3,l4)

% start = rotating_obstacles(start,c,t);
start(1,end)=start(1,end)-1;
for i=3:4
    plot(start(1,i),start(2,i),'p','MarkerSize',13,...
    'MarkerEdgeColor',Ss(i),...
    'MarkerFaceColor',Ss(i),'LineWidth',1.5)
    hold on
     navigate_exp(start(:,i),T,2,Sp(i),l1,l2,l3,l4);
%     deformed_navigation(start(:,i),k_path(i).idx,T,l1,l2,l3,l4,2,Sp(i))
end
grid on
xlabel('x')
ylabel('y')
axis equal