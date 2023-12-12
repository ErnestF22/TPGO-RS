% clear 
% close all
% [obstacle_list,~,~,~,~,~,~,~,l1,l2,l3,l4] = paper_simulation_data;
% w = load('solved_paper_tree_rnd13_new.mat');
% T = w.T;
% T = fitLandmarks(T);
% t = -15;
% c =[130;110];
% figure(1)
% for i=1:size(obstacle_list,2)
%     plot_circle(obstacle_list(1,i),obstacle_list(2,i)...
%         ,obstacle_list(3,i),'b')
% end
% plot_landmarks(l1,'blue','s')
% hold on
% plot_landmarks(l2,'blue','s')
% hold on
% plot_landmarks(l3,'blue','s')
% hold on
% plot_landmarks(l4,'blue','s')
% hold on
% 
% Ss = ["r","m","g","b","k"];
% Sp = ["r.-","m.-","g.-","b.-","k.-"];
start = [146 48;146 86;109 90;153 145]';
% for i=1:4
%     plot(start(1,i),start(2,i),'p','MarkerSize',13,...
%     'MarkerEdgeColor',Ss(i),...
%     'MarkerFaceColor',Ss(i),'LineWidth',1.5)
%     hold on
%     k_path(i).idx = navigate(start(:,i),T,obstacle_list,1,Sp(i),l1,l2,l3,l4);
% end
% axis([0 185 0 170])
% grid on 
% xlabel('x')
% ylabel('y')

obstacle_list2 = rotating_obstacles(obstacle_list,c,t);
figure(2)
for i=1:size(obstacle_list2,2)
    plot_circle(obstacle_list2(1,i),obstacle_list2(2,i)...
        ,obstacle_list2(3,i),'b')
end
L1 = rotating_obstacles(l1,c,t)+[0 -1 -1;0 -3 -3];
L2 = rotating_obstacles(l2,c,t);
L3 = rotating_obstacles(l3,c,t)+[6 6 0;-4 -3 0];
L4 = rotating_obstacles(l4,c,t)+[-2 2 -1;0 0 0];%[3.3 3 0;0 0 0];
L1(:,1) =[0;0];                
L2(:,1) =[175;0];
L3(:,3)=[150;165];%[128;154]
start = rotating_obstacles(start,c,t);
plot_landmarks(L1,'b','s')
hold on
plot_landmarks(L2,'b','s')
hold on
plot_landmarks(L3,'b','s')
hold on
plot_landmarks(L4,'b','s')
hold on

Ss = ["r","m","g","b","k"];
Sp = ["r.-","m.-","g.-","b.-","k.-"];
for i=1:4
    plot(start(1,i),start(2,i),'p','MarkerSize',13,...
    'MarkerEdgeColor',Ss(i),...
    'MarkerFaceColor',Ss(i),'LineWidth',1.5)
    hold on
%     navigate(start(:,i),T,obstacle_list,2,Sp(i),l1,l2,l3,l4);
    deformed_navigation(start(:,i),k_path(i).idx,T,L1,L2,L3,L4,2,Sp(i))
end
axis([0 185 0 170])
grid on
xlabel('x')
ylabel('y')
axis equal