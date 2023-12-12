% Lab global reference! 
clear
close all
w = load('solved_experiment_tree_global_6.mat');
T = w.T;
org = T(1).position;
[~,obs,obs_offset,l1,l2,l3,l4] = experiment_tree_global;

t = -25;
c =[90;-240];
figure(1)
plot_exp(T,obs,obs_offset,l1,l2,l3,l4)
plot(org(1),org(2),'r*');
start = [130 -200;130 -315;110 -380;215 -300]';%
Ss = ["r","m","g","b","k"];
Sp = ["r.-","r.-","g.-","b.-","k.-"];
% for i=1:4
%         '`MarkerEdgeColor',Ss(i),...
%     plot(start(1,i),start(2,i),'p','MarkerSize',13,...
%     'MarkerFaceColor',Ss(i),'LineWidth',1.5)
%     hold on
%     k_path(i).idx = navigate_exp(start(:,i),T,1,Sp(i),l1,l2,l3,l4);
%                                  
% end

grid on 
xlabel('x')
ylabel('y')

obsd_deformed = obs; 
obsd_deformed(:,3:6)= obsd_deformed(:,3:6)+[0 0 0 0;100 100 100 100];
obsd_deformed(:,end-3:end)= obsd_deformed(:,end-3:end)+[0 0 0 0;10 10 10 10];
obs_offset2 = obs_offset;
obs_offset2(:,3:6)=obs_offset2(:,3:6)+[0 0 0 0;100 100 100 100];
obs_offset2(:,end-3:end)= obs_offset2(:,end-3:end)+[0 0 0 0;10 10 10 10];
obsd_deformed = rotating_obstacles(obsd_deformed,c,t);
obs_offset2 = rotating_obstacles(obs_offset2,c,t);

deformed_l1 = l1+[-15 0 0 0;0 0 50 50];
deformed_l1 = rotating_obstacles(deformed_l1,c,t);
deformed_l2 = rotating_obstacles(l2,c,t);
deformed_l3 = rotating_obstacles(l3,c,t)+[0 0 0 -50;0 0 0 30];
deformed_l4 = l4;
deformed_l4(:,2:3) = deformed_l4(:,2:3)+[0 0;100 100];
deformed_l4 = rotating_obstacles(deformed_l4,c,t);

deformed_l4(:,1)=deformed_l4(:,1)+[100;100];

deformed_l1(:,1) = l1(:,1)+[5;65];%for 1:2 +[5;45]
deformed_l1(:,4) = l1(:,4)+[-15;-70];
deformed_l1(:,3) = l1(:,3)+[-5;40];
deformed_l1(:,2) = l1(:,2)+[30;0];

deformed_l2(:,1) = deformed_l2(:,1)+[0;70];
deformed_l2(:,2) = deformed_l2(:,2)+[17;40];
deformed_l2(:,3) = deformed_l2(:,3)+[0;5];
 
deformed_l3(:,2) = deformed_l3(:,2)+[20;-100];
deformed_l3(:,1) = deformed_l3(:,1)+[30;40]; 
deformed_l3(:,4) = deformed_l3(:,4)+[50;0]; 

% deformed_l3(:,1) = deformed_l3(:,1)+[30;40]; 
% deformed_l3(:,2) = deformed_l3(:,2)+[20;-100];
% % 
% % deformed_l3(:,2) = deformed_l3(:,2)+[-30;-60];
% % deformed_l3(:,1) = deformed_l3(:,1)+[30;40]; 

axis equal
figure(2)
plot_exp(T,obsd_deformed,obs_offset2,deformed_l1,deformed_l2,deformed_l3,deformed_l4)
plot(org(1),org(2),'r*');
start = rotating_obstacles(start,c,t);
start(:,4)= start(:,4)+[50;100];
return 
k_path(2).idx = [7,6,4,3];
for i=2:2
    plot(start(1,i),start(2,i),'p','MarkerSize',13,...
    'MarkerEdgeColor',Ss(i),...
    'MarkerFaceColor',Ss(i),'LineWidth',1.5)
    hold on
    deformed_navigation(start(:,i),k_path(i).idx,T,deformed_l1,deformed_l2,deformed_l3,deformed_l4,2,Sp(i))
end
grid on
xlabel('x')
ylabel('y')
axis equal