clear all
close all

% [T,obs,obs_offset,l1,l2,l3,l4] = experiment_tree;
[T,obs,obs_offset,l1,l2,l3,l4] = experiment_tree_global;
% T = fitLandmarks_exp(T);

T = voronoi(T);
T = nearestObstacle(T, obs);
T = CLF_constraint(T);
T = CBF_constraint(T, obs);
T = find_controller(T,l1,l2,l3,l4);
save('solved_experiment_tree_global_6.mat')
% plot_exp(T,obs,obs_offset,l1,l2,l3,l4)

% start = [245 155;240 265;263 360;165 265 ]';%
% start = [240 265]';
% Ss = ["r","m","g","b","k"];
% Sp = ["r.-","m.-","g.-","b.-","k.-"];
% for i=1:4
%     plot(start(1,i),start(2,i),'p','MarkerSize',13,...
%     'MarkerEdgeColor',Ss(i),...
%     'MarkerFaceColor',Ss(i),'LineWidth',1.5)
%     hold on
%     navigate_exp(start(:,i),T,1,Sp(i),l1,l2,l3,l4);
%                             
% end