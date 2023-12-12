% this function simulates the exmaple for the paper:
rng(0)
rng(13)
clear
close all
[obstacle_list,start,goal,max_itr,lim_x,lim_y,expand_dis,...
    path_resolution,l1,l2,l3,l4] = paper_simulation_data;

[~,obs,T,~,t,~]= main_test_RRTs(obstacle_list,start,goal,...
    max_itr,lim_x,lim_y,expand_dis,...
    path_resolution);

T = fitLandmarks(T);

T = voronoi(T);
figure(3)
plot_bisector_segment(T)
T = nearestObstacle(T, obs);
T = nearest_angle_Obstacle(T, obs);
T = CLF_constraint(T);
T = CBF_constraint(T,obs);
T = find_controller(T,l1,l2,l3,l4);
% save('solved_paper_tree_rnd13_new.mat')

% figure(2)
% start = [113;130];
% navigate(start,T,obstacle_list,1,'r*');
% figure(2)
% % plot_controllers(T,'b');
