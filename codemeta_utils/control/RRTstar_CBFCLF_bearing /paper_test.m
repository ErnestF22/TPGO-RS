% this function simulates the exmaple for the paper:

rng(11)
clear
close all

[obstacle_list,start,goal,max_itr,lim_x,lim_y,expand_dis,...
    path_resolution,l1,l2,l3,L] = iros22_paper_lab_data;
% 
% 
% [obstacle_list,start,goal,max_itr,lim_x,lim_y,expand_dis,...
%     path_resolution,l1,l2,l3,l4] = iros22_paper_simulation_data;
% % 
[~,obs,~,~,~]= main_test_RRTs(obstacle_list,start,goal,...
   max_itr,lim_x,lim_y,expand_dis,...
    path_resolution,0);
T = lab_test_iros2022;
% T =  test_iros2022;

% figure(1)
% plot_simplified_tree(T,'b',1.5,1); 
% figure(1)
% I = imread('Lab_Cspace.JPG');
% h = image([0 500],[-15 600],I); 
% uistack(h,'bottom')
% hold on
% axis equal
% axis ([0 500 -15 600])
% hold on
% T = lab_test_iros2022;
% % T =  test_iros2022;
% return 
figure(1)
plot_simplified_tree(T,'b',1.5,2); 
% hold on 
% I = imread('plan2.jpg');
% h = image([0 100],[0 100],I); 
% uistack(h,'bottom')
% hold on
% axis ([0 100 0 100])



% 
% % T = fitLandmarks(T);
% % 
T = voronoi(T);
T = nearestObstacle(T, obs);
T = nearest_angle_Obstacle(T, obs);
T = CLF_constraint(T);
T = CBF_constraint(T,obs);
T = find_controller(T,l1,l2,l3,L,0);
save('new_solved_paper_lab_test_with_offset_tree_rnd11_change_idx_correct_index_3March.mat')
return
% load ('solved_paper_lab_test_with_offset_tree_rnd11.mat')

% load ('new_solved_paper_lab_test_with_offset_tree_rnd11_change_idx.mat')
load('new_solved_paper_lab_test_with_offset_tree_rnd11_change_idx_correct_index_3March.mat')

figure(2)
start_1 = [50;40];
start_1 = [100;550];
k_path_1 = navigate(start_1,T,obstacle_list,2,'r*',l1,l2,l3,L,0);
hold on 

start_3 = [44;42];
start_3 = [300;375];
k_path_3 = navigate(start_3,T,obstacle_list,2,'m*',l1,l2,l3,L,0);

start_4 = [90;90];
start_4 = [220;550];
k_path_4 = navigate(start_4,T,obstacle_list,2,'c*',l1,l2,l3,L,0);

start_2 = [10;90];
start_2 = [500;590];
k_path_2 = navigate(start_2,T,obstacle_list,2,'g*',l1,l2,l3,L,0);
hold off

save('new_solved_paper_lab_test_with_offset_tree_rnd11_change_idx_correct_index_3March.mat')