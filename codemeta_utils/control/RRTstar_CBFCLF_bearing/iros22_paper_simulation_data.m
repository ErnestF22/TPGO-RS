%this function generates inputs for paepr simulation
function [obs,start,goal,max_itr,lim_x,lim_y,expand_dis,path_resolution,l1,l2,l3,l4] = iros22_paper_simulation_data
start = [0;0];
goal = [100;100];
max_itr = 300;
lim_x = [0 100];
lim_y = [0 100];
expand_dis = 50;
path_resolution = 1;

%obstacle = [x y radius]
obs = [0 60 20
       100 0 30
       30 25 15
       61 95 11
       61 85 11
       61 75 11
       61 65 11
       61 55 11]';

l1 = [73 95;72 85;73.5 75;74 65;72.5 55;62 44;100 30;80 22;72.5 10]';
l2 = [52 95;51.5 85;51 75;50 65;50.5 55;0 80;17.5 70;20 60;19.37 55]';
l3 = [71 0;54 47.47;60 0;40 0; 0 0;39 12;41.55 34.5;20 13.5]';
l4 = [0.1 0;0.2 10;0 20;0.3 30;0.5 40;20 36;31 40]';

end