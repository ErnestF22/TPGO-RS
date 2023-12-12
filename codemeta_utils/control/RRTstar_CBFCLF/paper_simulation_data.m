%this function generates inputs for paepr simulation
function [obs,start,goal,max_itr,lim_x,lim_y,expand_dis,path_resolution,l1,l2,l3,l4] = paper_simulation_data
start = [0;0];
goal = [160;130];
max_itr = 500;
lim_x = [0 180];
lim_y = [0 170];
expand_dis = 60;
path_resolution = 1;
%obstacle = [x y radius]
obs = [0 60 20
       100 0 30
       30 25 15
       62 95 11
       62 85 11
       62 75 11
       62 65 11
       62 55 11]';

l1 = [0 0;115 40;115 95]';
l2 = [180 0;160 55;145 55]';
l3 = [160 85;145 85;158 160]';
l4 = [40 125;90 125;115 140]';
end