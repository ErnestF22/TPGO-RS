%this function generates inputs for paepr simulation
function [obs,start,goal,max_itr,lim_x,lim_y,expand_dis,path_resolution,l1,l2,l3,landmarks] = iros22_paper_lab_data
start = [43;89];
goal = [475 ;352];
max_itr = 400;
lim_x = [43; 475];
lim_y = [98; 590];
expand_dis = 100;
path_resolution = 1;

landmarks = table2array(readtable('landmark_coordinates_transformed.csv'));
landmarks = [landmarks(:,1) (100*landmarks(:,2:3))+[100 100]];
landmarks = landmarks';
landmarks (2:3,:) = landmarks (2:3,:)+10*rand(2,32);
obstacles = table2array(readtable('obstacle_LABexperiment.csv')); 
obstacles = (100*obstacles(:,1:2))'+[100;100];


%obstacle = [x y radius]
obs = [275 475 25;
       275 450 25;
       275 425 25;
       275 400 25;
       275 375 25;
       275 350 25;
       132 200 63;
       350 12 92;
       -57 320 59]';
obs = obs+[100;100;24]; %24 is the offset for the jackel

l1=[];
l2 = [];
l3 =[];
l4 = [];% 
figure(2)
plot(landmarks(2,:),landmarks(3,:),'m*')
hold on 
plot(obstacles(1,:),obstacles(2,:),'g*')
hold on

for i=1:size(obs,2)
    plot_circle(obs(1,i),obs(2,i),obs(3,i),'b')
    hold on
end
end