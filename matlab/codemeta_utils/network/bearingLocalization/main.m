clear all
close all
clc

addpath util

%% parameters: adjacency matrix + noise 
 
% if Adj(i,j) =1, then i sees j, that is t_ij^(i) is available


Adj = [0 1 1 0 0 0 0 0 0 1 1;...
       1 0 1 1 0 0 0 0 0 0 1;...
       1 1 0 1 1 0 0 0 0 0 0;...
       0 1 1 0 1 1 0 0 0 0 0;...
       0 0 1 1 0 1 1 0 0 0 0;...
       0 0 0 1 1 0 1 1 0 0 0;...
       0 0 0 0 1 1 0 1 1 0 0;...
       0 0 0 0 0 1 1 0 1 1 0;...
       0 0 0 0 0 0 1 1 0 1 1;...
       1 0 0 0 0 0 0 1 1 0 1;...
       1 1 0 0 0 0 0 0 1 1 0];


%Adj = ones(25);
%Adj(1:size(Adj,1)+1:end) = 0;
   
sigma_in_deg = 2; % in degrees

 


%% generate the network and the measurements


E = E_from_Adj(Adj); % edge list
F = get_higher_order_edges(E,Adj); % higher order edges 
[Rgt,Tgt] = generate_poses(size(Adj,1)); % groundtruth poses
[tijgt,tijigt,tiji,tiji_mtrx] = generate_baerings(Rgt,Tgt,E,size(Adj,1),sigma_in_deg*pi/180); % true and noisy bearings
M = get_marix_M(E,F,size(Adj,1),tiji,tiji_mtrx); % create matrix M

%% optimization 

[X] = solve_manopt(Adj,F,M); %  manopt

[Rb,cost] = EstimateRotationsRGD(Adj,F,M, 1000); % maxIter = 200; %  distributed gradient descent

%% compare objectives
cost_gt   =  evaluate_cost(Adj,F,M,Rgt);
cost_mnpt =  evaluate_cost(Adj,F,M,X);
cost_dgd  =  evaluate_cost(Adj,F,M,Rb);

display('======================================================')
display(sprintf('Cost for groundtruth solution: % 5.2f',cost_gt))
display(sprintf('Cost for centralized solution: % 5.2f',cost_mnpt))
display(sprintf('Cost for distributed solution: % 5.2f',cost_dgd))

 

%% rotation errors gradient descent vs groundtruth

[Yhat_al_3d,err_deg,err_rad] = align_Rotations_Procrustes(Rgt,Rb);

display('======================================================')
display('Errors in degrees')
err_deg
mean(err_deg)
std(err_deg)
 

 


