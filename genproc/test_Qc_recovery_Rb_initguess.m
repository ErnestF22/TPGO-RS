clc;
clear;
close all;


%% case 1
load('Qbnn_data/Qbnn_R.mat', 'R')
load('Qbnn_data/Qbnn_T.mat', 'T')
load('Qbnn_data/Qbnn_edges.mat', 'edges')
% load('Qbnn_data/Qbnn_Q.mat', 'Q')
% load('Qbnn_data/Qbnn_x_out.mat', 'x_out')
% load('Qbnn_data/Qbnn_T_diffs.mat', 'T_dsiffs')
% load('Qbnn_data/Qbnn.mat') %whole workspace
% Tijs = problem_data.Tijs; %loaded from ws

% load('Qbnn_data/x_rs.mat', 'x_rs')

load('Qbnn_data/Tijs.mat', 'Tijs')
node_degrees = [2,3,3,3,2];
nrs = size(R, 1);
d = size(R, 2);
N = size(R, 3);
sz = [nrs,d,N];
[R_real_out, ~] = ...
    Qc_recovery_Rb_initguess(sz, edges, R, T, Tijs, node_degrees);

load('Qbnn_data/testdata.mat', 'testdata')
R_gt = G2R(testdata.gitruth);

disp('Printing real rotations for nodes with degree == 2')
nodes_low_deg = node_degrees == 2;
for ii = 1:N
    if nodes_low_deg(ii)
        %frm (full recovery method)
        fprintf("ii %g\n", ii);
    %     disp(R_out(:,:,ii) * inv(R_gt(:,:,ii))); %slower
        disp(R_real_out(:,:,ii) / (R_gt(:,:,ii)));
    end
end


%% case 2
clc;
clear;
close all;


N = 6;
edges = [1 2; 1 6; 2 1; 2 3; 3 2; 3 6; 3 4; 3 5; 4 3; 4 5; 4 6;
    5 6; 5 3; 5 4; 6 3; 6 1; 6 5; 6 4];
adj_mat = make_adj_mat_from_edges(edges, N);
testdata = testNetwork_adj(3, adj_mat, 'banded', 3);
G = graph(make_adj_mat_from_edges(testdata.E,N));
figure(1000)
plot(G)
title('graph')

load('Qbnn_data/test_transl_recovery.mat','R')
load('Qbnn_data/test_transl_recovery.mat','T')
load('Qbnn_data/test_transl_recovery.mat','Tijs')

% nrs = size(R, 1);
d = size(R, 2);
% N = size(R, 3);

node_degrees = [2,2,4,3,3,4];
nodes_high_deg = node_degrees == 3 | node_degrees == 4;
% nodes_low_deg = ~nodes_high_deg;

R_stacked_high_deg = matStackH(R(:,:,nodes_high_deg));
Qa_high_deg = POCRotateToMinimizeLastEntries(R_stacked_high_deg);
R_stacked_high_deg_poc = Qa_high_deg * R_stacked_high_deg;

disp('R_stacked_high_deg_poc')
disp(R_stacked_high_deg_poc)


T_edges = make_T_edges(T, edges);

RT_stacked_high_deg = [matStackH(R(:,:,nodes_high_deg)) , T_edges];
Qa_high_deg = POCRotateToMinimizeLastEntries(RT_stacked_high_deg);
RT_stacked_high_deg_poc = Qa_high_deg * RT_stacked_high_deg;

disp('RT_stacked_high_deg_poc')
disp(RT_stacked_high_deg_poc)

nrs = size(R, 1);
sz = [nrs,d,N];
[R_real_out, T_real_out] = ...
    Qc_recovery_Rb_initguess(sz, edges, R, T, Tijs, node_degrees);

disp('Printing real rotations for nodes with degree == 2')
R_gt = G2R(testdata.gitruth);
nodes_low_deg = node_degrees == 2;
for ii = 1:N
    if nodes_low_deg(ii)
        %frm (full recovery method)
        fprintf("ii %g\n", ii);
    %     disp(R_out(:,:,ii) * inv(R_gt(:,:,ii))); %slower
        disp(R_real_out(:,:,ii) / (R_gt(:,:,ii)));
    end
end