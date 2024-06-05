clc;
clear;
close all;

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

nrs = 4;
d = 3;
N = 5;
sz = [nrs,d,N];
Qc_recovery_Rb_initguess(sz, edges, R, T, Tijs);



