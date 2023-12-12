close all
clear
addpath('foundamental functions/')
addpath('nnet/')
addpath('old_version_functions')
addpath('find_plane_v2')

load('/Users/danyangli/Documents/mip/Ab_set.mat');
load('/Users/danyangli/Documents/mip/x_mip.mat');
load('/Users/danyangli/Documents/mip/y_mip.mat');
load('/Users/danyangli/Documents/mip/z_mip.mat');
load('x_input.mat');
z_mip = z_mip';

L = size(Ab_set,2);
[y,z_simplex] = forward(x_input,L,Ab_set); % get z with input
[As,basic] = getdualAmatrix(z_simplex,Ab_set);
[basic,Ar,P] = dual_simplex_xpn(As,basic);
if ~isempty(basic)
    [Q,~,ZSet] = find_plane_P(basic,As,Ar,P,z_simplex,x_input,Ab_set)
else
    disp('No vertices found');
end

find_match_flag = false;
d = size(As,2)-size(As,1);
z = Q(1).z;
while ~find_match_flag
    [Q,ZSet,z,find_match_flag] = find_match_z(Q,Ab_set,d,z_simplex,z_mip,ZSet,find_match_flag)
end