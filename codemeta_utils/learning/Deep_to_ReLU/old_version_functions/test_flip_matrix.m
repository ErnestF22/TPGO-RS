close all
clear
addpath('danyang/')
addpath('old_version_functions/')
% this function test simple network, tranform multiple times

t1 = -30;
t2 = 45;

A1 = [cosd(t1) -sind(t1);sind(t1) cosd(t1)];
b1 = [1;1];

A2 = [cosd(t2) -sind(t2);sind(t2) cosd(t2)];
b2 = [1;1];

A3 = [1 1];
b3 = 1;

% initialize
L = 3; % number of layers
n = 0;
for i=1:L
    A = eval(['A' num2str(i)]);
    n = n+size(A,1);
end
% generate A_out and b_out for z, y=A_out*x+b_out
Ab_set = struct; % Ab_set store A's and b's at each layer
for i=1:L
    Ab_set(i).A = eval(['A' num2str(i)]);
    Ab_set(i).b = eval(['b' num2str(i)]);
end

% test the matrix after flip
x = [1;1];
% zset=[1,1,1,1,1; 1,1,0,1,1; 1,0,1,1,1; 0,1,0,1,1; 0,1,1,1,1; 0,0,1,1,1];
% zset=[1,1,1,1,1; 1,1,0,1,1; 0,1,0,1,1; 0,1,1,1,1; 0,0,1,1,1; 1,0,1,1,1]';

z = [1,1,1,1,1]';

[As,basic] = getdualAmatrix(z,Ab_set);
[basic,Ar,P] = dual_simplex_xpn(As,basic);
if ~isempty(basic)
    [Q,ACT] = find_vertices_P(basic,Ar,P,z,x);
end

i = 1;
j = 1;
% C = struct;
while ~isempty(Q) % if there are unchecked points
    [Q,ACT] = find_neighbor_vertices(Q,1,ACT,x,Ab_set); % find the neighboring vertices of this vertex
    if i==1
        C = Q; % put it to C
        C(2:end) = [];
    else
        C(i) = Q(1); % put it to C
    end
    Q(1)=[]; % pop out from Q
    i = i+1;
end