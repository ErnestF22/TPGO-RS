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
% zset=[1,1,1,1,1; 1,1,0,1,1; 0,1,0,1,1; 0,1,1,1,1; 0,0,1,1,1; 1,0,1,1,1];

z = [1,1,1,1,1]';
[As,basic] = getdualAmatrix(z,Ab_set);
[basic,Ar,P] = dual_simplex_xpn(As,basic);
[Q1,~] = find_vertices_P(basic,Ar,P,z,x);

% z = [0,1,1,1,1]
zflip = 1;
znew = z;
znew(zflip) = ~znew(zflip); % z after flipping
A = Q1(1).A;
P = Q1(1).P;
basic = Q1(1).basic;
[Arf,Pnew] = find_tableau_after_flip_P(z,zflip,P,A,Ab_set);
Q2 = pivot_flip_matrix(basic,Arf,Pnew);

% z = [0,1,1,1,1]
z = znew;
zflip = 3;
znew = z;
znew(zflip) = ~znew(zflip); % z after flipping
A = Q2(2).A;
P = Q2(2).P;
basic = Q2(2).basic;
[Arf,Pnew] = find_tableau_after_flip_P(z,zflip,P,A,Ab_set);
Q3 = pivot_flip_matrix(basic,Arf,Pnew);
[Q3check,ACTset] = checkA(z,zflip,Ab_set,[1;1]); % use dual simplex to find vertices
sumall = abs(sum((Q3(1).A-Q3check(1).A),'all'));
if sumall<=1e-3
    disp('match!');
else
    disp('not match!');
end