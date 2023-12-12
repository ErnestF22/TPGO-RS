close all
clear
addpath('foundamental functions/')
addpath('find_plane_v2')

t1 = -30;
t2 = 45;

A1 = [cosd(t1) -sind(t1);sind(t1) cosd(t1)];
b1 = [1;1];

A2 = [cosd(t2) -sind(t2);sind(t2) cosd(t2)];
b2 = [1;1];

A3 = [-10 1];
b3 = 0.5;

nbColors=32;
colors=parula(nbColors);

% number of layers
L = 3;

% generate all possible z
n = 0;
for i=1:L
    A = eval(['A' num2str(i)]);
    n = n+size(A,1);
end
z_set=dec2bin(0:2^n-1)-'0';

for i=1:L
    Ab_set(i).A = eval(['A' num2str(i)]);
    Ab_set(i).b = eval(['b' num2str(i)]);
end

% Ab_set = {};
% Ab_set(1).A = [];
% Ab_set(1).b = [];
% L = 4;
% for i =1:L
%     d = 2;
%     limL = -15; limU = 15;
%     A1 = limL + (limU-limL) * rand([d,d]);
%     b1 = limU*rand(2,1);
%     Ab_set(end+1).A = A1;
%     Ab_set(end).b = b1;
% end
% Ab_set(1) = [];

x_input = [0;0];
% x = [-0.912426098337680;-0.419631639527302];
[~,z] = forward(x_input,L,Ab_set); % get z with input
save(['/Users/danyangli/Documents/mip/' 'Ab_set']);
n = size(z,1);
z_set=dec2bin(0:2^n-1)-'0';

zset = [];
h_input_set = [];
vset = [];
ACTset = [];
Qset = {};
x_input = [0;0];
for z = z_set'
    [As,basic] = getdualAmatrix(z,Ab_set);
    [basic,Ar,P] = dual_simplex_xpn(As,basic);
    if ~isempty(basic)
%         [Q,ACT] = find_vertices_P(basic,Ar,P,z,x_input,plane);
        [Q,ACT] = find_plane_P(basic,As,Ar,P,z,x_input,Ab_set);
        Qset = [Qset;Q];
        if isempty(ACTset)
            ACTset = ACT;
            h_plane_set = extractfield(Q,'h_plane');
            vset = reshape(extractfield(Q,'point_v'),[size(Q,2),size(Q(1).v,2)]);
            continue
        end
        for i = 1:size(ACT,1)
            act = ACT(i,:);
            ACTset = [ACTset;act];
            h_plane_set = [h_plane_set;Q(i).h_plane];
            vset = [vset;Q(i).point_v];
        end
    end
end
[min_distance,min_loc] = min(h_plane_set);
v_simplex = vset(min_loc,:);
save('v_simplex.mat','v_simplex');
save('x_input.mat','x_input');