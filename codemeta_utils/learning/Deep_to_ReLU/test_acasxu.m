close all
clear
addpath('foundamental functions/')
addpath('nnet/')
addpath('find_point')
addpath('find_plane_v2')

% input [rho, theta, phi, v_own, v_int]
% output [COC, WL, WR, SL, SR]

%1-3 1-5
i = 1;
j = 6;

% the score for COC is at most 1500
A_add = [1 0 0 0 0];
b_add = -1500;
% ρ ≥ 55947.691, vown ≥ 1145, vint ≤ 60.
x_input = [55947.691;0;0;1145;60];
z_desired = 0;

% % the score for COC is the minimal
% A_add = [1 -1 0 0 0; % COC-WL<=0 => COC<=WL
%          1 0 -1 0 0;
%          1 0 0 -1 0;
%          1 0 0 0 -1];
% b_add = [0 0 0 0]';
% % 1500 ≤ ρ ≤ 1800, −0.06 ≤ θ ≤ 0.06, ψ ≥ 3.10, vown ≥ 980,vint ≥ 960.
% x = [1800;0.06;3.1;980;960];
% % x = [0;0;0;0;0];
% % x = [1700;0;5;1000;1000];
% z_desired = [0;0;0;0];
network = ['nnet/ACASXU_run2a_' num2str(i) '_' num2str(j) '_batch_2000.nnet'];
Ab_set = load_network(network);
mip = 1;
if mip
    Ab_set(3:end-1) = [];
end
Ab_set(end+1).A = A_add;
Ab_set(end).b = b_add;
save(['/Users/danyangli/Documents/mip/' 'Ab_set']);
checkpoint = false;
if checkpoint
    [Q,violate_loc] = find_closest_point(Ab_set,z_desired,x_input);
    z_simplex = Q(violate_loc).z;
    v_simplex = Q(violate_loc).v;
else
    [Q,violate_loc] = find_closest_plane(Ab_set,z_desired,x_input);
    z_simplex = Q(violate_loc).z;
    disp(Q(violate_loc).h_plane);
end
save('z_simplex.mat','z_simplex');
save('x_input.mat','x_input');
% save('v_simplex.mat','v_simplex');