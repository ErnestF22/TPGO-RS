
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

clear;
close all;
clc;

resetRands(2);


d = 3;
nrs = d+1;
N = 5;

%graph random init
e = 7;
G = graph(true(N), 'omitselfloops'); % Alternative without self-loops
p = randperm(numedges(G), e);
G = graph(G.Edges(p, :));
edges = table2array(G.Edges);
num_edges = e;

%testdata random init
% testdata = testNetwork_som(3); %4 is the default option
% edges = (testdata.E);
% num_edges = testdata.NEdges;

R_globalframe = make_rand_stiefel_3d_array(nrs, d, N);
T_globalframe = 10 * rand(nrs, N);
% Tijs_vec = 10 * rand(d, num_edges);
Tijs_vec = zeros(d, num_edges);

R_transp = matStack(multitransp(R_globalframe));

P = make_p(T_globalframe, Tijs_vec, edges);
frct = compute_step1_fixed_cost_term(T_globalframe, Tijs_vec, edges);

cost_rot = trace(R_transp*P) + frct;
disp("cost_rot");
disp(cost_rot);

% cost_rot_noloops = trace(R_transp*P_noloops) + frct_noloops;
% disp("cost_rot_noloops");
% disp(cost_rot_noloops);

[A,b] = make_Aconst_b(R_globalframe, Tijs_vec, edges);


cost_1 = 0.0;
for e = 1:num_edges
    ii = edges(e,1);
    jj = edges(e,2); 
    R_i = R_globalframe(:,:,ii);
    T_j = T_globalframe(:, jj);
    T_i = T_globalframe(:, ii);
    cost_e = norm(R_i * Tijs_vec(:,e) - T_j + T_i);
    cost_1 = cost_1 + cost_e^2; %squared!
end

cost_transl = norm(vec(A .* T_globalframe) + b)^2;

disp("[cost_1, cost_rot, cost_transl]")
disp([cost_1, cost_rot, cost_transl])



