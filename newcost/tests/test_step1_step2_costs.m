
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

clear;
close all;
clc;

resetRands(2);

d = 3;
nrs = 8;
N = 80;

%graph random init
e = 29;
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
Tijs_vec = 10 * rand(d, num_edges);

R_transp = matStack(multitransp(R_globalframe));

P = make_p(T_globalframe, Tijs_vec, edges);
frct = compute_step1_fixed_cost_term(T_globalframe, Tijs_vec, edges);

cost_rot = trace(R_transp*P) + frct;
% disp("cost_rot");
% disp(cost_rot);

P_noloops = make_p_noloops(T_globalframe, Tijs_vec, edges);
frct_noloops = compute_step1_fct_noloops(T_globalframe, Tijs_vec, edges);

cost_rot_noloops = trace(R_transp*P_noloops) + frct_noloops;
% disp("cost_rot_noloops");
% disp(cost_rot_noloops);


cost_base = 0.0;
for e = 1:num_edges
    ii = edges(e,1);
    jj = edges(e,2); 
    R_i = R_globalframe(:,:,ii);
    T_j = T_globalframe(:, jj);
    T_i = T_globalframe(:, ii);
    cost_e = norm(R_i * Tijs_vec(:,e) - T_j + T_i);
    cost_base = cost_base + cost_e^2; %squared!
end

[LR,PR,BR] = make_LR_PR_BR(R_globalframe, Tijs_vec, edges);
cost_prb = trace(T_globalframe * LR * T_globalframe') + ...
    trace(T_globalframe * PR) + trace(BR);

[LR_noloops,PR_noloops,BR_noloops] = make_LR_PR_BR_noloops( ...
    R_globalframe, Tijs_vec, edges);
cost_prb_noloops = trace(T_globalframe * LR_noloops * T_globalframe') + ...
    trace(T_globalframe * PR_noloops) + trace(BR_noloops);



disp("[cost_1, " + ...
    "cost_rot, cost_rot_noloops, " + ...
    "cost_prb, cost_prb_noloops]")
disp([cost_base, ...
    cost_rot, cost_rot_noloops, ...
    cost_prb, cost_prb_noloops])




