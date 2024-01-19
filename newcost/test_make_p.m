clear;
close all;
clc;

resetRands(2);


d = 3;
% nrs = d;
% % nrs_next = nrs + 1;
% N = 5;


testdata = testNetwork_som(3); %4 is the default option

N = testdata.NNodes;
edges = (testdata.E);
num_edges = testdata.NEdges;

% Tijs_vec = G2T(testdata.gijtruth);
% R_globalframe = G2R(testdata.gitruth);
% T_globalframe = G2T(testdata.gitruth);

R_globalframe = randrot(d, N);
T_globalframe = 10 * rand(d, N);
Tijs_vec = 10 * rand(d, num_edges);

%make P matrix with and without loops and check that the result is the same

P = make_p(R_globalframe, T_globalframe, Tijs_vec, edges);
P_noloops = make_p_noloops(R_globalframe, T_globalframe, Tijs_vec, edges);


% fixed rot cost term
frct = compute_step1_fixed_cost_term(T_globalframe, Tijs_vec, edges);
R_transp = matStack(multitransp(R_globalframe));


cost_2 = trace(R_transp*P) + frct;

disp("cost_2");
disp(cost_2);



