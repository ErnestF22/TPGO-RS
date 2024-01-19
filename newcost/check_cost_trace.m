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



cost_1 = 0.0;
for e = 1:num_edges
    ii = edges(e,1);
    jj = edges(e,2); 
    R_i = R_globalframe(:,:,ii);
    T_j = T_globalframe(:, jj);
    T_i = T_globalframe(:, ii);
    cost_e = norm(R_i * Tijs_vec(:,e) - T_j + T_i);
    cost_1 = cost_1 + cost_e^2;
end


disp("cost_1");
disp(cost_1);


cost_2 = 0.0;
for e = 1:num_edges
    ii = edges(e,1);
    jj = edges(e,2); 
    R_i = R_globalframe(:,:,ii);
    T_j = T_globalframe(:, jj);
    T_i = T_globalframe(:, ii);
    T_ij = Tijs_vec(:,e);
    P_ij = 2.*T_i * T_ij' - 2.*T_j * T_ij';
    c = T_i * T_i' + T_j * T_j' - T_i * T_j' - T_j * T_i';
    d = T_ij * T_ij';
    cost_2 = cost_2 + trace(R_i' * P_ij + c + d);
end


disp("cost_2");
disp(cost_2);

