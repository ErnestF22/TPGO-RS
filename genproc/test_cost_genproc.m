close all
clear
clc

resetRands(2);

d = 3;
nrs = 8;
N = 80;

% %graph random init
% e = 29;
% G = graph(true(N), 'omitselfloops'); % Alternative without self-loops
% p = randperm(numedges(G), e);
% G = graph(G.Edges(p, :));
% edges = table2array(G.Edges);
% num_edges = e;

%testdata random init
testdata = testNetwork_som(3); %4 is the default option
edges = (testdata.E);
num_edges = testdata.NEdges;

R_globalframe = make_rand_stiefel_3d_array(nrs, d, N);
T_globalframe = 10 * rand(nrs, N);
Tijs_vec = 10 * rand(d, num_edges);

%...




