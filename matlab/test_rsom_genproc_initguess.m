clc;
clear;
close all;

nrs = 4;
d = 3;
N = 5;



num_edges = 8;
G = graph(true(N), 'omitselfloops'); % Alternative without self-loops
p = randperm(numedges(G), num_edges);
G = graph(G.Edges(p, :));
edges = table2array(G.Edges);


Tijs = 10 * rand(d, num_edges);
problem_data.Tijs = Tijs;
problem_data.edges = edges;




problem_data.Tijs = Tijs;
problem_data.edges = edges;




tuple.R = stiefelfactory(nrs, d, N);
tuple.T = euclideanfactory(nrs, N);
M = productmanifold(tuple);

% Setup the problem structure with manifold M and cost+grad functions.
problem.M = M;
problem_data.sz = [nrs, d, N];
problem.cost = @(x) cost_genproc(x, problem_data);
problem.egrad = @(x) egrad_genproc(x, problem_data);
problem.ehess = @(x, u) ehess_genproc(x, u, problem_data);


figure(1)
checkgradient(problem);

figure(2)
checkhessian(problem);

R_start = zeros(nrs, d, N);
R_start(:,:,1) = [...
-8.6111316201e-01	-2.4231807714e-01	2.2041231019e-01;
2.9905667201e-01	-3.2397439207e-01	8.9628847729e-01;
4.0690794745e-01	-3.9424233375e-01	-2.3628533360e-01;	
5.8950415130e-02	8.2516393829e-01	3.0373445658e-01];	

R_start(:,:,2) = [ ...
-4.1059440931e-01	-1.1437917202e-01	2.7183095287e-01;
3.5386645467e-01	-8.6801045842e-01	-3.1338931219e-01;
6.3411042845e-01	4.7706138888e-01	-3.7143849333e-01;	
5.5144784689e-01	-7.6731221225e-02	8.3061935791e-01];	

R_start(:,:,3) = [ ...
-2.1430690851e-01	-2.4534556283e-01	-5.3386741786e-01;
-6.2161433491e-01	-6.9622980720e-01	3.1324520402e-01;	
2.0106036221e-01	-3.6483516939e-01	-7.2147848594e-01;	
7.2611493466e-01	-5.6741951813e-01	3.1037367259e-01];	

R_start(:,:,4) = [ ...
-7.1037398226e-01	-5.0679527971e-01	2.7375993632e-01;
3.0531086862e-01	-5.0948717294e-01	-7.0995604558e-01;	
4.2847851248e-02	5.9294792062e-01	2.1311271115e-02;	
6.3270699417e-01	-3.6330996201e-01	6.4850885910e-01];	

R_start(:,:,5) = [ ...
-1.5666248790e-01	5.8716493319e-01	-1.7441611935e-01;	
-5.0524314380e-01	-7.1781280960e-01	1.0920867765e-01;	
-9.6019087979e-02	-1.5821933671e-01	-9.7547959060e-01;	
-8.4318833322e-01	3.3904093130e-01	7.8051587764e-02];	

% T_start = zeros(nrs, N);
T_start = [...
3.2774819472e-01	-1.1292397076e+00	1.3968437082e-01	3.8852371422e-02	2.9319206102e-01;	
7.6930246442e-02	9.2838890252e-02	-1.8538640042e-02	-1.9448161338e-01	1.0684563714e+00;	
2.6546361981e+00	4.6517108840e-02	3.7439802807e-01	1.4314678759e+00	1.0700675098e-01;	
-1.1277842920e+00	-6.4363522916e-02	-3.2549888398e-01	8.7314119356e-02	4.7163430369e-01];	

X_start.R = R_start;
X_start.T = T_start;

[X_out, cost_out, ~, ~] = trustregions(problem, X_start);

T_out = X_out.T;
R_out = X_out.R;

disp("R_out")
disp(R_out)

disp("T_out")
disp(T_out)

disp("cost_out")
disp(cost_out)