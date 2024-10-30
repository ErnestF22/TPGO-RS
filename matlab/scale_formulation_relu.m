function scale_formulation_relu

d = 3;
nrs = 4;
N = 5;

%graph random init
num_edges = 8;
G = graph(true(N), 'omitselfloops'); % Alternative without self-loops
p = randperm(numedges(G), num_edges);
G = graph(G.Edges(p, :));
edges = table2array(G.Edges);


lambdas = 2 * ones(num_edges, 1);
Tijs_vec = 10 * rand(d, num_edges);

R_globalframe = make_rand_stiefel_3d_array(nrs, d, N);
T_globalframe = 10 * rand(nrs, N);


rho = 10 * rand(1,1);

cost_lambda = 0.0;
for ee = 1:num_edges
    ii = edges(ee, 1);
    jj = edges(ee, 2);
    lambda_e = lambdas(ee);
    tij_e = Tijs_vec(:, ee);
    T_i = T_globalframe(:, ii);
    T_j = T_globalframe(:, jj);
    R_i = R_globalframe(:, :, ii);
    a = T_i - T_j;
    b = R_i * tij_e;
    cost_lambda_ee = trace(a' * a + 2 * lambda_e * (a' * b) + lambda_e^2 * (b' * b)); 
    cost_relu_ee = relu_som(lambda_e - 1);
    cost_lambda = cost_lambda + cost_lambda_ee + rho * cost_relu_ee;
end



cost_base = 0.0;
for ee = 1:num_edges
    ii = edges(ee,1);
    jj = edges(ee,2); 
    R_i = R_globalframe(:,:,ii);
    T_j = T_globalframe(:, jj);
    T_i = T_globalframe(:, ii);
    lambda_e = lambdas(ee);
    tij_e = Tijs_vec(:, ee);
    cost_e = norm(R_i * lambda_e * tij_e - T_j + T_i);
    cost_relu_ee = relu_som(lambda_e - 1);
    cost_base = cost_base + cost_e^2 + cost_relu_ee; %squared!
end

disp("cost_base")
disp(cost_base)

disp("cost_lambda")
disp(cost_lambda)




end %file function
