function check_grad_hess_cpp_matlab

folder = "/home/rimlab/workspace/matlab_ws/som/matlab/2026_06_24_18_42_01/";

nrs = 3;
d = 3;
N = 5;

if isfile(folder+'X_vec.csv')
    
    X_vec = readmatrix(folder+'X_vec.csv');
    X = convertXtoRTLambdas(X_vec, nrs, d, N);
end

if isfile(folder+'X_tg_vec.csv')
    
    X_tg_vec = readmatrix(folder+'X_tg_vec.csv');
    X_tg = convertXtoRTLambdas(X_tg_vec, nrs, d, N);
end


if isfile(folder+'rho.csv')
    rho = readmatrix(folder + 'rho.csv');
end
problem_data.rho = rho;

if isfile(folder+'z.csv')
    z = readmatrix(folder + 'z.csv');
end
problem_data.z = z;

if isfile(folder+'y.csv')
    y = readmatrix(folder + 'y.csv');
end
problem_data.y = y;

if isfile(folder+'mu.csv')
    mu = readmatrix(folder + 'mu.csv');
end
problem_data.mu = mu;

if isfile(folder+'tijs.csv')
    tijs = readmatrix(folder + 'tijs.csv');
end
problem_data.tijs = tijs;

if isfile(folder+'edges.csv')
    edges = readmatrix(folder + 'edges.csv');
end
problem_data.edges = edges;

if isfile(folder+'a.csv')
    a = readmatrix(folder + 'a.csv');
end
problem_data.a = a;

%% compare

problem.cost = @(x) lsom_cost(x, problem_data);
% problem.egrad = @(x) ssom_egrad(x, problem_data);
problem.grad = @(x) lsom_rgrad(x, problem_data);
% problem.ehess = @(x, u) ssom_ehess_genproc(x, u, problem_data);
problem.hess = @(x, u) lsom_rhess_genproc(x, u, problem_data);


grad_X = problem.grad(X);
% if exist('grad_X','var')
%     grad_X_vec = vectorizeXrtlambdas(grad_X);
%     readmatrix(grad_X_vec, fullfile(outDir, 'grad_X_vec.csv'));
% end

hess_X_U = problem.hess(X, X_tg);
% if exist('hess_X_U','var')
%     hess_X_U_vec = vectorizeXrtlambdas(hess_X_U);
%     readmatrix(hess_X_U_vec, fullfile(outDir, 'hess_X_U_vec.csv'));
% end


end %file function
