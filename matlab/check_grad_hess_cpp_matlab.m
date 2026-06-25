function check_grad_hess_cpp_matlab

folder = "/home/rimlab/workspace/matlab_ws/som/matlab/2026_06_25_12_13_43/";

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

%read cpp grad
% Read C++ gradient if applicable
if isfile(folder + 'grad_X_cpp.csv')
    grad_X_cpp = readmatrix(folder + 'grad_X_cpp.csv');
end




grad_X_matlab = problem.grad(X);
% if exist('grad_X','var')
%     grad_X_vec = vectorizeXrtlambdas(grad_X);
%     readmatrix(grad_X_vec, fullfile(outDir, 'grad_X_vec.csv'));
% end
disp("[grad_X_matlab, grad_X_cpp]")
disp([vectorizeXrtlambdas(grad_X_matlab), grad_X_cpp])

disp("max(abs(grad_X_cpp - grad_X_matlab))")
disp(max(abs(grad_X_cpp - vectorizeXrtlambdas(grad_X_matlab))))

% Read C++ Hessian if applicable
if isfile(folder + 'hess_X_U_cpp.csv')
    hess_X_U_cpp = readmatrix(folder + 'hess_X_U_cpp.csv');
end
hess_X_U_matlab = problem.hess(X, X_tg);
% if exist('hess_X_U','var')
%     hess_X_U_vec = vectorizeXrtlambdas(hess_X_U);
%     readmatrix(hess_X_U_vec, fullfile(outDir, 'hess_X_U_vec.csv'));
% end
disp("[hess_X_U_matlab, hess_X_U_cpp]")
disp([vectorizeXrtlambdas(hess_X_U_matlab), hess_X_U_cpp])

disp(max(abs(hess_X_U_cpp - vectorizeXrtlambdas(hess_X_U_matlab))))


end %file function
