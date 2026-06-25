function generate_grad_hess_cpp_matlab
timestamp = datetime('now','TimeZone','local','Format','yyyy_MM_dd_HH_mm_ss');
    
disp(['Function executed at: ', char(timestamp)]);
   

nrs = 3;
d = 3;
N = 5;

mindeg = 3;

testdata = testNetwork_params(3, N, 'banded', mindeg); 

problem_data.a = 2.0;

num_edges = size(testdata.E, 1);
tijs = G2T(testdata.gij);



tuple.R = stiefelfactory(nrs, d, N);
tuple.T = euclideanfactory(nrs, N);
tuple.lambda = euclideanfactory(num_edges, 1);
M = productmanifold(tuple);

% Setup the problem structure with manifold M and cost+grad functions.
problem.M = M;


problem_data.edges = testdata.E;
problem_data.tijs = tijs;

problem_data.mu = 0.1;
problem_data.z = ones(num_edges, 1); % better to initialize this as lambdas initguess
problem_data.y = zeros(num_edges, 1);

problem_data.rho = 0.0;



% checkgradient(problem);
% tmp.R = make_rand_stiefel_3d_array(nrs, d, N);
% tmp.R = eye3d(nrs, d, N);
% tmp.T = rand(nrs, N);
% tmp.lambda = rand(num_edges, 1);
% tmpU.R = zeros(nrs, d, N);
% tmpU.T = rand(nrs, N);
% tmpU.T = normalize(tmpU.T);
% tmpU.lambda = rand(num_edges, 1);
% tmpU.lambda = normalize(tmpU.lambda);

% checkgradient(problem, tmp);
X = M.rand();
X.lambda = 2 * problem_data.a + rand(num_edges, 1); %!! function not defined on all lambdas
X_tg = M.randvec(X);
% X_tg.lambda = 2 * problem_data.a + rand(num_edges, 1);
% checkgradient(problem, X,X_tg)

% checkhessian(problem, tmp);
% checkhessian(problem, X, X_tg)

timestamp_str = string(timestamp);
outDir = fullfile(pwd, timestamp_str);         % change pwd to a parent path if needed
[status, msg, msgID] = mkdir(outDir);
if ~status
    error('Failed to create output folder "%s": %s (%s)', outDir, msg, msgID);
end

if exist('X','var')
    X_vec = vectorizeXrtlambdas(X);
    writematrix(X_vec, fullfile(outDir, 'X_vec.csv'));
end

if exist('X_tg','var')
    X_tg_vec = vectorizeXrtlambdas(X_tg);
    writematrix(X_tg_vec, fullfile(outDir, 'X_tg_vec.csv'));
end

rho = problem_data.rho;
if exist('X_tg','var')
    writematrix(rho, fullfile(outDir, 'rho.csv'));
end

z = problem_data.z;
if exist('z','var')
    writematrix(z, fullfile(outDir, 'z.csv'));
end

y = problem_data.y;
if exist('y','var')
    writematrix(y, fullfile(outDir, 'y.csv'));
end

mu = problem_data.mu;
if exist('mu','var')
    writematrix(mu, fullfile(outDir, 'mu.csv'));
end

num_edges = size(testdata.E, 1);
if exist('num_edges','var')
    writematrix(num_edges, fullfile(outDir, 'num_edges.csv'));
end

% tijs = problem_data.tijs;
if exist('tijs','var')
   writematrix(tijs, ...
        convertStringsToChars(strcat(outDir, "/tijs.csv")), 'Delimiter', ',')
end

edges = problem_data.edges;
if exist('edges','var')
    writematrix(edges, fullfile(outDir, 'edges.csv'));
end

a = problem_data.a;
if exist('a','var')
    writematrix(a, fullfile(outDir, 'a.csv'));
end

end %file function