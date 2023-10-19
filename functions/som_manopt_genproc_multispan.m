function transf_curr = som_manopt_genproc_multispan(T_globalframe, Tijs_vec, edges, params, transf_init_guess)
%SOM_GENPROC Computes transformation through Procrustes pipeline in only
%one step (rotation and translation are not separated).
%The genproc name comes from the Generalized Procrustes class that is used
%by the Manopt called procedure (opposed to the stepone, steptwo procedure
%tested initially)
%Inputs can obviously be noisy

d = params.d;
N = params.N;
% transf_end_thresh = params.transf_end_thresh;
% max_icp_iterations = params.max_icp_iterations;
% initguess_is_available = params.initguess_is_available;

Tijs_mat = tijs_vec_2_tijs_mat(Tijs_vec, edges, N);

disp("NOTE: At least one ICP iteration is done by default");


% Construct a manifold structure representing the product of groups of
% rotations with the Euclidean space for A. We optimize simultaneously
% for the reference cloud and for the rotations that affect each of the
% measured clouds. Notice that there is a group invariance because
% there is no way of telling which orientation the reference cloud
% should be in.
tuple.R = rotationsfactory(d, N);
tuple.A = euclideanfactory(d, N);
M = productmanifold(tuple);

% Define the cost function here. Points on the manifold M are
% structures with fields X.A and X.R, containing matrices of sizes
% respectively nxm and nxnxN. The store structure (the caching system)
% is used to keep the residue matrix E in memory, as it is also used in
% the computation of the gradient and of the Hessian. This way, we
% prevent redundant computations.
function [f, store] = cost(X, store)
    if ~isfield(store, 'E')
        R = X.R;
        A = X.A;
        store.E = multiprod(R, A);
    end
    %L,P,const_cost_term_tij
    [L, P] = make_LT_PT_noloops(T_globalframe, Tijs_vec, edges, params);
    cost_const_term_tij = compute_fixed_cost_term(Tijs_vec, d); 
    f.R = trace((matStack(R))'*L*(matStack(R)) + (matStack(R))'*P) + cost_const_term_tij;
    %A,b
    %Note: variable name A is already used for the translations!
    [Amat,b] = make_A_b_noloops(R, T_globalframe, Tijs_vec, Tijs_mat, edges, params);
    f.A = (Amat*A(:) + b)'*(Amat*A(:) + b);
end

% Riemannian gradient of the cost function.
function [g, store] = grad(X, store)
    R = X.R;
    A = X.A;
    if ~isfield(store, 'E')
        [~, store] = cost(X, store);
    end
    E = store.E;
    % Compute the Euclidean gradient of the cost wrt the rotations R
    % and wrt the cloud A,
    %make L,P
    [L, P] = make_LT_PT_noloops(T_globalframe, Tijs_vec, edges, params);
    egrad.R = matUnstack(L*matStack(R) + (L')*matStack(R) + P);
    %A,b
    %Note: variable name A is already used for the translations!
    [Amat,b] = make_A_b_noloops(R, T_globalframe, Tijs_vec, Tijs_mat, edges, params);
    egrad.A = 2 * (Amat' * Amat) * A(:) + 2*Amat'*b;
    % then transform this Euclidean gradient into the Riemannian
    % gradient.
    g = M.egrad2rgrad(X, egrad);
    store.egrad = egrad;
end


% Setup the problem structure with manifold M and cost+grad functions.
problem.M = M;
problem.cost = @cost;
problem.grad = @grad;
% problem.hess = @hess;

% An alternative way to compute the gradient and the hessian is to use 
% automatic differentiation provided in the deep learning toolbox (slower)
% problem.cost = @cost_AD;
%    function f = cost_AD(X)
%        R = X.R;
%        A = X.A;
%        E = multiprod(R, A) - A_measure;
%        f = (E(:)'*E(:))/(2*N);
%    end
% call manoptAD to prepare AD for the problem structure
% problem = manoptAD(problem);

% For debugging, it's always nice to check the gradient a few times.
% checkgradient(problem);
% pause;
% checkhessian(problem);
% pause;

% Call a solver on our problem. This can probably be much improved if a
% clever initial guess is used instead of a random one.
options.verbosity = 0;
X = trustregions(problem, [], options); %[] is in place of initguess
A = X.A;
R = X.R;

% To evaluate the performance of the algorithm, see how well Atrue (the
% reference cloud) matches A (the found cloud). Since the recovery is
% up to rotation, apply Kabsch algorithm (or standard Procrustes),
% i.e., compute the polar factorization to best align Atrue and A.
if exist('Atrue', 'var')
    [U, ~, V] = svd(Atrue*A');
    Ahat = (U*V')*A;
    fprintf('Registration error: %g.\n', norm(Atrue-Ahat, 'fro'));
end

% Plot the spectrum of the Hessian at the solution found.
% Notice that the invariance of f under a rotation yields dim SO(n),
% that is, n*(n-1)/2 zero eigenvalues in the Hessian spectrum at the
% solution. This indicates that critical points are not isolated and
% can theoretically prevent quadratic convergence. One solution to
% circumvent this would be to fix one rotation arbitrarily. Another
% solution would be to work on a quotient manifold. Both can be
% achieved in Manopt: they simply require a little more work on the
% manifold description side.
if M.dim() <= 512
    stairs(sort(hessianspectrum(problem, X)));
    title('Spectrum of the Hessian at the solution found.');
    xlabel('Eigenvalue number (sorted)');
    ylabel('Value of the eigenvalue');
end

end %file function







% It is not necessary to define the Hessian of the cost. We do it
% mostly to illustrate how to do it and to study the spectrum of the
% Hessian at the solution (see further down).
% % function [h, store] = hess(X, Xdot, store)
% %     R = X.R;
% %     A = X.A;
% %     % Careful: tangent vectors on the rotation group are represented as
% %     % skew symmetric matrices. To obtain the corresponding vectors in
% %     % the ambient space, we need a little transformation. This
% %     % transformation is typically not needed when we compute the
% %     % formulas for the gradient and the Hessian directly in Riemannian
% %     % form instead of resorting the egrad2rgrad and ehess2rhess. These
% %     % latter tools are convenient for prototyping but are not always
% %     % the most efficient form to execute the computations.
% %     Rdot = tuple.R.tangent2ambient(R, Xdot.R);
% %     Adot = Xdot.A;
% %     if ~isfield(store, 'egrad')
% %         [~, store] = grad(X, store);
% %     end
% %     E = store.E;
% %     egrad = store.egrad;
% %     
% %     ehess.R = multiprod(multiprod(Rdot, A) + multiprod(R, Adot), A') + ...
% %               multiprod(E, Adot');
% %     ehess.R = ehess.R / N;
% %     ehess.A = Adot-mean(multiprod(multitransp(Rdot), A_measure), 3);
% %     
% %     h = M.ehess2rhess(X, egrad, ehess, Xdot);
% % end

