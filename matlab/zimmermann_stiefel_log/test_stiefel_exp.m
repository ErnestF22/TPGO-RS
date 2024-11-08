% From https://epubs.siam.org/doi/epdf/10.1137/16M1074485
% Section 2

% Every tangent vector ∆ ∈ TU St(n, p) may be written as
% ∆ = U A + (I − U U T )T, A ∈ Rp×p skew, T ∈ Rn×p arbitrary


% St(n,p)
n = 4;
p = 3;

sz_hopefully_6 = n*p - 0.5 *p*(p + 1);

M = stiefelfactory(n,p);

%% Stiefel_Log_supp

A = make_rand_stiefel_3d_array(n,p,1);
A_tang = M.proj(A, make_rand_stiefel_3d_array(n,p,1));
tau = 1e-5;

isTg = checkIsInStiefelTangentSpace(A, A_tang);
disp('isTg')
disp(isTg)

[Delta, k, conv_hist, norm_logV0] = Stiefel_Log_supp(A, A_tang, tau);
disp('Delta');
disp(Delta);
disp('k');
disp(k);
disp('conv_hist');
disp(conv_hist);
disp('norm_logV0');
disp(norm_logV0);

%% Stiefel_exp_supp

A_exp = Stiefel_Exp_supp(A, A_tang);
disp('A')
disp(A)
disp('A_tang')
disp(A_tang)
disp('A_exp')
disp(A_exp)


%%
% From paper appendix

%
function [Delta, k, conv_hist, norm_logV0] = ...
Stiefel_Log_supp(U0, U1, tau)
%-------------------------------------------------------------
%@author: Ralf Zimmermann, IMADA, SDU Odense
%
% Input arguments
% U0, U1 : points on St(n,p)
% tau : convergence threshold
% Output arguments
% Delta : Log^{St}_U0(U1),
% i.e., tangent vector such that Exp^St_U0(Delta) = U1
% k : iteration count upon convergence
% supplementary output
% conv_hist : convergence history
% norm_logV0 : norm of matrix log of first iterate V0
%-------------------------------------------------------------
% get dimensions
[~,p] = size(U0); % !! number of rows of U0 is not used
% store convergence history
% conv_hist = 0; % was originally = [0]
% step 1
M = U0'*U1;
% step 2
[Q,N] = qr(U1 - U0*M,0); % thin qr of normal component of U1
% step 3
[V, ~] = qr([M;N]); % orthogonal completion
% "Procrustes preprocessing"
[D,~,R] = svd(V(p+1:2*p,p+1:2*p));
V(:,p+1:2*p) = V(:,p+1:2*p)*(R*D');
V = [[M;N], V(:,p+1:2*p)]; % |M X0|
% now, V = |N Y0|
% just for the record
norm_logV0 = norm(logm(V),2);
num_iter_max = 10000;
conv_hist = zeros(num_iter_max,1); % preallocating to disable warnings 
% step 4: FOR-Loop
for k = 1:num_iter_max
    % step 5
    [LV, ~] = logm(V); % ~ was originally exitflag
    % standard matrix logarithm
    % |Ak -Bk'|
    % now, LV = |Bk Ck |
    C = LV(p+1:2*p, p+1:2*p); % lower (pxp)-diagonal block
    % steps 6 - 8: convergence check
    normC = norm(C, 2);
    conv_hist(k) = normC;
    if normC<tau
        disp(['Stiefel log converged after ', num2str(k),...
        ' iterations.']);
        break;
    end
    % step 9
    Phi = expm(-C); % standard matrix exponential
    % step 10
    V(:,p+1:2*p) = V(:,p+1:2*p)*Phi; % update last p columns
end
% prepare output |A -B'|
% upon convergence, we have logm(V) = |B 0 | = LV
% A = LV(1:p,1:p); B = LV(p+1:2*p, 1:p)
% Delta = U0*A+Q*B
Delta = U0*LV(1:p,1:p) + Q*LV(p+1:2*p, 1:p);
return;

% Note: The performance of this method may be enhanced by computing expm,
% logm via a Schur decomposition.

end


function bool_out = checkIsInStiefelTangentSpace(U, Delta)
% CHECKISINTANGENTSPACE Checks whether Delta is in the tangent space
% to Stiefel manifold element U

    val = U' * Delta + Delta' * U;
    if max(abs(val), [], "all") > 1e-6
        bool_out = false;
    else
        bool_out = true;
    end
end