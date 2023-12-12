function [ m, flag ] = contraction3D_LMIopt( eMin, eMax, kvIn, kdIn, betaIn )
% Use the LMI as the objective and constraint. NOTE this method has not
% been shown to have max and min contraction eigenvalues at the max and min
% eigenvalues of the Hessisan. However, base on results of
% contraction3D_plotEvals this method holds.
% INPUT:
%   eMin := min eigenvalue of the Hessian of the cost function
%   eMax := max eigenvalue of the Hessian of the cost function
%   kvIn := scalar velocity gain
%   kdIn := scalar position gain
%   betaIn := scalar for convergence rate
% OUTPUT:
%   m := [1 m12 m22] values
%   flag := boolean where true means the system is feasible

%% Default Parameters
MAX_ITER = 100;
alpha = 10;
eps = 1e-3;
ERR_TOL = 1e-5;
flag = false; %set to false by default
%% Optimize using LMI as constraints
cvx_begin sdp quiet
    variable m(2,2) semidefinite
% %     variable gammaVar
    Gmat = @(eVal) [2*m(1,2)+betaIn*m(2,2)-2*kvIn*m(2,2), betaIn*m(1,2)-kvIn*m(1,2)-eVal*kdIn*m(2,2)+m(1,1);...
        betaIn*m(1,2)-kvIn*m(1,2)-eVal*kdIn*m(2,2)+m(1,1), betaIn*m(1,1)-2*eVal*kdIn*m(1,2)];
    Gmax = Gmat(eMax);
    Gmax <= 0
    Gmin = Gmat(eMin);
    Gmin <= 0
    m(1,1) == 1
cvx_end

if  strcmp(lower(cvx_status),'solved')
    % double check to make sure the G matrices are negative semidefinite
    eGmax = eig(Gmax);
    eGmin = eig(Gmin);
    em = eig(m);
    if ( all(eGmin <= 0) && all(em > 0) && all(eGmax <= 0))
%         fprintf('System is feasible\n');
        flag = true;
    end
end

end

