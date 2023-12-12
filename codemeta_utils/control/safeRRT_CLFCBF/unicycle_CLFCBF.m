function [u] = unicycle_CLFCBF(x,xd)
% Compute control by solving CLF/CBF quadratic program
% INPUTS:
%   x := current state vector [3x1]
% OUTPUTS:
%   u := control vector [v;w] (linear, angular)

% Define QP parameters
TOL = 1e-6;
U_Bounds = [500;500];
[CLF,Grad_CLF] = compute_CLF(x,xd);
% Define objective as 1/2u'u
H = eye(2);
f = zeros(1,2);
A_clf = [Grad_CLF(1,1:2)*[cos(x(3));sin(x(3))] Grad_CLF(1,3)];

% Compose QP matrices for solver
A = [A_clf];
b = [-1];

% Solve QP
% Compute control if CLF is not zero
if (CLF > TOL)
    options = optimoptions('quadprog',...
        'Algorithm','interior-point-convex','Display','off','MaxIterations',50);
    [u,fval,exitflag,output,lambda] = ...
       quadprog(H,f,A, b);
    % ...
    %    [],[],-U_Bounds,U_Bounds,[],options);

    if (isempty(u))
        error('No feasible u found');
    end
else
    u = zeros(2,1); % No control require, [x1;x2]-[x1s;x2s]=0
end
    
end

