function [nablaY_bar] = eulBundle_covar(X, Y, u, m, n)
% Compute the connection on the tangent bundle of R^n of Y wrt X
% NOTE: Assumes that X,Y are a vertical vector s.t. two tangent vectors on 
% R^n are stacked as [R^n; R^n].
% INPUTS:
%   X(p,u), Y(p,u) := vector fields on TR^n evaluated at (p,u)
%   u(p) := velocity vector function handle in R^n
%   m := an array of 3 values for the non natural metric on TR^n
%   n := scalar representing the dimension of the manifold
% OUTPUTS:
%   nablaY_bar(p) := the covarient derivative of Y wrt X on TR^n

%% Extract useful components
m1 = m(1); m2 = m(2); m3 = m(3);
Z = @(p) [p; u(p)]; % point to evaluate the vector fields
Xh = @(p) [eye(n) zeros(n)]*X(p,u(p));
Xv = @(p) [zeros(n) eye(n)]*X(p,u(p));
Yh = @(p) [eye(n) zeros(n)]*Y(p,u(p));
Yv = @(p) [zeros(n) eye(n)]*Y(p,u(p));

% Components of the connection using induced metric/covar
C1 = 1/(m1*m3-m2^2); % a common constant
Del_Xh_Yh = eul_covar(Xh, Yh);
Del_Xh_Yv = eul_covar(Xh, Yv);
% NOTE: there is no curvature on R^n

%% Define the horizontal components
Del_bar_X_Y_h = @(p) [C1*(m1*m3-m2^2)*Del_Xh_Yh(p);zeros(n,1)];

%% Define the vertical components
Del_bar_X_Y_v = @(p) [zeros(n,1);C1*(m1*m3-m2^2)*Del_Xh_Yv(p)];

%% Return the connection Dbar_{X}Y on TR^n
nablaY_bar = @(p) Del_bar_X_Y_h(p) + Del_bar_X_Y_v(p);

end

