function [Ft, dFt, appder, matder] = eulBundle_metricCompatibility(p, dp, u, du, X, Y, m, n, t)
% Check if d/dt<X,Y> = <D_{Z_dot}X,Y> + <X,D_{Z_dot}Y> evaluated at a given
% t
% INPUTS:
%   p(t) := Geodesic on R^n
%   dp(p) := The vector field induced by p(t) as a function of p
%   u(p,t) := Velocity component of the point Z(p,u) st. [p,u] is a path 
%       on TR^n
%   du := The numerical time derivative of u at (p,u) and time p
%   X(p,u), Y(p,u) := Vector fields on TR^n evaluated at (p,u)
%   m := an array of 3 values for the non natural metric on TR^n
%   n := scalar representing the dimension of the manifold
%   t := time to evaluate the metric
% OUTPUTS:
%   Ft := Numerical evaluation of <X,Y>
%   dFt := Analytical evaluation of <D_{Z_dot}X,Y> + <X,D_{Z_dot}Y>
%   appder := Nuerical evaluation of d/dt<X,Y>


%% Extract useful components
Z = @(t) [p(t);u(p(t),t)];
u_p = @(p) u(p,t); % u evaluated at p and time t
Z_dot_lower = eul_covar(dp, u_p);
Z_dot = @(p,u) [dp(p); Z_dot_lower(p) + du];
Y = Z_dot; % Assume Y == Z_Dot

%% Define the LHS
metric_X_Y = @(t) eulBundle_metric(Z(t), ...
    X(p(t),u(p(t),t)), ...
    Y(p(t),u(p(t),t)), m, n);

%% Define the RHS
D_Zdot_X = eulBundle_covar(Z_dot,X,u_p,m,n);
D_Zdot_Y = eulBundle_covar(Z_dot,Y,u_p,m,n);

D_1 = @(t) eulBundle_metric(Z(t),D_Zdot_X(p(t)), Y(p(t),u_p(p(t))), m, n);
D_2 = @(t) eulBundle_metric(Z(t),D_Zdot_Y(p(t)), X(p(t),u_p(p(t))), m, n);

Right = @(t) D_1(t) + D_2(t);

%% Perform numeric test
[Ft, dFt, appder] = funCheckDer(metric_X_Y, Right, t, 'nodisplay');

%% Find result using contraction matrix result 
M = eulBundle_metric_contractionMat(5,2,m,n);
zeta = dp(p(t)); eta = Z_dot_lower(p(t)) + du;
matder = [zeta' eta']*M*[zeta;eta];
end

